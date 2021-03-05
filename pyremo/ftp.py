#! /usr/bin/python
# coding: utf-8
#

"""This module is for ftp functions.

It provides mostly functions needed to upload data to and download data
from the tape archive at dkrz. They are part of the pre- and post-processing
scripts for REMO.
"""
import logging
import subprocess
import os

__all__ = [
    'ftp_call', 'ftp_tape', 'check_ftp_path', 'make_ftp_path',
    'rename_ftp_path', 'upload_file', 'download_file', 'delete_file']


def ftp_call(server, ftp_command, ignore_errors=False):
    """Performs an ftp-command on a ftp-server.

    Args:
	server:
	    ftp-server to use the command on
        ftp_command:
            Command that is to run on the tape archive
        ignore_errors:
            Ignore ftp-errors

    Returns:
        ftp_output:
            The output message of the ftp-server

    Written by Claas Teichmann
    """
    logging.info('FTP Commands:')
    for i in ftp_command.splitlines():
        logging.info(i)
    if server == 'tape':
        ftp_output = ftp_tape(ftp_command, ignore_errors=ignore_errors)
    else:
        raise Exception
    return ftp_output


def ftp_tape(ftp_command, ignore_errors=False):
    """Performes an ftp-command on the dkrz tape archive.

    Args:
        ftp_command:
            Command that is to run on the tape archive.
        ignore_errors:
            Ignore ftp-errors.

    Returns:
        ftp_output:
            The output message of the ftp-server.

    Written by Claas Teichmann
    """
    pftp = subprocess.Popen(['pftp'], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, stdin=subprocess.PIPE,
                            shell=False)
    ftp_output = pftp.communicate(ftp_command.encode('utf-8'))
    #ftp_output[0] = ftp_output[0].decode()
    #ftp_output[1] = ftp_output[1].decode()
    logging.debug('ftp_output (stdout):')
    log_ftp_output(ftp_output)
    if 'Error' in ftp_output[0].decode() + ftp_output[1].decode() and not ignore_errors:
        raise Exception
    return ftp_output


def log_ftp_output(ftp_output, logging_function=logging.debug):
    """Passes ftp_output to the logging-function.

    Args:
        ftp_output:
            output from an ftp-command
    """
    logging_function('ftp_output (stdout)')
    for line in ftp_output[0].splitlines():
        logging_function(line)
    logging.debug('ftp_output (stderr)')
    for line in ftp_output[1].splitlines():
        logging_function(line)


def check_ftp_path(server, ftp_path):
    """Returns True if path exists.

    Logs into an ftp-server and checks whether the given path exists.

    Args:
        server:
            Name of the server
        ftp_path:
            Path on the server to check

    Returns:
        path_exists:
            True if path exits
    """
    ftp_command = ' ls {0} \n quit \n'.format(ftp_path)
    ftp_output = ftp_call(server, ftp_command)
    logging.debug('Results of checking directory on {0}'.format(server))
    log_ftp_output(ftp_output)
    if 'Could not stat' in ftp_output[0]:
        # this means the directory exits
        path_exists = False
    else:
        path_exists = True
    logging.debug('Results of checking directory on {0}'.format(server))
    return path_exists


def check_ftp_path_list(server, ftp_path_list):
    """Returns True if path exists.

    Logs into an ftp-server and checks whether the given list of pathes exists.

    Args:
        server:
            Name of the server
        ftp_path_list:
            List of Paths on the server to check

    Returns:
        path_exists:
            List of True/False if path exits
    """
    ftp_commands = ''
    for ftp_path in ftp_path_list:
       ftp_commands += ' ls {0} \n'.format(ftp_path)
    ftp_output = ftp_call(server, ftp_commands)
    logging.debug('Results of checking directory on {0}'.format(server))
    log_ftp_output(ftp_output)
    path_exists = []
    for line in ftp_output[0].splitlines()[1:]:
      if 'ftp>' not in line.decode(): # ignore lines that show the ftp command
        if 'Could not stat' in line.decode():
            # this means the directory/file does not exit
            path_exists.append((line,False))
        else:
            path_exists.append((line,True))
    logging.debug('Results of checking directory on {0}'.format(server))
    return path_exists


def make_ftp_path(server, ftp_path,
                  rename_if_path_exists=False, path_suffix='_backup'):
    """Creates a directory on the server.

    This routine checks for the existance of a directory on the
    ftp-server. If the directory exits, it raises an Exception.
    The existing diretory can also be renamed (in case rename_if_exits == True).

    Args:
        server:
            Name of the server
        ftp_path:
            Path on the server to check
        rename_if_path_exists:
            Move the directory if it exists
        path_suffix:
            Adding this suffix when renaming the directory [default: '_backup']

    Written by Claas Teichmann
    """
    path_exists = check_ftp_path(server, ftp_path)
    if path_exists:
        if rename_if_path_exists:
            new_ftp_path = ftp_path + path_suffix
            rename_ftp_path(server, ftp_path, new_ftp_path)
        else:
            raise Exception
    else:
        ftp_command = (' mkdir  {0} \n '.format(ftp_path))
    ftp_output = ftp_call(server, ftp_command)
    logging.debug('Results of making directory on {0}'.format(server))
    log_ftp_output(ftp_output)


def rename_ftp_path(server, ftp_path, new_ftp_path):
    """Rename a directory on the server.

    A directory on the ftp-server is renamed. If the new path already exists
    or the original path not, an exception is raised.

    Args:
        server:
            Name of the server.
        ftp_path:
            Path on the server to be renamed.
        new_ftp_path:
            New Path name.

    Written by Claas Teichmann
    """
    #
    # checking whether the new path and the source path exist
    #
    new_path_exists = check_ftp_path(server, new_ftp_path)
    if new_path_exists:
        raise Exception
    path_exists = check_ftp_path(server, ftp_path)
    if path_exists:
        ftp_cmd = (' rename  {0} {1} \n '.format(ftp_path,
                                                 new_ftp_path))
    else:
        raise Exception
    ftp_output = ftp_call(server, ftp_cmd)
    logging.debug('Results of renaming directory on {0}'.format(server))
    log_ftp_output(ftp_output)


def upload_file(server, sourcepath, datafile, destdir,
                remove_source=False):
    """Upload a file to an ftp-server.

    Args:
        server:
            Name of the server.
        sourcepath:
            Source path where of the file to be uploaded
        datafile:
            Name of the source file to be uploaded
        destdir:
            Destination path on the ftp-server where the file will be stored
        remove_source:
            Remove the source file after uploading

    Written by Claas Teichmann
    """
    cwd = os.getcwd()
    logging.info('changing to ' + sourcepath)
    os.chdir(sourcepath)
    ftp_command = ' cd {0} \n put {1} \n quit \n'.format(destdir,
                                                         datafile)
    ftp_output = ftp_call(server, ftp_command)
    if ftp_output[0].strip().endswith(b'quit'):
        logging.info('Upload successfull')
        if remove_source:
            os.remove(datafile)
    else:
        log_ftp_output(ftp_output, logging_function=logging.warning)
        raise Exception
    logging.info('changing back to '+ cwd)
    os.chdir(cwd)


def download_file(server, sourcepath, datafile, destdir):
    """Download a file from an ftp-server.

    If the file exists locally, the download is skipped.

    Args:
        server:
            Name of the server
        sourcepath:
            Source path of the file in the achive
        datafile:
            Name of the source file in the archive
        destdir:
            Destination path where the file will be stored locally

    Written by Claas Teichmann
    """
    cwd = os.getcwd()
    logging.info('changing to ' + destdir)
    os.chdir(destdir)
    if datafile in os.listdir('.'):
        logging.warning('File already exists: {0}'.format(datafile))
        logging.warning('Skipping download!')
    else:
        ftp_command = 'cd {0} \n get {1} \n quit \n'.format(sourcepath,
		                                            datafile)
        ftp_output = ftp_call(server, ftp_command)
        if ftp_output[0].strip().endswith(b'quit'):
            logging.info('Successfully downloaded {0}'.format(datafile))
        else:
            log_ftp_output(ftp_output, logging_function=logging.warning)
            raise Exception
    logging.info('changing back to ' + cwd)
    os.chdir(cwd)


def delete_file(server, ftp_path, datafile):
    """Deletes a file on an ftp-server.

    Args:
        server:
            Name of the server
        ftp_path:
            Path on the server to the file
        datafile:
            Name of the source file in the archive

    Written by Kevin Sieck
    """
    file_exists = check_ftp_path(server, os.path.join(ftp_path, datafile))
    if file_exists:
        ftp_command = 'cd {0} \n delete {1} \n quit \n'.format(ftp_path,
                                                               datafile)
        ftp_output = ftp_call(server, ftp_command)
        if ftp_output[0].strip().endswith(b'quit'):
            logging.info('Successfully deleted {0}'.format(datafile))
        else:
            log_ftp_output(ftp_output, logging_function=logging.warning)
            raise Exception
    else:
        raise Exception


def ftp_chmod(server, ftp_path, chmod_arg, datafile=None):
    """Changes access rights of a file/directory on an ftp-server.

    Args:
        server:
            Name of the server
        ftp_path:
            Path on the server to a directory
        chmod_arg:
            argument for chmod.
        datafile:
            Optional name of a file in the archive

    Written by Kevin Sieck
    """
    if datafile == None:
        datafile = '.'
        path_exists = check_ftp_path(server, ftp_path)
    else:
        path_exists = check_ftp_path(server, os.path.join(ftp_path, datafile))

    if path_exists:
        ftp_command = 'cd {0} \n chmod {1} {2} \n quit \n'.format(ftp_path,
                                                                  chmod_arg,
                                                                  datafile)
        ftp_output = ftp_call(server, ftp_command)
        if ftp_output[0].strip().endswith(b'quit'):
            logging.info('Successfully changed access rights of {0} to {1}'.format(
                    datafile, chmod_arg))
        else:
            log_ftp_output(ftp_output, logging_function=logging.warning)
