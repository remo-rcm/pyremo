#! /usr/bin/python
# coding: utf-8
#
"""DataHandler
Classes and methods in :mod:`DataHandler` should help retrieving and
archiving raw REMO data.
"""

import logging
import os
import tarfile

from hpc_scheduler import Scheduler

from .file_conventions import REMO_2015_DISK_ARCHIVE, REMO_2015_TAPE_ARCHIVE
from .ftp import check_ftp_path_list, download_file

# template for executing a PyREMO command in a batch script
JOB_TEMPLATE = """
python << END
#!/usr/bin/env python
from pyremo.Ftp import download_file
from pyremo.DataHandler import DataHandler
import tarfile
dh = DataHandler(parallel=False)
{}
END
"""


class DataHandler:
    """Handles REMO disk and archive data

    Attributes
    ----------
        sys:
            The description of the scheduler used for parallel processing.
        parallel:
            True if parallel processing should be used.
        header_dict:
            a dictionory containing template replacements for the jobfile header.
        job_dir:
            Directory where job scripts should be written.
        log_dir:
            Directory where log files   should be written.

    Written by Lars Buntemeyer
    """

    def __init__(
        self, sys="SLURM", parallel=False, header_dict={}, job_dir="", log_dir=""
    ):
        self.sys = sys
        self.parallel = parallel
        self.job_dir = job_dir
        self.log_dir = log_dir
        self.header_dict = header_dict

    def _check_disk_path_list(self, file_list):
        check_list = []
        for filename in file_list:
            exists = True if os.path.isfile(filename) else False
            check_list.append((filename, exists))
        return check_list

    def _create_job(self, scheduler, jobname, commands, job_dir, log_dir):
        fill = {}  # SLURM_DEFAULTS
        fill["log_dir"] = log_dir
        fill["job_name"] = jobname
        jobscript = os.path.join(job_dir, jobname + ".sh")
        scheduler.create_job(
            jobname=jobname,
            jobscript=jobscript,
            commands=commands,
            header_dict=fill,
            write=True,
        )

    def _check_file_list(self, file_list, loc):
        check_list = []
        if loc == "tape":
            check_list = check_ftp_path_list("tape", file_list)
        elif loc == "disk":
            check_list = self._check_disk_path_list(file_list)
        else:
            raise Exception(
                "Unknown location defined in check_archive_files: " + str(loc)
            )
        return check_list

    def extract_archive(
        self,
        expid,
        datapath,
        years=[1979],
        months=range(1, 13),
        types=[],
        filters=[],
        codes=[],
        loc="",
        conv=REMO_2015_DISK_ARCHIVE(),
        destdir="",
        log_dir="",
        job_dir="",
        sc_logfile="",
    ):
        code_filters = []

        if not loc:
            loc = conv.loc
        if not destdir:
            destdir = os.getcwd()
        if codes:
            code_filters = ["_c%03d_" % int(code) for code in codes]

        job_dir = self.job_dir if not (job_dir and self.job_dir) else destdir
        log_dir = self.log_dir if not (log_dir and self.log_dir) else destdir

        # command_tpl = "files = {}                  "  + \
        #               "tar = tarfile.open(\'{}\')\n"  + \
        #               "tar.extractall(path=\'{}\')\n" + \
        #               "tar.extract(f,path=\'{}\')\n" + \
        #               "tar.close()\n"
        command_tpl = (
            "files_to_extract = {} \n"
            + "tar = tarfile.open('{}')\n"
            + "for f in files_to_extract: tarfile.extract(f,'{}')\n"
            + "tar.close()\n"
        )

        logging.info("Extracting files " + str(expid) + " archive files: " + str(types))
        logging.info("Directory: : " + datapath)

        if not self.parallel:
            for year in years:
                file_list = self.get_file_list(
                    expid, datapath, [year], months, types, conv
                )
                logging.info("extracting, year: " + str(year))
                logging.info("file list: " + str(file_list))
                for filename in file_list:
                    tar = tarfile.open(filename)
                    # members_to_extract = tar.getmembers()
                    files_to_extract = tar.getnames()
                    if code_filters:
                        # members_to_extract = [f for f in members_to_extract if any (c for c in code_filters if c in f.name )]
                        files_to_extract = [
                            f
                            for f in files_to_extract
                            if any(c for c in code_filters if c in f)
                        ]
                    if filters:
                        # members_to_extract = [f for f in members_to_extract if any (str(x) for x in filters if str(x) in f.name )]
                        files_to_extract = [
                            f
                            for f in files_to_extract
                            if any(str(x) for x in filters if str(x) in f)
                        ]
                    # logging.info('extracting files: '+str([member.name for member in members_to_extract]))
                    logging.info("extracting files: " + str(files_to_extract))
                    # tar.extractall(members=members_to_extract, path=destdir)
                    for f in files_to_extract:
                        tar.extract(f, path=destdir)
                    tar.close()
        else:
            if not sc_logfile:
                sc_logfile = os.path.join(
                    log_dir, "DataHandler_" + expid + "_extract.jobids.ini"
                )
            sc = Scheduler(
                self.sys,
                name="DataHandler",
                header_dict=self.header_dict,
                logfile=sc_logfile,
            )
            for year in years:
                file_list = self.get_file_list(
                    expid, datapath, [year], months, types, conv
                )
                py_command = ""
                for filename in file_list:
                    tar = tarfile.open(filename)
                    files_to_extract = tar.getnames()
                    if code_filters:
                        # members_to_extract = [f for f in members_to_extract if any (c for c in code_filters if c in f.name )]
                        files_to_extract = [
                            f
                            for f in files_to_extract
                            if any(c for c in code_filters if c in f)
                        ]
                    if filters:
                        # members_to_extract = [f for f in members_to_extract if any (str(x) for x in filters if str(x) in f.name )]
                        files_to_extract = [
                            f
                            for f in files_to_extract
                            if any(str(x) for x in filters if str(x) in f)
                        ]
                    # logging.info('extracting files: '+str([member.name for member in members_to_extract]))
                    py_command += command_tpl.format(
                        files_to_extract, filename, destdir
                    )
                job_command = JOB_TEMPLATE.format(py_command)
                jobname = "DataHandler_extract_" + expid + "_" + str(year)
                self._create_job(
                    sc,
                    jobname=jobname,
                    commands=job_command,
                    job_dir=job_dir,
                    log_dir=log_dir,
                )
                # sc.submit(jobname)

    def retrieve_archive_from_tape(
        self,
        expid,
        datapath,
        years=[1979],
        months=range(1, 13),
        types=[],
        loc="",
        conv=REMO_2015_TAPE_ARCHIVE(),
        destdir="",
        log_dir="",
        job_dir="",
        sc_logfile="",
    ):
        if not loc:
            loc = conv.loc
        if not destdir:
            destdir = os.getcwd()
        job_dir = self.job_dir if not job_dir and self.job_dir else os.getcwd()
        log_dir = self.log_dir if not log_dir and self.log_dir else os.getcwd()

        command_tpl = "download_file('tape','{}','{}','{}')\n"

        logging.info("Retrieve files " + str(expid) + " archive files: " + str(types))
        logging.info("Archive Directory: : " + datapath)
        logging.info("Destination Directory: : " + destdir)
        logging.info("Log Directory: : " + log_dir)
        logging.info("Job Directory: : " + job_dir)

        # create file list depending on file conventions
        if not self.parallel:
            for year in years:
                out_dir = os.path.join(destdir, str(year))
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                file_list = self.get_file_list(
                    expid, datapath, [year], months, types, conv
                )
                logging.info("downloading from hpss, year: " + str(year))
                logging.info("file list: " + str(file_list))
                # download files
                for filename in file_list:
                    download_file(
                        "tape",
                        os.path.dirname(filename),
                        os.path.basename(filename),
                        out_dir,
                    )
        else:
            if not sc_logfile:
                sc_logfile = os.path.join(
                    log_dir, "DataHandler_" + expid + "_extract.jobids.ini"
                )
            sc = Scheduler(
                self.sys,
                name="DataHandler",
                header_dict=self.header_dict,
                logfile=sc_logfile,
            )
            for year in years:
                out_dir = os.path.join(destdir, str(year))
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                file_list = self.get_file_list(
                    expid, datapath, [year], months, types, conv
                )
                py_command = ""
                for filename in file_list:
                    py_command += command_tpl.format(
                        os.path.dirname(filename), os.path.basename(filename), out_dir
                    )
                job_command = JOB_TEMPLATE.format(py_command)
                jobname = "DataHandler_retrieve_from_tape_" + expid + "_" + str(year)
                self._create_job(
                    sc,
                    jobname=jobname,
                    commands=job_command,
                    job_dir=job_dir,
                    log_dir=log_dir,
                )
                sc.submit(jobname)

    def check_archive_files(
        self,
        expid,
        datapath,
        years=[1979],
        months=range(1, 13),
        types=[],
        loc="",
        conv=REMO_2015_TAPE_ARCHIVE(),
        logfile="",
    ):
        """Checks for yearly archived data either on tape or disk.

        Parameters
        ----------
            expid:
                String containgin the experiment id.
            datapath:
                Path to the data in the archive (without the yearly subdirectory).
                (default: None)
            years:
                List of years to check
                (default: [1979])
            months:
                List of months to check
                (default: range(1,13))
            types:
                List of strings denoting filetypes
                (default: None)
            loc:
                String defining the location of data (tape, disk)
                (default: None)
            conv:
                File naming convention
                (default: REMO_2015_TAPE_ARCHIVE)

        Returns:
            missing:
                Boolean flag to indicate if all files have been found.
        """

        if not loc:
            loc = conv.loc
        # if not logfile:
        #  logfile = _create_check_file_logger(expid,expid+'_check_'+str(years[0])+'-'+str(years[-1])+'.log')

        logging.info("Checking " + str(expid) + " archive files: " + str(types))
        logging.info("Directory: : " + datapath)

        # create file list depending on file conventions
        file_list = self.get_file_list(expid, datapath, years, months, types, conv)

        # check if files exist
        check_list = self._check_file_list(file_list, loc)

        # evaluate check list
        missing = False
        for check in check_list:
            if check[1]:
                logging.info("Check of " + str(check[0]) + ": " + str(check[1]))
            else:
                missing = True
                logging.warning("Check of " + str(check[0]) + ": " + str(check[1]))

        if missing:
            logging.warning("Found missing files!")
            logging.warning("Check logfile for details: " + logfile)
        else:
            logging.info("No missing files, Congratulations!")
            logging.warning(
                "Please note, that the check does not include checks of file size!"
            )
            logging.warning(
                "Also note, that these checks are not fail save, this tool is not responsible."
            )

        return missing

    def get_file_list(
        self,
        expid,
        datapath="",
        years=[],
        months=range(1, 13),
        types=[],
        conv=REMO_2015_TAPE_ARCHIVE(),
    ):
        """Creates a list of filenames containing the datapath.

        Args:
            expid:
                String containgin the experiment id.
            datapath:
                Path to the data in the archive (without the yearly subdirectory).
                (default: None)
            years:
                List of years to check
                (default: [1979])
            months:
                List of months to check
                (default: range(1,13))
            types:
                List of strings denoting filetypes
                (default: None)
            conv:
                File naming convention
                (default: REMO_2015_TAPE_ARCHIVE)

        Returns:
            file_list:
                List of filenames.
        """

        if not types:
            types = conv.known_types
        filenames = conv.get_archive_filenames(expid, years, months, types)
        file_list = [os.path.join(datapath, f) for f in filenames]
        return file_list
