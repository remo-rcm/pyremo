
import os
import logging

from PyRemo.FileConventions import REMO_2015, REMO_2015_TAPE_ARCHIVE, REMO_2015_DISK_ARCHIVE
from PyRemo.DataHandler import DataHandler
from PyRemo.Scheduler import Scheduler


logging.basicConfig(level=logging.DEBUG)


SLURM_HEADER = {  'partition'   :'shared',
                  'ntasks'      :'1',
                  'mem_per_cpu' :'1280',
                  'mail_type'   :'FAIL',
                  'account'     :'ch0636',
                  'time'        :'08:00:00'}


expid = '061074'

tape_archive = REMO_2015_TAPE_ARCHIVE()
disk_archive = REMO_2015_DISK_ARCHIVE()

disk_path = '/work/ch0636/g300046/remo_results_061074'
tape_path = '/hpss/arch/ch0636/CORDEX/CORE-ATLAS/remo/exp061074'
dest_path = os.path.join(os.environ['SCRATCH'], 'remo_results_061074')

dh = DataHandler(parallel=False)#,sys='SLURM',header_dict=SLURM_HEADER)

#dh.extract_archive(expid, disk_path, years=range(2006,2007), types=['efiles'],\ 
#                   codes = [167,201], \
#                   conv=disk_archive, destdir=dest_path)

dh.extract_archive(expid, disk_path, years=range(2010,2012), types=['pfiles'],\
                   codes = [156], \
                   conv=disk_archive, destdir=dest_path)

#dh.retrieve_archive_from_tape(expid, tape_path, years=range(2006,2100), types=['efiles','pfiles'], \
#                              conv=tape_archive, destdir=dest_path)


#my_scheduler = Scheduler(sys='SLURM',logfile='DataHandler_061074_extract.jobids.ini')

#my_scheduler.log_jobs_acct() 
