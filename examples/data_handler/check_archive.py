
from pyremo.file_conventions import REMO_2015_TAPE_ARCHIVE, REMO_2015_DISK_ARCHIVE
from pyremo.data_handler import DataHandler

# define experiment id for file conventions
expid = '061074'

# path to hpss tape archive 
tape_path = '/hpss/arch/ch0636/CORDEX/CORE-ATLAS/remo/exp061074'
# path to disk archive
disk_path = '/work/ch0636/g300046/remo_results_061074'

# create a DataHandler
dh = DataHandler()

#
# check files on tape archive
#
# check all files depending on naming convention
dh.check_archive_files(expid, tape_path, years=[2005], conv=REMO_2015_TAPE_ARCHIVE())
# check only a-files
dh.check_archive_files(expid, tape_path, years=[2006], conv=REMO_2015_TAPE_ARCHIVE(), types=['afiles'])
# check only e-files for one month
dh.check_archive_files(expid, tape_path, years=[2006], months=[1], conv=REMO_2015_TAPE_ARCHIVE(), types=['afiles'])
# check only p-files for the whole range
dh.check_archive_files(expid, tape_path, years=range(2006,2100), conv=REMO_2015_TAPE_ARCHIVE(), types=['pfiles'])
#
#
# check files on disk archive
#
# check only p-files 
#dh.check_archive_files(expid, disk_path, years=range(2006,2100), conv=REMO_2015_DISK_ARCHIVE(), types=['mfiles'])
# check all for one year and print a warning if a file does not exist
check_list = dh.check_archive_files(expid, disk_path, years=[2006], conv=REMO_2015_DISK_ARCHIVE())
for check in check_list:
  if not check[1]: print('Warning, could not find file: '+check[0])


#print(REMO_2015_TAPE_ARCHIVE().FileConventions)
