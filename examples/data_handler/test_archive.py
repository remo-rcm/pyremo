
from pyremo.file_conventions import REMO_2015_TAPE_ARCHIVE, REMO_2015_DISK_ARCHIVE
from pyremo.data_handler import DataHandler

# define experiment id for file conventions
expid = '061074'

# path to hpss tape archive 
tape_path = '/hpss/arch/ch0636/CORDEX/CORE-ATLAS/remo/exp061074'
# path to disk archive
disk_path = '/work/ch0636/g300046/remo_results_061074'


conv = REMO_2015_DISK_ARCHIVE()

print(conv.known_types)

arch_filenames = conv.get_archive_filenames(expid, types=['mfiles'])
#filenames      = conv.get_filenames(expid,types=['efiles','pfiles'])

print(arch_filenames)
#print(filenames)
#print(filenames)
