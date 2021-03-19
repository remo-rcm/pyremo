

import pandas as pd
from pathlib import Path
from tqdm import tqdm

import parse

file_pattern = 'e{usr_nr:3d}{exp_nr:3d}{type:1}{date}'
efile_pattern = 'e{usr_nr:3d}{exp_nr:3d}{type:1}_c{code:3d}_{date}'

date_patterns = ['%Y', '%Y%m', '%Y%m%d', '%Y%m%d', '%Y%m%d%H', '%Y%m%d%H%M']

patterns = [efile_pattern, file_pattern]

def _find_files(dir, file_template='e*', exclude=None):
    """search a directory for remo output."""
    pattern = file_template
    return [path.absolute() for path in Path(dir).rglob(pattern)]


def _parse_file(file):
    for pattern in patterns:
        finfo = parse.parse(pattern, file)
        if finfo:
            return finfo.named


def get_data_catalog(dir, parse_dates=False, exclude=None):
    filepathes = _find_files(dir, exclude)
    df_all = pd.DataFrame()
    for path in tqdm(filepathes):
        finfo = _parse_file(path.stem)
        if finfo is None:
            print('could not parse: ', path.stem)
            continue
        finfo['path'] = str(path)
        finfo['suffix'] = path.suffix
        df = pd.DataFrame(finfo, index=[0])
        if parse_dates:
            for dpat in date_patterns:
                try:
                    df['date'] = pd.to_datetime(df.date, format=dpat)
                    break
                except:
                    pass
        df_all = df_all.append(df, ignore_index=True)
    try:
        df_all['code'] = df_all.code.astype(pd.Int64Dtype())
    except:
        pass
    return df_all


def move_output_data(dir):
    pass


class RemoExpData():

    def __init__(self):
        pass




