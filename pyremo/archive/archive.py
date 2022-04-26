import os
import subprocess
import tarfile
from pathlib import Path
from warnings import warn

import pandas as pd
import parse

file_pattern = "e{usr_nr:3d}{exp_nr:3d}{type:1}{date}"
efile_pattern = "e{usr_nr:3d}{exp_nr:3d}{type:1}_c{code:3d}_{date}"

date_patterns = ["%Y", "%Y%m", "%Y%m%d", "%Y%m%d", "%Y%m%d%H", "%Y%m%d%H%M"]

patterns = [efile_pattern, file_pattern]

try:
    from tqdm import tqdm
except Exception:

    def tqdm(x):
        return x


def cdo_call(options="", op="", input="", output=None):
    call = "cdo {} {} {} {}".format(options, op, input, output)
    subprocess.Popen(call, shell=True, stdout=subprocess.PIPE).stdout.read()
    return output
    # return subprocess.run("cdo {} {} {} {}".format(options, op, input, output), shell=True)


def convert_with_cdo(f):
    """Convert a single file into NetCDF format."""
    file = os.path.basename(f)
    path = os.path.dirname(f)
    # cdo = Cdo()
    # return cdo.copy(options='-f nc', input=f, output=os.path.join(path,'nc',file+'.nc'))
    return cdo_call(
        options="-f nc",
        op="copy",
        input=f,
        output=os.path.join(path, "nc", file + ".nc"),
    )


def convert_files(files, with_dask=False):
    """Convert many files into NetCDF format."""
    results = []
    if with_dask:
        from dask import delayed
    else:

        def delayed(x):
            return x

    for f in files:
        if os.path.isfile(f):
            nc = delayed(convert_with_cdo)(f)
            results.append(nc)
    return results


def get_filename_from_archive(tar, pattern=None):
    return next(f for f in tar.getnames() if pattern in f)


def extract_file(filename, path, pattern):
    open_tar = tarfile.open(filename, "r")
    filename = get_filename_from_archive(open_tar, pattern)
    # dest = os.path.join(scratch, filename)
    open_tar.extract(filename, path=path)
    return os.path.join(path, filename)


def extract_files(tars, pattern=None, path=None, delayed=False, compute=False):
    if not isinstance(tars, list):
        tars = [tars]
    if path is None:
        path = os.path.join(os.environ["SCRATCH"], "_archive_extract")
    filenames = []
    if delayed is True:
        from dask import delayed as delay
    else:

        def delay(x):
            return x

    for tar in tars:
        if pattern is not None:
            filenames.append(delay(extract_file)(tar, path, pattern))
        else:
            filenames.append(delay(extractall)(tar, path))
    if delayed is True and compute is True:
        import dask

        return list(dask.compute(*filenames))
    return filenames


def extractall(filename, path=None):
    open_tar = tarfile.open(filename, "r")
    open_tar.extractall(path=path)
    return [os.path.join(path, f) for f in open_tar.getnames()]


def _find_files(dir, file_template="e*", exclude=None):
    """search a directory for remo output."""
    pattern = file_template
    return [path.absolute() for path in Path(dir).rglob(pattern)]


def _parse_file(file):
    file = Path(file).stem
    for pattern in patterns:
        finfo = parse.parse(pattern, file)
        if finfo:
            return finfo.named


def _parse_filepathes(filepathes, parse_dates=True):
    finfos = []
    for f in tqdm(filepathes):
        path = Path(f)
        finfo = _parse_file(path.stem)
        if finfo is None:
            warn("could not parse: ", path.stem)
            continue
        finfo["path"] = str(path)
        finfo["suffix"] = path.suffix
        if parse_dates is True:
            for dpat in date_patterns:
                try:
                    finfo["date"] = pd.to_datetime(finfo["date"], format=dpat)
                    break
                except Exception:
                    pass
        finfos.append(finfo)
    return finfos


def get_data_catalog(dir, parse_dates=False, exclude=None):
    filepathes = _find_files(dir, exclude=exclude)
    df = pd.DataFrame(_parse_filepathes(filepathes, parse_dates=parse_dates))
    try:
        df["code"] = df.code.astype(pd.Int64Dtype())
    except Exception:
        pass
    return df


def _time_range(df, time_range):
    return df[(df["date"] >= time_range[0]) & (df["date"] <= time_range[1])]


def _search(df, time_range=None, **kwargs):
    for key, value in kwargs.items():
        if value is not None:
            df = df[df[key] == value]
    if time_range is not None:
        return _time_range(df, time_range)
    return df


class RemoCatalog:
    def __init__(self, path, force_parse=False):
        self.path = path
        self.csv = os.path.join(self.path, "data_catalog.csv")
        if os.path.isfile(self.csv) and force_parse is False:
            self.df = pd.read_csv(self.csv)
        else:
            self.df = self.parse(path, parse_dates=True)
            # self.to_csv()

    def __add__(self, other):
        df = pd.concat([self.df, other.df])
        cat = RemoCatalog(self.path)
        cat.df = df
        return cat

    def parse(self, path=None, parse_dates=True):
        if path is None:
            path = self.path
        return get_data_catalog(path, parse_dates=parse_dates)

    def to_csv(self, path=None):
        self.df.to_csv(self.csv, index=False)

    def search(self, time_range=None, **kwargs):
        """arbitrary search in the catalog

        search for arbitrary keyword arguments.

        """
        return _search(self.df, time_range=time_range, **kwargs).sort_values("date")

    def archive(self, type, suffix=".tar", time_range=None, **kwargs):
        return _search(
            self.df, type=type, suffix=suffix, time_range=time_range, **kwargs
        ).sort_values("date")

    def update(self, filepathes):
        """update catalog with filepathes

        drops duplicates.

        """
        self.df = (
            pd.concat([self.df, pd.DataFrame(_parse_filepathes(filepathes))])
            .drop_duplicates()
            .reset_index(drop=True)
        )


class RemoArchive:
    def __init__(self, path, extract_path=None, force_parse=True):
        self.path = path
        if extract_path is None:
            self.extract_path = os.path.join(os.environ["SCRATCH"], "_archive_extract")
        else:
            self.extract_path = extract_path
        self.catalog = RemoCatalog(self.path, force_parse) + RemoCatalog(
            self.extract_path, force_parse=True
        )
        # self.df = get_data_catalog(self.extract_path, parse_dates=True)

    def _extractall(self, type="t", time_range=None, **kwargs):
        tars = list(self.catalog.archive(type, time_range=time_range, **kwargs).path)
        return tars
        # filepathes = extractall()

    def _extract_code(self, code, time_range=None, parallel=False, **kwargs):
        pattern = "c{:03d}".format(code)
        tars = list(self.catalog.archive("e", time_range=time_range, **kwargs).path)
        filepathes = extract_files(
            tars,
            pattern=pattern,
            path=self.extract_path,
            delayed=parallel,
            compute=parallel,
        )
        self.catalog.update(filepathes)
        return filepathes

    def first_tfile(self):
        tar = list(self.catalog.archive("t").path)[0]
        return extract_files(tar, pattern="010106")[0]

    def monthly(self, time_range=None):
        return list(
            self.catalog.search(type="m", suffix=".nc", time_range=time_range).path
        )

    def hourly(self, code=None, time_range=None):
        return list(
            self.catalog.search(
                type="e", code=code, suffix=".nc", time_range=time_range
            ).path
        )

    def tfiles(self, time_range=None):
        return list(
            self.catalog.search(type="t", suffix=".nc", time_range=time_range).path
        )
