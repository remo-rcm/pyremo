#! /usr/bin/python
# coding: utf-8
#

"""File Conventions

Classes and methods in :mod:`FileConventions` define a number of
file naming conventions that can be used by the DataHandler to
store and retrieve REMO data from the archive or the filesystem.

"""


import os

FILE_FREQS = ["1h", "6h", "day", "mon"]
ARCH_FREQS = ["mon", "year"]


class FileArchiveConvention:
    def __init__(self, file_list=[], freq="month", filename=""):
        self.file_list = file_list
        self.freq = freq
        self.filename = filename

    def tar():
        pass

    def untar(location=""):
        pass


class FileConvention:
    def __init__(
        self,
        file_freq="",
        arch_freq="mon",
        codes=[],
        file_formatter="",
        arch_formatter="",
    ):
        self.file_freq = file_freq
        self.arch_freq = arch_freq
        self.arch_formatter = arch_formatter
        self.file_formatter = file_formatter
        self.codes = codes


class REMO_2015:
    """Class to hold classical REMO file conventions

    Written by Lars Buntemeyer

    Last modified: 06.03.2019
    """

    def __init__(
        self,
        file_suffix="",
        arch_suffix=".tar",
        year_path="year%04d",
        pcodes=[0],
        ecodes=[0],
    ):
        # suffix of REMO data archive files
        self.arch_suffix = arch_suffix
        self.file_suffix = file_suffix
        # subdir convention for storing yearly data
        self.year_path = year_path

        self.ecodes = ecodes
        self.pcodes = pcodes

        self.afile_prefix = "a%s"
        self.efile_prefix = "e%s"

        self.year_conv = "%04d"
        self.month_conv = "%02d"
        self.day_conv = "%02d"
        self.hour_conv = "%02d"

        self.hourly_conv = "%04d%02d%02d%02d"  # year, month, day, hour
        self.daily_conv = "%04d%02d%02d"  # year, month, day
        self.monthly_conv = "%04d%02d"  # year, month
        self.yearly_conv = "%04d"  # year

        self.afile_base = self.afile_prefix + "a"
        self.tfile_base = self.efile_prefix + "t"
        self.efile_base = self.efile_prefix + "e"
        self.ffile_base = self.efile_prefix + "f"
        self.mfile_base = self.efile_prefix + "m"
        self.nfile_base = self.efile_prefix + "n"
        self.sfile_base = self.efile_prefix + "s"
        self.pfile_base = self.efile_prefix + "p"
        self.mnsfile_base = self.efile_prefix + "mns"

        self.arch_suffix = ".tar"

        arch_formatter = self.afile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.afile_base + self.hourly_conv + self.file_suffix
        self.afiles = FileConvention(
            file_freq="6h",
            arch_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.tfile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.tfile_base + self.hourly_conv + self.file_suffix
        self.tfiles = FileConvention(
            file_freq="6h",
            arch_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.efile_base + self.monthly_conv + self.arch_suffix
        file_formatter = (
            self.efile_base + "_c%03d_" + self.monthly_conv + self.file_suffix
        )
        self.efiles = FileConvention(
            file_freq="mon",
            arch_freq="mon",
            codes=self.ecodes,
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.pfile_base + self.yearly_conv + self.arch_suffix
        file_formatter = (
            self.pfile_base + "_c%03d_" + self.monthly_conv + self.file_suffix
        )
        self.pfiles = FileConvention(
            file_freq="mon",
            arch_freq="year",
            codes=self.pcodes,
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        self.FileConventions = {
            "afiles": self.afiles,
            "tfiles": self.tfiles,
            "efiles": self.efiles,
            "pfiles": self.pfiles,
        }

        self.known_types = self.FileConventions.keys()

    def get_archive_filenames(
        self, expid, years=[1979], months=range(1, 12 + 1), types=[]
    ):
        if not types:
            types = self.known_types
        filenames = []
        for ftype in types:
            conv = self.FileConventions[ftype]
            for year in years:
                year_path = self.year_path % year
                if conv.arch_freq == "year":
                    fill = (expid, year)
                    filenames.append(
                        os.path.join(
                            year_path, self.FileConventions[ftype].arch_formatter % fill
                        )
                    )
                elif conv.arch_freq == "mon":
                    for month in months:
                        fill = (expid, year, month)
                        filenames.append(
                            os.path.join(
                                year_path,
                                self.FileConventions[ftype].arch_formatter % fill,
                            )
                        )
        return filenames

    def get_filenames(
        self,
        expid,
        years=[1979],
        months=range(1, 12 + 1),
        types=[],
        codes=[],
        path="",
        year_in_path=False,
    ):
        if not types:
            types = self.known_types
        filenames = []
        for ftype in types:
            conv = self.FileConventions[ftype]
            if codes:
                for code in codes:
                    for year in years:
                        year_path = self.year_path % year
                        if year_in_path:
                            prefix_path = os.path.join(path, year_path)
                        else:
                            prefix_path = path
                        if conv.file_freq == "year":
                            fill = (expid, code, year)
                            filenames.append(
                                os.path.join(
                                    prefix_path,
                                    self.FileConventions[ftype].file_formatter % fill,
                                )
                            )
                        elif conv.file_freq == "mon":
                            for month in months:
                                fill = (expid, code, year, month)
                                filenames.append(
                                    os.path.join(
                                        prefix_path,
                                        self.FileConventions[ftype].file_formatter
                                        % fill,
                                    )
                                )
            else:
                for year in years:
                    year_path = self.year_path % year
                    if year_in_path:
                        prefix_path = os.path.join(path, year_path)
                    else:
                        prefix_path = path
                    if conv.file_freq == "year":
                        fill = (expid, year)
                        filenames.append(
                            os.path.join(
                                prefix_path,
                                self.FileConventions[ftype].file_formatter % fill,
                            )
                        )
                    elif conv.file_freq == "mon":
                        for month in months:
                            fill = (expid, year, month)
                            print(fill)
                            print(conv.file_formatter)
                            filenames.append(
                                os.path.join(
                                    prefix_path,
                                    self.FileConventions[ftype].file_formatter % fill,
                                )
                            )
        return filenames


class REMO_2015_DISK_ARCHIVE(REMO_2015):
    """Class to hold REMO convetions for archiving output on the disk.

    Written by Lars Buntemeyer

    Last modified: 06.03.2019
    """

    def __init__(
        self,
        file_suffix="",
        arch_suffix=".tar",
        year_path="%04d",
        pcodes=[0],
        ecodes=[0],
    ):
        REMO_2015.__init__(self, file_suffix, arch_suffix, year_path, pcodes, ecodes)
        self.file_suffix = ""
        self.loc = "disk"

        # define special formats for data archived on disk

        arch_formatter = self.ffile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.ffile_base + self.hourly_conv + self.file_suffix
        self.ffiles = FileConvention(
            file_freq="mon",
            arch_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.mfile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.mfile_base + self.monthly_conv + self.file_suffix
        self.mfiles = FileConvention(
            file_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.sfile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.sfile_base + self.monthly_conv + self.file_suffix
        self.sfiles = FileConvention(
            file_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.nfile_base + self.monthly_conv + self.arch_suffix
        file_formatter = self.nfile_base + self.monthly_conv + self.file_suffix
        self.nfiles = FileConvention(
            file_freq="mon",
            arch_freq="mon",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        self.FileConventions["ffiles"] = self.ffiles
        self.FileConventions["sfiles"] = self.sfiles
        self.FileConventions["mfiles"] = self.mfiles
        self.FileConventions["nfiles"] = self.nfiles

        self.known_types = self.FileConventions.keys()


class REMO_2015_TAPE_ARCHIVE(REMO_2015):
    """Class to hold REMO convetions for archiving output in a tape archive.
       File conventions includes an mns archive file and the yearly archived
       restart file convention.

    Written by Lars Buntemeyer

    Last modified: 06.03.2019
    """

    def __init__(
        self,
        file_suffix="",
        arch_suffix=".tar",
        year_path="year%04d",
        pcodes=[0],
        ecodes=[0],
    ):
        REMO_2015.__init__(self, file_suffix, arch_suffix, year_path, pcodes, ecodes)

        self.loc = "tape"

        arch_formatter = self.mnsfile_base + self.yearly_conv + self.arch_suffix
        file_formatter = self.mnsfile_base + self.monthly_conv + self.file_suffix
        self.mnsfiles = FileConvention(
            file_freq="mon",
            arch_freq="year",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        arch_formatter = self.ffile_base + self.yearly_conv + self.arch_suffix
        file_formatter = self.ffile_base + self.hourly_conv + self.file_suffix
        self.ffiles = FileConvention(
            file_freq="1h",
            arch_freq="year",
            arch_formatter=arch_formatter,
            file_formatter=file_formatter,
        )

        self.FileConventions["mnsfiles"] = self.mnsfiles
        self.FileConventions["ffiles"] = self.ffiles

        self.known_types = self.FileConventions.keys()


# REMO_2015 = REMO_2015()
# REMO_2015_DISK_ARCHIVE = REMO_2015_DISK_ARCHIVE()
# REMO_2015_TAPE_ARCHIVE = REMO_2015_TAPE_ARCHIVE()

CONVENTIONS = {
    "REMO_2015": REMO_2015(),
    "REMO_2015_DISK_ARCHIVE": REMO_2015_DISK_ARCHIVE(),
    "REMO_2015_TAPE_ARCHIVE": REMO_2015_TAPE_ARCHIVE(),
}
()


def get_convention(name):
    return CONVENTIONS[name]
