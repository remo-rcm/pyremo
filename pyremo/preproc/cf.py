import os
import subprocess
import tempfile
from warnings import warn

import pandas as pd
from cdo import Cdo

cdo_exe = "cdo"


def convert_to_datetime(time):
    try:
        return pd.to_datetime(str(int(time)))
    except Exception:
        return time


def search_df(df, **kwargs):
    condition_list = []
    for key, item in kwargs.items():
        if isinstance(item, list):
            cond = "(df['{0}'].isin({1}))".format(key, repr(item))
        else:
            cond = "(df['{0}'] == {1})".format(key, repr(item))
        condition_list.append(cond)
    conditions = " & ".join(condition_list)
    return df[eval(conditions)]


def get_var_by_time(df, datetime=None, **kwargs):
    df = search_df(df, **kwargs)
    if datetime is not None:
        df = df[(datetime >= df.time_min) & (datetime <= df.time_max)]
    return df


class CFModelSelector:
    def __init__(self, df, scratch=None, **kwargs):
        df = df.copy()
        if kwargs:
            df = search_df(df, **kwargs)
        try:
            self.df = self._update_time(df)
        except Exception:
            self.df = df
            warn("could not parse times in dataframe")
        self.tempfiles = []
        if scratch is None:
            try:
                scratch = os.path.join(os.environ["SCRATCH"], ".cf-selector")
            except Exception:
                pass
        self.cdo = Cdo(tempdir=scratch)

    def __repr__(self):
        return repr(self._group())

    def _repr_html_(self):
        return self._group()._repr_html_()

    def _group(self):
        return self.df.groupby("source_id")[
            [
                "variable_id",
                "experiment_id",
                "member_id",
                "institution_id",
                "table_id",
                "activity_id",
            ]
        ].agg(["unique"])

    def _update_time(self, df):
        df["time_min"] = df["time_min"].apply(convert_to_datetime)
        df["time_max"] = df.time_max.apply(convert_to_datetime)
        return df

    def get_file(self, datetime=None, **kwargs):
        sel = get_var_by_time(self.df, datetime=datetime, **kwargs)
        if len(sel.index) != 1:
            raise Exception("file selection is not unique")
        return sel.iloc[0].path

    def extract_timestep(self, datetime=None, **kwargs):
        file = self.get_file(datetime=datetime, **kwargs)
        if datetime is None:
            return file
        return self.cdo.seldate(datetime, input=file)

    # own cdo call here due to https://github.com/Try2Code/cdo-bindings/issues/34
    def _cdo_call(self, options="", op="", input="", output="temp", print_command=True):
        if output is None:
            output = ""
        elif output == "temp":
            output = tempfile.TemporaryDirectory(dir=self.scratch).name
            self.tempfiles.append(output)
        if isinstance(input, list):
            input = " ".join(input)
        call = "{} {} {} {} {}".format(cdo_exe, options, op, input, output)
        if print_command is True:
            print(call)
        stdout = subprocess.Popen(
            call, shell=True, stdout=subprocess.PIPE
        ).stdout.read()
        if output:
            return output
        return stdout
