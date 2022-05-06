import pandas as pd


def convert_to_datetime(time):
    try:
        return pd.to_datetime(str(int(time)))
    except Exception:
        return time


def search_df(df, **kwargs):
    conditions = " & ".join(
        ["(df['{0}'] == {1})".format(col, repr(cond)) for col, cond in kwargs.items()]
    )
    return df[eval(conditions)]


def get_var_by_time(df, datetime=None, **kwargs):
    df = search_df(df, **kwargs)
    if datetime is None:
        return df
    return df[(datetime >= df.time_min) & (datetime <= df.time_max)]


class CFModelSelector:
    def __init__(self, df):
        try:
            self.df = self._update_time(df)
        except Exception:
            self.df = df

    def _update_time(self, df):
        df["time_min"] = df.time_min.apply(convert_to_datetime)
        df["time_max"] = df.time_max.apply(convert_to_datetime)
        return df

    def get_file(self, datetime=None, **kwargs):
        return get_var_by_time(self.df, datetime=datetime, **kwargs)
