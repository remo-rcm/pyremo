from .tables import vc

_vc_tables = vc 
table = _vc_tables


def tables():
    return _vc_tables


def table(name):
    return _vc_tables[name]


def names():
    return list(_vc_tables.keys())


class VerticalCoordinate:
    def __init__(self, name):
        self.table = table(name)

    @property
    def ak(self):
        return self.table["ak"].values

    @property
    def bk(self):
        return self.table["bk"].values

    @property
    def akbk(self):
        return self.ak, self.bk


def vc(name):
    """top level function to create a vertical coordinate object."""
    if name not in names():
        raise Exception("unknown coordinate table name, choose from {}".format(names()))
    return VerticalCoordinate(name)
