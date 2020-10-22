from .tables import vc as _vc_tables

def tables():
    """Get all available vc tables.
    """
    return _vc_tables


def table(name):
    """Top level function to get a vertical coordinate table."""
    return _vc_tables[name]


def names():
    """Returns names of all available vc tables.
    """
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
    """Top level function to create a vertical coordinate object."""
    if name not in names():
        raise Exception("unknown coordinate table name, choose from {}".format(names()))
    return VerticalCoordinate(name)
