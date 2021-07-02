
import pandas as pd

from ._resources import read_remo_domain_tables, read_remo_code_tables, read_remo_vc_tables


class read_cls():
    
    def __init__(self, reader):
        self.reader = reader
    
    @property
    def tables(self):
        return self.reader()
    
    @property
    def table(self):
        return pd.concat(self.tables.values())
    
    def __getattr__(self, table):
        return self.tables[table]
    

class codes_cls():

    @property
    def tables(self):
        return read_remo_code_tables()


class vc_cls():

    @property
    def tables(self):
        return read_remo_vc_tables()

    
domains = read_cls(read_remo_domain_tables)
codes = read_cls(read_remo_code_tables)
vc = read_cls(read_remo_vc_tables)
