
import pkg_resources
import pandas as pd

from ._address import _table_address_dict

#code_table = 'code-list.csv'
code_table = 'https://raw.githubusercontent.com/remo-rcm/tables/master/code-list/code-list.csv'


domain_tables = {'CORDEX': 'cordex-domains.csv',
        'AUX': 'aux-domains.csv',
        'CORDEX-FPS': 'cordex-fps.csv'}

domain_tables_external = {'CORDEX': 'https://raw.githubusercontent.com/remo-rcm/tables/master/domains/cordex-domains.csv',
        'AUX': 'https://raw.githubusercontent.com/remo-rcm/tables/master/domains/aux-domains.csv',
        'CORDEX-FPS': 'https://raw.githubusercontent.com/remo-rcm/tables/master/domains/cordex-fps.csv'}


def read_pkg_table(table, index_col=None):
    """reads a csv table from the package resource.
    """
    filename = pkg_resources.resource_stream('PyRemo.tables', table)
    return pd.read_csv(filename, index_col=index_col)

def read_table(table, index_col=None):
    """reads a csv table from the package resource.
    """
    filename = table
    return pd.read_csv(filename, index_col=index_col)


def tables():
    tables = {}
    tables['domain'] = {}
    tables['vc'] = {}
    tables['code'] = {}
    for group,address in _table_address_dict['code'].items():
        tables['code'][group] = read_table(address, index_col='code')
    for group,address in _table_address_dict['domain'].items():
        tables['domain'][group] = read_table(address, index_col='short_name')
    for group,address in _table_address_dict['vc'].items():
        tables['vc'][group] = read_table(address)
    return tables


tables = tables()
