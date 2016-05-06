import pandas as pd
import requests

# http://stackoverflow.com/questions/22604564/how-to-create-a-pandas-dataframe-from-string
import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

from .urls import CODE_TABLE_ADDRESS

def load_tcga_tabfile(path,
                      skiprows=0):
    columns = pd.read_csv(path,
                          sep='\t',
                          skiprows=skiprows,
                          nrows=10).columns
    df = pd.read_csv(path,
                    sep='\t',
                    skiprows=skiprows + 1,
                    header=0,
                    names=columns,
                    na_values='[Not Available]')

    return df

def load_studies():
    payload = {'exportType': 'csv',
               'dir': 'undefined',
               'sort': 'undefined',
               'codeTablesReport': 'bcrBatchCode'}
    r = requests.post(CODE_TABLE_ADDRESS, payload)
    df = pd.DataFrame.from_csv(StringIO(r.text))
    return df[['Study Abbreviation', 'Study Name']] \
            .drop_duplicates() \
            .sort_values(by="Study Abbreviation") \
            .reset_index(drop=True)
