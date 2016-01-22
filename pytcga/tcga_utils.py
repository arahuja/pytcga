import pandas as pd

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