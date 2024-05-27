import pandas as pd
from typing import List, Union, Dict
import argparse
import functools
import sys


def verbose(func):
    @functools.wraps(func)
    def r(df: pd.DataFrame, col: str, *args, **kwargs) -> pd.DataFrame:
        result = func(df, col, *args, **kwargs)
        filtered = df.loc[lambda df: ~df.index.isin(result.index.values), :]
        print(f'Filtered for "{func.__name__}" values in {col} column: {len(filtered.groupby("rec_id"))} rows')
        filtered.to_csv(sys.stdout, sep='\t')
        print('==================================')
        return result
    return r


def count_proteins(df: pd.DataFrame, col: str) -> pd.DataFrame:
    ncol = df[col].str.split(';').str.len()
    ndf = df.copy()
    ndf['count_' + col] = ncol
    return ndf


@verbose
def filter_none(df: pd.DataFrame, col: str) -> pd.DataFrame:
    return df.groupby('rec_id').filter(lambda x: "none" not in x[col].values)


def get_adjacent(domains: List[str]) -> List[str]:
    def get_name(n: int) -> str:
        return 'domain' + str(n)

    def get_num(name: str) -> int:
        return int(name[6])

    dms = [1,2,3]
    minors = [get_num(i) for i in domains]
    rest = [d for d in dms if d not in minors]
    if len(domains) == 1:
        m = minors[0]
        if m == 2:
            return [get_name(n) for n in [1, 3]]
        else:
            return [get_name(n) for n in [2]]
    elif len(domains) == 2:
        return [get_name(n) for n in rest]
    else:
        return domains


def congruency(df: pd.DataFrame, col: str) -> bool:
    minors = df.loc[lambda df: df['parent_type'] == 'minor', 'domain'].values
    adj = get_adjacent(minors)
    unidomain = adj if len(adj) == 1 and len(minors) > 1 else minors
    pairs = []
    if unidomain[0] == 'domain1':
        pairs.append(['domain1', 'domain2'])
    elif unidomain[0] == 'domain3':
        pairs.append(['domain2', 'domain3'])
    else:
        pairs.append(['domain1', 'domain2'])
        pairs.append(['domain2', 'domain3'])
    con = True
    for pair in pairs:
        vals = df.loc[lambda df: df['domain'].isin(pair), col].values
        if vals[0] != vals[1]:
            con = False
    return con


@verbose
def filter_noncongruent(df: pd.DataFrame, col: str) -> pd.DataFrame:
    return df.groupby('rec_id').filter(lambda x: congruency(x, col))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to filter events based on the phylogeny correction")
    parser.add_argument('-i', '--ifile', dest='ifn', help='file with corrected events', type=argparse.FileType('r'))
    parser.add_argument('-n', '--negative', dest='neg', help='Return filtered out entities', action='store_true')
    parser.add_argument('-o', '--ofile', dest='ofn', help='file for output', type=argparse.FileType('w'))

    args = parser.parse_args()

    (pd.read_csv(args.ifn.name, sep='\t', index_col='rec_id')
     .pipe(count_proteins, col='Rec')
     .pipe(count_proteins, col='parent')
     .pipe(filter_none, col='Rec')
     .pipe(filter_none, col='parent')
     .pipe(filter_noncongruent, col='Rec')
     .to_csv(args.ofn.name, sep='\t')
     )
