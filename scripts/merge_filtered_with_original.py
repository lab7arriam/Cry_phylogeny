import pandas as pd
from typing import List
import argparse


def insert_parents(ori: pd.DataFrame, f: pd.DataFrame) -> pd.DataFrame:
    def process_row(row: pd.Series) -> pd.Series:
        event_id = row.name
        rec_id = f.loc[f['event_id'] == event_id, 'rec_id'].drop_duplicates().values[0]
        rc = row.copy()
        for d in range(1, 4):
            dom = f"domain{d}"
            par = f.loc[(f['rec_id'] == rec_id) & (f['domain'] == dom), 'parent'].values[0]
            if f.loc[(f['rec_id'] == rec_id) & (f['domain'] == dom), 'unknown'].values[0] == 'yes':
                par = 'unknown'
            dev = f.loc[(f['rec_id'] == rec_id) & (f['domain'] == dom), 'event_id'].values[0]
            if f"{dev}" != f"{event_id}":
                par = f"event_{dev}({par})"
            rc[f'parent_{d}'] = par
        return rc

    return ori.apply(process_row, axis=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge original table with filtered results')
    parser.add_argument('-i', '--ifile', dest='ifile', help='original table with events', type=argparse.FileType('r'))
    parser.add_argument('-f', '--ffile', dest='ffile', help='filtered events', type=argparse.FileType('r'))
    parser.add_argument('-o', '--ofile', dest='ofile', help='output file', type=argparse.FileType('w'))

    args = parser.parse_args()

    filtered = pd.read_csv(args.ffile, sep='\t')
    original = pd.read_csv(args.ifile, sep='\t')

    (original[original.index.isin(filtered['event_id'].values)]
     .pipe(insert_parents, f=filtered)
     .to_csv(args.ofile, sep='\t'))
