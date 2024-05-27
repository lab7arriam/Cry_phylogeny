#!/usr/local/bin/python3

import argparse
import pandas as pd
from typing import List, Dict


def resort_column(df: pd.DataFrame, col: str) -> pd.DataFrame:
    df[col] = df[col].str.split(';').apply(lambda x: sorted(x)).apply(lambda x: ';'.join(x))
    return df


def major_domains(dom: str) -> List[str]:
    domains = [i for i in range(1, 4)]
    dn = int(dom[len('domain')])
    domains.remove(dn)
    return ['domain' + str(x) for x in domains]


def prepare_major(df: pd.DataFrame) -> pd.DataFrame:
    df['domain'] = df['domain'].apply(major_domains)
    df1 = df.copy()
    df2 = df.copy()
    df2['domain'] = df2['domain'].str.get(1)
    df1['domain'] = df1['domain'].str.get(0)
    frames = [df1, df2]
    return pd.concat(frames, ignore_index=True, axis=0)


def filter_partial(df: pd.DataFrame, exclude: bool) -> pd.DataFrame:
    return df[~df['Type'].str.contains('partial')] if exclude else df[df['Type'].str.contains('partial')]


def group_by_rec(df: pd.DataFrame) -> pd.DataFrame:
    def row_transform(row: pd.Series, rdict: Dict[str, int]) -> pd.Series:
        new_row = row.copy()
        rec = row['Rec'].unique()[0]
        new_row['rec_id'] = rdict[rec]
        return new_row

    recs = df['Rec'].unique()
    recs_dict = {}
    for i, r in enumerate(recs):
        recs_dict[r] = i
    return df.groupby('Rec').apply(lambda grp: row_transform(grp, recs_dict))


def prepare_df(df: pd.DataFrame) -> pd.DataFrame:
    df['event_id'] = df.index
    df['domain'] = df['Type'].apply(lambda x: x[0:7])
    df['unknown'] = 'no'
    min_df = df.copy()
    maj_df = df.copy()
    min_df['parent'] = min_df['Min_par']
    min_df['parent_type'] = 'minor'
    maj_df['parent_type'] = 'major'
    maj_df['parent'] = min_df['Maj_par']
    min_df.loc[min_df['Unknown_flag'] == 'minor', 'unknown'] = 'yes'
    maj_df.loc[maj_df['Unknown_flag'] == 'major', 'unknown'] = 'yes'
    for d in [min_df, maj_df]:
        d.drop(['Min_par', 'Maj_par'], axis=1, inplace=True)
    maj_df = prepare_major(maj_df)
    return pd.concat([min_df, maj_df], axis=0)


def set_similarity(df: pd.DataFrame) -> pd.DataFrame:
    def process_row(row: pd.Series) -> pd.Series:
        domain = row['domain'][6]
        parent = 'maj' if row['parent_type'] == 'major' else 'min'
        row['similarity'] = row['r_{0}_{1}'.format(parent, domain)]
        row['similarity'] = -row['similarity'] if row['unknown'] == 'yes' else row['similarity']
        return row

    return df.apply(process_row, axis=1)


def sort_df(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(by=['parent_type', 'similarity'], axis=0, ascending=[False, False])


def grouping(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(['Rec', 'domain'])


def subset_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df[['rec_id', 'event_id', 'Rec', 'parent', 'parent_type', 'unknown', 'domain', 'similarity']]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Restructure the recombination table')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--ofile', dest='ofile', help='output file', type=argparse.FileType('w'))
    parser.add_argument('-p', '--partials', dest='partials', help='select only partials', action='store_false', default=True)

    args = parser.parse_args()

    (pd.read_csv(args.recomb_file, sep='\t', usecols=range(0, 41))
     .pipe(resort_column, col='Rec')
     .pipe(resort_column, col='Min_par')
     .pipe(resort_column, col='Maj_par')
     .pipe(filter_partial, exclude=args.partials)
     .pipe(prepare_df)
     .pipe(set_similarity)
     .pipe(sort_df)
     .pipe(grouping)
     .pipe(group_by_rec)
     .pipe(subset_columns)
     .to_csv(args.ofile.name, index=False, sep='\t')
     )
