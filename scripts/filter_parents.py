#!/usr/local/bin/python3

from ete3 import Tree
import argparse
import pandas as pd
import filter_tree as ft
from typing import List


def gca(tree: Tree, nodes: List[str]) -> Tree:
    no = [tree & n for n in nodes]
    if len(no) == 1:
        return no[0]
    return tree.get_common_ancestor(no)


def gnd(tree, nodes):
    return [tree & n for n in nodes]


def compactize(tree, nodes):
    if len(nodes) == 1:
        return True, nodes
    st = gca(tree, nodes)
    no = gnd(tree, nodes)

    if len(st.get_leaves()) == len(nodes):
        return True, nodes
    else:
        children = {ch for ch in st.get_children()}
        nd = {n for n in no}
        result = {n for n in nodes}
        if nd.issubset(children):
            return True, nodes
        else:
            df = nd.difference(children)
            for nc in df:
                while nc not in children:
                    nc = nc.up
                result.update(nc.get_leaf_names())
            return False, list(result)


def check_sister(tree, nodes1, nodes2):
    no1, no2 = gca(tree, nodes1), gca(tree, nodes2)
    return no2 in no1.get_sisters()


def sisterize(tree, nodes1, nodes2, limit = 1):
    no1, no2 = gca(tree, nodes1), gca(tree, nodes2)
    gca_both = gca(tree, nodes1 + nodes2)
    in_limit = True
    def prepare_sets(gca, n1, n2):
        ns1, ns2 = set(n1), set(n2)
        nsf1, nsf2 = set(), set()
        for c in gca.get_children():
            cset = set(c.get_leaf_names())
            if len(ns1.intersection(cset)) > 0 and len(ns2.intersection(cset)) > 0:
                return prepare_sets(c, n1, n2)
            if len(ns1.intersection(cset)) > 0:
                nsf1.update(cset)
                continue
            if len(ns2.intersection(cset)) > 0:
                nsf2.update(cset)
        return nsf1, nsf2


    for n in [no1, no2]:
        visit = False
        current = n
        for _ in range(limit):
            if current == gca_both:
                visit = True
                break
            current = current.up
        if not visit:
            in_limit = False
    if not in_limit and no1.up is not None:
        superset = set(no1.up.get_leaf_names())
        ns1, ns2 = set(no1.get_leaf_names()), set(no2.get_leaf_names())
        return ns1, ns2.intersection(superset).difference(ns1)
    else:
        return prepare_sets(gca_both, nodes1, nodes2)


def select_tree(row: pd.Series, trees: List[Tree]):
    rtype = str(row['domain'])
    if rtype.startswith('domain1'):
        return trees[0]
    elif rtype.startswith('domain2'):
        return trees[1]
    elif rtype.startswith('domain3'):
        return trees[2]
    else:
        return None


def get_process_row(trees: List[Tree], limit: int):
    def process_row(row: pd.Series):
        rec = row['Rec'].split(';')
        rec = [r.strip() for r in rec]
        tree = select_tree(row, trees)
        brec, nrec = compactize(tree, rec)
        result = row.copy()
        if row['parent'] != 'unknown' and row['parent'] != 'none':
            min = row['parent'].split(';')
            min = [m.strip() for m in min]
            bmin, nmin = compactize(tree, min)
            cnmin = list(set(nmin) - set(nrec))
            if len(cnmin) > 0:
                nrec, nmin = sisterize(tree, nrec, cnmin, limit)
            else:
                nmin = cnmin
            nmin = list(nmin)
            nmin.sort()
            result['parent'] = ';'.join(nmin) if len(nmin) > 0 else 'none'
        nrec = list(nrec)
        nrec.sort()
        result['Rec'] = ';'.join(nrec) if len(nrec) > 0 else 'none'
        return result

    return process_row


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter recombination by position on the tree')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-t', '--tree', dest='tree_file', help='trees for domains', type=argparse.FileType('r'),
                        action='append')
    parser.add_argument('-s', '--support', dest='support', help='threshold for support', type=float)
    parser.add_argument('-l', '--limit', dest='limit', help='limit sisterize union', type=int)
    parser.add_argument('-o', '--ofile', dest='ofile', help='output file', type=argparse.FileType('w'))

    args = parser.parse_args()

    trees = []

    for tf in args.tree_file:
        t = Tree(tf.name)
        trees.append(ft.filter(t, args.support))

    recomb = pd.read_csv(args.recomb_file, sep='\t')

    nrecomb = recomb.apply(get_process_row(trees, args.limit), axis=1)
    nrecomb.to_csv(args.ofile.name, sep='\t', index=False)
