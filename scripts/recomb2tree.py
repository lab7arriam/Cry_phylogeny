#!/usr/local/bin/python3

from ete3 import Tree, TreeNode
import argparse
import pandas as pd
import filter_tree as ft
import filter_parents as fr
from typing import List

min_dist = 0

def compactize_tree(tree: Tree, nodes_list: List[str]) -> Tree:
    nca = fr.gca(tree, nodes_list)
    if len(nca.get_leaves()) == len(nodes_list):
        return tree
    nca_children = nca.children.copy()
    nodes_set = set(nodes_list)
    new_node = TreeNode()
    for ch in nca_children:
        ch_set = {t.name for t in ch.get_leaves()}
        if ch_set <= nodes_set:
            chd = ch.dist
            chs = ch.support
            ch.detach()
            new_node.add_child(ch, dist=chd, support=chs)
    nca.add_child(new_node, dist=min_dist, support=100)
    return tree


def poly2duo(tree: Tree, recomb: List[str], par: List[str]) -> Tree:
    tree = compactize_tree(tree, recomb)
    if len(par) == 0:
        return tree
    tree = compactize_tree(tree, par)
    rca, pca = fr.gca(tree, recomb), fr.gca(tree, par)
    if rca.up != pca.up:
        return tree
    aca = rca.up
    if len(aca.children) == 2:
        return tree
    rd, rs = rca.dist, rca.support
    pcd, pcs = pca.dist, pca.support
    rca.detach()
    pca.detach()
    new_node = TreeNode()
    new_node.add_child(rca, dist=rd, support=rs)
    new_node.add_child(pca, dist=pcd, support=pcs)
    aca.add_child(new_node, dist=min_dist, support=100)
    return tree


def process_tree(tree: Tree, df: pd.DataFrame) -> Tree:
    for _, row in df.iterrows():
        rec = row['Rec'].split(';')
        rec = [r.strip() for r in rec]
        min = row['Min_par'].split(';')
        min = [m.strip() for m in min]
        if len(min) == 1 and (min[0] == 'none' or min[0] == 'unknown'):
            min = []
        tree = poly2duo(tree, rec, min)
    return tree


def subset_df(df: pd.DataFrame, domain: int) -> pd.DataFrame:
    dname = 'domain' + str(domain)
    return df.loc[lambda df: df['Type'].str.startswith(dname), ['Rec', 'Min_par']]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add recombination information to phylogenetic trees')
    parser.add_argument('-r', '--recomb', dest='recomb_file', help='path to the recombination table',
                        type=argparse.FileType('r'))
    parser.add_argument('-t', '--tree', dest='tree_file', help='tree for domains', type=argparse.FileType('r'))
    parser.add_argument('-s', '--support', dest='support', help='threshold for support', type=float)
    parser.add_argument('-d', '--domain', dest='domain', help='set domain name for recomb table', type=int)
    parser.add_argument('-o', '--ofile', dest='ofile', help='output file', type=argparse.FileType('w'))

    args = parser.parse_args()

    t = Tree(args.tree_file.name)
    tree = ft.filter(t, args.support)

    recomb = pd.read_csv(args.recomb_file, sep='\t')

    new_tree = process_tree(tree, subset_df(recomb, args.domain))

    args.ofile.write(tree.write())
