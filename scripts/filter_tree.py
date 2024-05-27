#!/usr/local/bin/python3

from ete3 import Tree, TreeStyle
import argparse

def filter(tree, threshold):
    for l in tree.iter_descendants():
        if not l.is_leaf() and l.support < threshold:
            l.delete(prevent_nondicotomic=False)
    return tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add supports to the tree')
    parser.add_argument('-t', '--tree', dest='tree_file', help='path to the tree file', type=argparse.FileType('r'))
    parser.add_argument('-s', '--support', dest='support', help='threshold for support', type=float)
    parser.add_argument('-o', '--ofile', dest='ofile', help='file for output', type=argparse.FileType('w'), default='-')
    parser.add_argument('-d', '--draw', dest='draw', help='draw tree', action='store_true')

    args = parser.parse_args()

    print(args.tree_file.name)
    t = Tree(args.tree_file.name)
    tree = filter(t, args.support)

    if args.draw:
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        ts.show_branch_support = True
        t.render(args.ofile.name, w=210, units='mm', tree_style=ts)
    else:   
        args.ofile.write(tree.write())
