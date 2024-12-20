'''
TODOS:
MAXIMUM LIKELIHOOD BRANCH OPTIMIZATION - PhyML
    -- iterate over all E for loop --
    4. compute conditional likelihoods for each U and V possible (3 in total)
    5. determine if La, Lb, or Lc is better and compute score if swapping
    -- end loop --
    -- while loop till convergence --
    6. rank possible swaps by score, if two swaps are from adjacent branches, keep only the one with the higher score
    7. start with lambda = 1 and perform lambda swaps and set unswapped branches to l = l + lambda(la-l)
    8. divide lambda by 2 and start all over again if likelihood goes down
    -- end loop --
    9. repeat process until convergence
if we want to do more, we can do IQ-TREE which involves using a candidate set of 100 possible trees and doing the same shit
while pruning the candidate set (I don't want to do this)
'''

import numpy as np
import argparse
import tree_helpers as th
import sequence_helpers as sh


''' 
Evaluates P(b|a, t) under the Jukes-Cantor model

Arguments:
    b: descendant base (string)
    a: ancestral base (string)
    t: branch length (float)
    u: mutation rate (float, defaults to 1)
Returns:
    prob: float probability P(b|a, t)
'''
def jcm(b, a, t, u = 1.0):
    ''' Complete this function. '''
    if b == a:
        return 0.2 * (1 + (4 * np.exp(-5 * u * t)))
    return 0.2 * (1 - np.exp(-5 * u * t))


''' 
Given a distance matrix of sequences, generate a phylogeny using
neighbor-joining and output the post-order representation of the tree.

Arguments:
    dist_m: the distance matrix

Returns:
    ordering: the post-order representation of the tree
'''
def make_tree(sequences, size):
    dist_m, mapping = sh.get_distances(sequences)
    og = 0
    E, uD = th.neighbor_join(dist_m)
    for l, r in E:
        u, v = th.find_children(l, r, E)
        subtrees = {}
        for i in range(len(u)):
            subtrees[i] = th.get_ordering(u[i], th.assemble_rooted_tree(u[i], l, E))
        for i in range(len(v)):
            subtrees[2 + i] = th.get_ordering(v[i], th.assemble_rooted_tree(v[i], r, E))
        break
    
    return


def main():
    parser = argparse.ArgumentParser(
        description='Maximum Likelihood phylogeny on a set of sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='data/geneious_msa/combo - all data - realigned - 2.fasta')
    args = parser.parse_args()
    seq_file = args.f
    sequences, size = sh.read_data(seq_file)
    make_tree(sequences, size)
    return


if __name__ == "__main__":
    main()