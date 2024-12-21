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
from scipy.optimize import minimize_scalar

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

#check num of keys in subtrees to see if external or internal
''' 
Given a distance matrix of sequences, generate a phylogeny using
neighbor-joining and output the post-order representation of the tree.

Arguments:
    dist_m: the distance matrix

Returns:
    ordering: the post-order representation of the tree
'''
def make_tree(sequences, bases, size):
    dist_m, mapping = sh.get_distances(sequences)
    og = 0
    E, uD = th.neighbor_join(dist_m)
    for l, r in E:
        u, v = th.find_children(l, r, E)
        subtrees = {}
        treemaps = {}
        for i in range(len(u)):
            #indices 0, 1 = W, X
            treemaps[i] = th.assemble_rooted_tree(u[i], l, E)
            subtrees[i] = th.get_ordering(u[i], treemaps[i])
        for i in range(len(v)):
            #indices 2, 3 = Y, Z
            treemaps[2 + i] = th.assemble_rooted_tree(v[i], r, E)
            subtrees[2 + i] = th.get_ordering(v[i], treemaps[2 + i])
        
        L = np.zeros((size, 40, 5))

        for i in subtrees:
            likelihood(sequences, size, subtrees[i], treemaps[i], mapping, uD, L)
        #internal edge
        U = np.zeros((size, 5))
        V = np.zeros((size, 5))
        if 0 and 2 in subtrees:
            w = 0
            x = 0
            y = 0
            z = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        #(a)
                        w += L[i][u[0]][k] * jcm(bases[k], bases[j], uD[l][u[0]])
                        x += L[i][u[1]][k] * jcm(bases[k], bases[j], uD[l][u[1]])

                        y += L[i][v[0]][k] * jcm(bases[k], bases[j], uD[r][v[0]])
                        z += L[i][v[1]][k] * jcm(bases[k], bases[j], uD[r][v[1]])     

                    U[i][j] = w * x
                    V[i][j] = y * z
            
            def total_likelihood(t):
                log_likelihood = 0
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            prob += 0.2 * U[i][j] * V[i][k] * jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_result = minimize_scalar(lambda x: -total_likelihood(x), bounds=(0,2*uD[r][l]))
        #external node V
        elif 0 in subtrees:
            w = 0
            x = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        #(a)
                        w += L[i][u[0]][k] * jcm(bases[k], bases[j], uD[l][u[0]])
                        x += L[i][u[1]][k] * jcm(bases[k], bases[j], uD[l][u[1]])
  
                    U[i][j] = w * x

            def total_likelihood(t):
                log_likelihood = 0  
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            h = 0
                            if sequences[mapping[r]][i] == bases[k]:
                                h = 1
                            prob += 0.2 * U[i][j] * h * jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_result = minimize_scalar(lambda x: -total_likelihood(x), bounds=(0,2*uD[r][l]))
        #external node U
        else:
            y = 0
            z = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        #(a)
                        y += L[i][v[0]][k] * jcm(bases[k], bases[j], uD[r][v[0]])
                        z += L[i][v[1]][k] * jcm(bases[k], bases[j], uD[r][v[1]])     
                    V[i][j] = y * z
            def total_likelihood(t):
                log_likelihood = 0   
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            h = 0
                            if sequences[mapping[l]][i] == bases[j]:
                                h = 1
                            prob += 0.2 * h * V[i][k] * jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_result = minimize_scalar(lambda x: -total_likelihood(x), bounds=(0,2*uD[r][l]))
        break
    return


#IMPORTED FROM A3, add sup fpr gaps
''' Computes the likelihood of the data given the topology specified by ordering

Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
    seqlen: length of sequences
    ordering: postorder traversal of our topology
    bp: branch probabilities for the given branches: 6x4x4 matrix indexed as
        branch_probs[branch_id][a][b] = P(b | a, t_branch_id)
Returns:
    total_log_prob: log likelihood of the topology given the sequence data
'''
def likelihood(data, seqlen, ordering, treemap, mapping, ud, l):
    bases = 'ACGT-'
    base_map = {b: i for i, b in enumerate(bases)}

    for index in range(seqlen):
        for node_idx in ordering:
            # leaf case
            if (treemap[node_idx] == []): 
                base = data[mapping[node_idx]][index]
                for i in range(5):
                    if base_map[base] == i:
                        l[index][node_idx][i] = 1.0 
                    else:
                        l[index][node_idx][i] = 0.0
            # node with 2 children
            else:
                for i in range(5): 
                    #init probs
                    left_prob = 0.0
                    right_prob = 0.0
                    for j in range(5):
                        #find left prob
                        left_prob += jcm(bases[i], bases[j], ud[i][j]) * l[index][treemap[node_idx][0]][j]
                        #find right prob
                        right_prob += jcm(bases[i], bases[j], ud[i][j]) * l[index][treemap[node_idx][1]][j]
                    #combine for prob at node
                    l[index][node_idx][i] = left_prob * right_prob
        


def main():
    parser = argparse.ArgumentParser(
        description='Maximum Likelihood phylogeny on a set of sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='MUSCLE/output/example.fasta')
    args = parser.parse_args()
    seq_file = args.f
    sequences, size = sh.read_data(seq_file)
    bases = 'ACGT-'
    make_tree(sequences, bases, size)
    return


if __name__ == "__main__":
    main()