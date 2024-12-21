'''
TODOS:
MAXIMUM LIKELIHOOD BRANCH OPTIMIZATION - PhyML
    -- iterate over all E for loop --
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
Given sequence data, generate a phylogeny using similar methods to
PhyML

Arguments:
    sequences: a dictionary of taxon to genetic sequence
    bases: a string of the possible bases in a sequence
    size: the length of each sequence

Returns:
    ordering: the post-order representation of the tree
'''
def make_tree(sequences, bases, size, og):
    dist_m, mapping = sh.get_distances(sequences)
    E, uD, fake_root = th.neighbor_join(dist_m)
    lamb = 1
    #only edges with no swaps going here: tuples of edge and opt_result
    noswaps = []
    #only edges with swaps going here: tuples of edge, score, swap, opt_result
    swaps = []
    # root the tree to calculate total likelihood
    th.root_tree(E, og, fake_root, uD)

    curr_likelihood = th.total_likelihood(fake_root, E, size, sequences, mapping, uD, bases)
    print(curr_likelihood)
    # unroot the tree to modify branch lengths and perform NNIs
    th.unroot_tree(E, fake_root)
    for l, r in E:
        print(l,r)
        # create subtrees
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
        
        L = np.zeros((size, fake_root + 1, 5))

        # find likelihoods of subtrees
        for i in subtrees:
            th.likelihood(sequences, size, subtrees[i], treemaps[i], mapping, uD, L, bases)
        # determine optimal branch lengths and swaps
        U = np.zeros((3, size, 5))
        V = np.zeros((3, size, 5))
        # internal branch
        if 0 in subtrees and 2 in subtrees:
            w = 0
            x = 0
            y = 0
            z = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        w += L[i][u[0]][k] * th.jcm(bases[k], bases[j], uD[l][u[0]])
                        x += L[i][u[1]][k] * th.jcm(bases[k], bases[j], uD[l][u[1]])

                        y += L[i][v[0]][k] * th.jcm(bases[k], bases[j], uD[r][v[0]])
                        z += L[i][v[1]][k] * th.jcm(bases[k], bases[j], uD[r][v[1]])

                    #(a)
                    U[0][i][j] = w * x
                    V[0][i][j] = y * z
                    #(b)
                    U[1][i][j] = w * y
                    V[1][i][j] = x * z
                    #(c)
                    U[2][i][j] = w * z
                    V[2][i][j] = y * x

            def total_likelihooda(t):
                log_likelihood = 0
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            prob += 0.2 * U[0][i][j] * V[0][i][k] * th.jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_resulta = minimize_scalar(lambda x: -total_likelihooda(x), bounds=(0,2*uD[r][l]))

            #change this to make it better in the future maybe
            def total_likelihoodb(t):
                log_likelihood = 0
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            prob += 0.2 * U[1][i][j] * V[1][i][k] * th.jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_resultb = minimize_scalar(lambda x: -total_likelihoodb(x), bounds=(0,2*uD[r][l]))

            def total_likelihoodc(t):
                log_likelihood = 0
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            prob += 0.2 * U[2][i][j] * V[2][i][k] * th.jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_resultc = minimize_scalar(lambda x: -total_likelihoodc(x), bounds=(0,2*uD[r][l]))
            
            #maximized likelihoods per swap
            opta = -opt_resulta.fun
            optb = -opt_resultb.fun
            optc = -opt_resultc.fun
            if max(opta, optb, optc) == opta:
                noswaps.append(((l, r), opt_resulta.x))
            elif max(opta, optb, optc) == optb:
                score = optb - opta
                swaps.append(((l, r), score, (u[1], v[0]), opt_resultb.x))
            else:
                score = optc - opta
                swaps.append(((l, r), score, (u[1], v[1]), opt_resultc.x))        
        # external node V
        elif 0 in subtrees:
            w = 0
            x = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        #(a)
                        w += L[i][u[0]][k] * th.jcm(bases[k], bases[j], uD[l][u[0]])
                        x += L[i][u[1]][k] * th.jcm(bases[k], bases[j], uD[l][u[1]])
  
                    U[0][i][j] = w * x

            def total_likelihood(t):
                log_likelihood = 0  
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            h = 0
                            if sequences[mapping[r]][i] == bases[k]:
                                h = 1
                            prob += 0.2 * U[0][i][j] * h * th.jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_result = minimize_scalar(lambda x: -total_likelihood(x), bounds=(0,2*uD[r][l]))
            noswaps.append(((l, r), opt_result.x)) 
        #external node U
        else:
            y = 0
            z = 0
            for i in range(size):
                for j in range(5):
                    for k in range(5):
                        #(a)
                        y += L[i][v[0]][k] * th.jcm(bases[k], bases[j], uD[r][v[0]])
                        z += L[i][v[1]][k] * th.assemble_rooted_treejcm(bases[k], bases[j], uD[r][v[1]])     
                    V[0][i][j] = y * z
            def total_likelihood(t):
                log_likelihood = 0   
                for i in range(size):
                    prob = 0
                    for j in range(5):
                        for k in range(5):
                            h = 0
                            if sequences[mapping[l]][i] == bases[j]:
                                h = 1
                            prob += 0.2 * h * V[0][i][k] * th.jcm(bases[j], bases[k], t)
                    log_likelihood += np.log(prob)
                return log_likelihood
            opt_result = minimize_scalar(lambda x: -total_likelihood(x), bounds=(0,2*uD[r][l]))
            noswaps.append(((l, r), opt_result.x))

    # apply swaps and update branch lengths
    swaps = sorted(swaps, key=lambda tup: tup[1], reverse=True)
    seen = []
    swap_count = int(lamb * len(swaps))
    for e, score, swap, length in swaps:
        if swap_count == 0:
            break
        l, r = e
        swap1, swap2 = swap
        if l not in seen and r not in seen:
            print(l, r)
            seen.append(l)
            seen.append(r)
            seen.append(swap1)
            seen.append(swap2)
            if (l, swap1) in E:
                E.remove((l, swap1))
                E.append((l, swap2))
            elif (swap1, l) in E:
                E.remove((swap1, l))
                E.append((swap2, l))
            if (r, swap2) in E:
                E.remove((r, swap2))
                E.append((r, swap1))
            elif (swap2, r) in E:
                E.remove((swap2, r))
                E.append((swap1, r))
        uD[l][r] = uD[l][r] + lamb * (length - uD[l][r])
        uD[r][l] = uD[r][l] + lamb * (length - uD[r][l])
        swap_count -= 1
    for e, length in noswaps:
        l, r = e
        uD[l][r] = uD[l][r] + lamb * (length - uD[l][r])
        uD[r][l] = uD[r][l] + lamb * (length - uD[r][l])
    th.root_tree(E, og, fake_root, uD)
    new_likelihood = th.total_likelihood(fake_root, E, size, sequences, mapping, uD, bases)
    # for i in range(size):
    #     for j in range(39):
    #         for k in range(5):
    #             if new_likelihood[i][j][k] != curr_likelihood[i][j][k]:
    #                 print(new_likelihood[i][j][k], curr_likelihood[i][j][k])
    print(curr_likelihood, new_likelihood)
    return


def main():
    parser = argparse.ArgumentParser(
        description='Maximum Likelihood phylogeny on a set of sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='MUSCLE/output/aligned_1000.fasta')
    args = parser.parse_args()
    seq_file = args.f
    sequences, size = sh.read_data(seq_file)
    bases = 'ACGT-'
    og = 0
    make_tree(sequences, bases, size, og)
    return


if __name__ == "__main__":
    main()