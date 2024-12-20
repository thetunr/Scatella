'''
TODOS:
MAXIMUM LIKELIHOOD BRANCH OPTIMIZATION - PhyML
    1. add kmer distance/pairwise distance matrix maker function
    2. reformat sequence alignment to either be all gaps or all 'N's whichever we want (might want to swap these tasks)
    3. future step: have geneious make NJ tree to save some time and we just parse a tree from what NJ gives us
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
from collections import Counter

def read_data(distances_file):
    with open(distances_file, "r") as f:
        lines = [l.strip().split() for l in f.readlines()]
        mapping = {i: s for i, s in enumerate(lines[0])}
        lines = [l[1:] for l in lines[1:]]
        D = {i: {} for i in range(len(lines))}
        for i, l in enumerate(lines):
            for j, sval in enumerate(l):
                D[i][j] = float(sval)
    return D, mapping

''' Performs the neighbor joining algorithm on the distance matrix and the index of the outgroup species.

Arguments:
    D: A dictionary of dictionaries, defining distances between all species, every key is a species index,
        and the corresponding value is a dictionary containing all species indexes as keys. (See the description
        in `read_data` for details).
    og: outgroup index, defining which species serves as an outgroup.
        A fake root should be inserted in the middle of the pendant edge
        leading to this outgroup node.

Returns:
    E : A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
        For example [(3,1),(3,2)] represents an unrooted NJ tree of two edges, 
        3<-->1 and 3<-->2, where 1 & 2 are indexes of leaf nodes in the tree,
        and 3 is the index of the internal node you added.
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes),
        it's of the same format as D, storing all edge lengths of the NJ tree whose topology is specified by E.
        For example, {1: {1: 0.0, 2: 1.0, 3: 1.5}, 2: {1: 1.0, 2: 0.0, 3: 2.0}, 3: {1: 1.5, 2: 2.0, 3: 0.0}}
        will fully specify the edge lengths for the tree represented by the E example ([(3,1),(3,2)]):
        Length(3<-->1) = 1.5, Length(3<-->2) = 2.0.
    fake_root: which node index represent the root
'''
def neighbor_join(D):
    # init
    E = []
    uD = {i: {j: val for j, val in D[i].items()} for i in D}
    n = len(D)
    z = len(D)

    # while building tree
    while n > 2:
        # Q-criteria: 
        x = 0; y = 0
        minQ = float('inf')

        for i in D:
            for j in D[i]:
                if i < j:
                    Q = (n - 2) * D[i][j] - sum(D[i].values()) - sum(D[j].values())
                    if Q < minQ:
                        minQ = Q
                        x = i; y = j

        # Branch lengths:
        # z = {x, y}
        D[x][z] = ( D[x][y] + ( (sum(D[x].values()) - sum(D[y].values())) / (n - 2) ) ) / 2
        D[y][z] = D[x][y] - D[x][z]

        # Update distance matrix:
        D[z] = {}; uD[z] = {}
        for k in D:
            if k != x and k != y:
                d_zk = (D[x][k] + D[y][k] - D[x][y]) / 2
                D[z][k] = d_zk; D[k][z] = d_zk   # update D
                uD[z][k] = d_zk; uD[k][z] = d_zk # update uD
        uD[z][x] = D[x][z]; uD[z][y] = D[y][z]
        uD[x][z] = D[x][z]; uD[y][z] = D[y][z]

        # Deletion of rows x, y
        del D[x]; del D[y]
        # Deletion of columns x, y
        for i in D:
            if x in D[i]: del D[i][x]
            if y in D[i]: del D[i][y]

        E.append((z, x)); E.append((z, y)) # append to E
        n = len(D)
        z += 1

    # insert last edge
    last_pair = list(D.keys())
    E.append((last_pair[0], last_pair[1]))

    # IMPORTANT: code that uses the outgroup, og to create a fake root, turning the tree into a rooted tree

    # for i in range(0, len(E)):
    #     x, y = E[i]
    #     if x == og or y == og:
    #         del E[i]
    #         E.append((z, x)); E.append((z, y))
    #         dist = uD[x][y] / 2
    #         uD[z] = {}
    #         uD[z][x] = dist; uD[z][y] = dist
    #         fake_root = z
    #         break
    
    return E, uD


''' Helper function for defining a tree data structure.
    First finds the root node and find its children, and then generates 
    the whole binary tree based on.

Arguments:
    E：A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
         (See the description in `neighbor_join` for details).
    fake_root: which node index represent the root
Returns:
    tree_map：A dictionary storing the topology of the tree, where each key is a node index and each value is a list of
              node indices of the direct children nodes of the key, of at most length 2. For example, {3: [1, 2], 2: [], 1: []}
              represents the tree with 1 internal node (3) and two leaves (1, 2). the [] value shows the leaf node status.
'''
def assemble_rooted_tree(fake_root, E):
    tree_map = {}
    def find_children(parent, E):
        for l, r in E:
            if l == parent or r == parent:
                if l == parent:
                    if l in tree_map:
                        tree_map[l].append(r)
                    else:
                        tree_map[l] = [r]
                    find_children(r, list(filter(lambda x: x != (l, r), E)))
                else:
                    if r in tree_map:
                        tree_map[r].append(l)
                    else:
                        tree_map[r] = [l]
                    find_children(l, list(filter(lambda x: x != (l, r), E)))
    find_children(fake_root, E)
    return tree_map

''' Finds the post-order traversal of the tree given by tree_map
    rooted at fake_root.

Arguments:
    fake_root: which node index represent the root
    tree_map：A dictionary storing the topology of the tree (See the description in `assemble_tree` for details).
Returns:
    ordering: the ordering of nodes given by the post order traversal of the tree
'''
def get_ordering(fake_root, tree_map):
    def post_order(parent):
        l, r = tree_map[parent]
        lefts = [l]
        rights = [r]
        if l in tree_map:
            lefts = post_order(l) + lefts
        if r in tree_map:
            rights = post_order(r) + rights
        if len(lefts) > len(rights):
            return lefts + rights
        return rights + lefts
    ordering = post_order(fake_root)
    return ordering



'''
Computes the k-mer distance score of two strings.
Arguments:
    k: the value of k
    x: the first string we're aligning
    y: the second string we're aligning
Returns:
    score: the score of the optimal sequence alignment
'''
def kmerdistance(k, x, y):
  # Step 1: Extract k-mers from each string
  def get_kmers(s, k):
      """Generate all k-mers from the input string s."""
      return [s[i:i+k] for i in range(len(s) - k + 1)]
  
  # Generate k-mers for both strings
  kmers_x = get_kmers(x, k)
  kmers_y = get_kmers(y, k)
  
  # Step 2: Count occurrences of each k-mer
  count_x = Counter(kmers_x)
  count_y = Counter(kmers_y)
  
  # Step 3: Compute shared k-mers (intersection of k-mers in x and y)
  shared_kmers = sum((count_x & count_y).values())  # Intersection counts
  
  # Step 4: Compute total unique k-mers across both strings
  total_kmers_x = sum(count_x.values())
  total_kmers_y = sum(count_y.values())
  
  # Total k-mers in both sequences
  total_kmers = total_kmers_x + total_kmers_y
  
  # Step 5: Compute the k-mer distance score
  # Distance score = total kmers - 2 * shared kmers
  score = total_kmers - 2 * shared_kmers
    
  return score
    

'''
Replaces all "N"'s in a given sequence with gaps.
Arguments:
    data: sequence data (dict: name of sequence owner -> sequence)
Returns:
    replaced_data: sequence data with no N's
'''
def n_remover(data):
    replaced_data = {}
    for name, sequence in data.items():
        replaced_sequence = sequence.replace("N", "-")
        replaced_data[name] = replaced_sequence
    
    return replaced_data





''' Returns a string of the Newick tree format for the tree, rooted at a pre-defined node (fake_root).

Arguments:
    fake_root: which node index represent the root
    tree_map：A dictionary storing the topology of the tree (See the description in `assemble_tree` for details).
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes)
        (See the description in `neighbor_join` for details)
    mapping: A dictionary mapping indices to species names. (See the description in `read_data` for details)
Returns:
    output: rooted tree in Newick tree format (string). The branch lengths should be in 6 decimal digits.
'''
def generate_newick(fake_root, tree_map, uD, mapping = None):
    ''' Complete this function. '''
    output = "("
    def newick_builder(parent):
        curr = ""
        for i in tree_map[parent]:
            if i == tree_map[parent][1]:
                curr += ", "
            if i in mapping:
                curr += str(mapping[i])
            else:
                curr += "(" + newick_builder(i) + ")"
            curr += '%s:%.6f' % ("", uD[parent][i])
        return curr
    output += newick_builder(fake_root) + ");"
    return output

''' Given a distance matrix of sequences, generate a phylogeny using
neighbor-joining and output the post-order representation of the tree.

Arguments:
    dist_m: the distance matrix

Returns:
    ordering: the post-order representation of the tree
'''
def make_tree(dist_m):
    D = dict(enumerate([dict(enumerate(dist_m[i])) for i in range(len(dist_m))]))
    mapping = dict(enumerate(range(len(D))))
    og = 0
    E, uD, fake_root = neighbor_join(D, og)
    tree_map = assemble_rooted_tree(fake_root, E)
    ordering = get_ordering(fake_root, tree_map)
    print(ordering)
    nwk_str = generate_newick(fake_root, tree_map, uD, mapping)
    print(nwk_str)
    return ordering

''' Evaluates P(b|a, t) under the Jukes-Cantor model

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

def main():
    distances = [[0,5,9,9,8],[5,0,10,10,9],[9,10,0,8,7],[9,10,8,0,3],[8,9,7,3,0]]
    make_tree(distances)
    return

if __name__ == "__main__":
    main()