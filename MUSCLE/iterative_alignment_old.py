import numpy as np
import copy
import argparse
import json
import time
from matrix_helpers import get_distance_matrix, matrix_to_dict, dict_to_matrix
from profile_alignment import profile_alignment, combine_profiles, get_profile_sequence
from msa_helpers import get_msa_sequences, merge_msa
from refinement import tree_bipartition, insert_gaps_sequences, subst, gap_cost, compute_gap_intervals, compute_sp_score
from print_helpers import print_dict_matrix, print_profiles, print_profile_matrix, print_profile_sequence

''' Default score matrix from BLASTZ. '''
score_matrix = {"A": {"A": 91, "C": -114, "T": -123, "G": -31}, "C": {"A": -114, "C": 100, "T": -31, "G": -125}, "T": {"A": -123, "C": -31, "T": 91, "G": -114}, "G": {"A": -31, "C": -125, "T": -114, "G": 100}}

''' Mapping from nucleotide to integer that represents it in profile. '''
nucleotide_mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3, '-': 4}

''' Finds the post-order traversal of the tree given by tree_map rooted at root.
Arguments:
    root: which node index represent the root
    tree_map：A dictionary storing the topology of the tree (See the description in `assemble_tree` for details).
Returns:
    ordering: the ordering of nodes given by the post order traversal of the tree
'''
def get_ordering(root, tree_map):
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
    ordering = post_order(root)
    return ordering


''' Creates a rooted tree using UPGMA.
Arguments:
    D: distance matrix
Returns:
    E:  list storing the edges chosen from the UPGMA algorithm in the form of tuples: (index, index). 
        For example [(3,1),(3,2)] represents an rooted UPGMA tree of two edges, 
        3<-->1 and 3<-->2, where 1 & 2 are indexes of leaf nodes in the tree,
        and 3 is the index of the internal node you added.
    uD: dictionary of dictionary, defining distances between all nodes (leaves and internal nodes),
        it's of the same format as D, storing all edge lengths of the UPGMA tree whose topology is specified by E.
        For example, {1: {1: 0.0, 2: 1.0, 3: 1.5}, 2: {1: 1.0, 2: 0.0, 3: 2.0}, 3: {1: 1.5, 2: 2.0, 3: 0.0}}
        will fully specify the edge lengths for the tree represented by the E example ([(3,1),(3,2)]):
        Length(3<-->1) = 1.5, Length(3<-->2) = 2.0.
    root: integer that represents which node index represent the root.
    cluster_counts: list storing leaf count for each node.
'''
def upgma(D):
    # init
    E = []
    uD = {i: {j: val for j, val in D[i].items()} for i in D}
    n = len(D)
    z = len(D)
    cluster_counts = [1] * len(D)

    # until D is 1x1
    while n > 1:
        # Finding minimum score: 
        x = 0; y = 0
        min_score = float('inf')
        for i in D:
            for j in D[i]:
                if i < j:
                    score = D[i][j]
                    if score < min_score:
                        min_score = score; x = i; y = j

        # Update distance matrix:
        D[z] = {  }; uD[z] = {  }
        x_count = cluster_counts[x]; y_count = cluster_counts[y]

        for k in D:
            if k != x and k != y and k != z:
                d_zk = (D[x][k] * x_count + D[y][k] * y_count) / (x_count + y_count)
                D[z][k] = d_zk; D[k][z] = d_zk   # update D
                uD[z][k] = d_zk; uD[k][z] = d_zk # update uD
        
        branch_length = D[x][y] / 2
        uD[z][x] = branch_length; uD[z][y] = branch_length
        uD[x][z] = branch_length; uD[y][z] = branch_length

        # Deletion of rows x, y
        del D[x]; del D[y]
        # Deletion of columns x, y
        for i in D:
            if x in D[i]: del D[i][x]
            if y in D[i]: del D[i][y]

        E.append((z, x)); E.append((z, y)) # append to E
        cluster_counts.append(x_count + y_count)

        n = len(D)
        z += 1
        
    root = z - 1
    return E, uD, root, cluster_counts


''' Helper function for defining a tree data structure.
    First finds the root node and find its children, and then generates 
    the whole binary tree based on.
Arguments:
    root: which node index represent the root
    E： A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
         (See the description in `neighbor_join` for details).
Returns:
    tree_map：A dictionary storing the topology of the tree, where each key is a node index and each value is a list of
              node indices of the direct children nodes of the key, of at most length 2. For example, {3: [1, 2], 2: [], 1: []}
              represents the tree with 1 internal node (3) and two leaves (1, 2). the [] value shows the leaf node status.
'''
def assemble_tree(root, E):
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
    find_children(root, E)
    return tree_map


''' Aligns sequences from guide tree progressively. 
Arguments:
    sequences: list of sequences to align
    root: root of guide tree
    guide_tree: guide tree for profile profile alignment
    ordering: ordering of guide tree to align sequences by
    cluster_counts: list storing leaf count for each node
Returns:
    msa: list of profiles after multiple sequence alignment
'''
def progressive_alignment(sequences, root, guide_tree, ordering, cluster_counts = []):
    ordering.append(root)
    msa = [[]] * (2 * len(sequences) - 1)
    msa_sequences = [[] for _ in range((2 * len(sequences) - 1))] # MSA SEQUENCES

    nodes_completed = 0
    # Aligning profiles up guide tree before root:
    for node in ordering:
        print("Completed ", nodes_completed, " out of ", len(ordering))
        nodes_completed += 1

        # Leaf: 
        if node < len(sequences):
            sequence = sequences[node]
            msa[node] = [[0.0 for _ in range(len(sequence))] for _ in range(len(nucleotide_mapping))]
            msa_sequences[node].append(node)
            for j, nucleotide in enumerate(sequence):
                if nucleotide in nucleotide_mapping:
                    msa[node][nucleotide_mapping[nucleotide]][j] = 1.0
        # Node: 
        else:
            leaves = guide_tree[node]
            if len(leaves) == 1:
                leaf = leaves[0]
                msa[node] = msa[leaf]
                msa_sequences[node] = msa_sequences[leaf]
            else:
                x_leaf = leaves[0]; y_leaf = leaves[1]            
                x_count = len(msa_sequences[x_leaf]); y_count = len(msa_sequences[y_leaf]) # MSA SEQUENCES            
                score, (p_x, p_y), x_gaps, y_gaps = profile_alignment(msa[x_leaf], msa[y_leaf], x_count, y_count, 60)
                msa[node] = combine_profiles(p_x, p_y, x_count, y_count)
                msa_sequences[node] = msa_sequences[x_leaf] + msa_sequences[y_leaf]

                gap_probs = [0.0 for _ in range(len(nucleotide_mapping) - 1)] + [1.0]
                for node in msa_sequences[x_leaf]:
                    for x_i in x_gaps:
                        for row, value in zip(msa[node], gap_probs):
                            row.insert(x_i, value)
                for node in msa_sequences[y_leaf]:
                    for y_i in y_gaps:
                        for row, value in zip(msa[node], gap_probs):
                            row.insert(y_i, value)

    return msa


''' Performs MUSCLE on list of sequences. 
Arguments: 
    sequences: list of sequences
Return:
    msa_seqs: list of aligned sequences
'''
def muscle(sequences):
    # MUSCLE Step 1: Draft Progressive
        print("MUSCLE Step 1 begin.")
        start_time = time.time()
        # 1.1 k-mer counting
        kmer_matrix = matrix_to_dict(get_distance_matrix(sequences, "kmers", 5))
        # 1.2 UPGMA
        E, uD, root, cluster_counts = upgma(kmer_matrix)
        tree1 = assemble_tree(root, E)
        # 1.3 progressive alignment
        post_ordering_tree1 = get_ordering(root, tree1)
        print("entering progressive_alignment")
        msa1 = progressive_alignment(sequences, root, tree1, post_ordering_tree1, cluster_counts)
        print("completed progressive_alignment")
        msa1_seqs = get_msa_sequences(sequences, msa1)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"MUSCLE Step 1 completed in {elapsed_time:.2f} seconds.\n---------------------")
    
    # MUSCLE Step 2: Improved progressive (CAN BE ITERATED)
        print("MUSCLE Step 2 begin.")
        start_time = time.time()
        # 2.1 kimura distance
        kimura_matrix = matrix_to_dict(get_distance_matrix(msa1_seqs, "kimura"))
        # 2.2 UPGMA
        E, uD, root, cluster_counts = upgma(kimura_matrix)
        tree2 = assemble_tree(root, E)
        # 2.3 progressive alignment
        post_ordering_tree2 = get_ordering(root, tree2)
        msa2 = progressive_alignment(sequences, root, tree2, post_ordering_tree2, cluster_counts)
        msa2_seqs = get_msa_sequences(sequences, msa2)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"MUSCLE Step 2 completed in {elapsed_time:.2f} seconds.\n---------------------")

    # MUSCLE Step 3: Refinement
        print("MUSCLE Step 3 begin.")
        start_time = time.time()
        best_msa = msa2_seqs
        sp_max = compute_sp_score(msa2_seqs, subst, gap_cost)
        tree2 = {key: tree2[key] for key in sorted(tree2)}
        for node in tree2: # For each node / two edges
            # For each leaf of the node
            for leaf in tree2[node]:
                # Split tree into two subtrees
                tree2A, tree2B, seqs2A, seqs2B = tree_bipartition(leaf, tree2, post_ordering_tree2)
                # Make profiles of each half of the tree
                msa2A = progressive_alignment(sequences, leaf, tree2A, seqs2A.copy())
                msa2B = progressive_alignment(sequences, node, tree2B, seqs2B.copy())
                msa2AB = merge_msa(msa2A, msa2B)
                msa2AB_seqs = get_msa_sequences(sequences, msa2AB)
                # Re-align profiles
                x_count = len([val for val in post_ordering_tree2 if val in seqs2A])
                y_count = len([val for val in post_ordering_tree2 if val in seqs2B])
                score, (p_x, p_y), x_gaps, y_gaps = profile_alignment(msa2A[leaf], msa2B[node], x_count, y_count, 60)
                msa2AB_seqs = insert_gaps_sequences(msa2AB_seqs, seqs2A, x_gaps, len(sequences))
                msa2AB_seqs = insert_gaps_sequences(msa2AB_seqs, seqs2B, y_gaps, len(sequences))
                # Accept or reject new alignment
                sp_score = compute_sp_score(msa2AB_seqs, subst, gap_cost)
                # PRINT STATEMENTS
                # print("node: ", node, ", leaf: ", leaf)
                # print(msa2AB_seqs)
                # print("SP score: ", sp_score)
                if sp_score > sp_max:
                    print("sp max changed")
                    sp_max = sp_score
                    best_msa = msa2AB_seqs
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"MUSCLE Step 3 completed in {elapsed_time:.2f} seconds.\n---------------------")
        return best_msa


''' Parses FASTA file and returns list of sequences, ignoring descriptions.
Arguments:
    faste_file: directory of the file
Returns:
    sequences: list of sequences as a list of strings
'''
def parse_fasta(fasta_file):
    sequences = []
    with open(fasta_file) as f:
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ""
            else:
                sequence += line
        if sequence: 
            sequences.append(sequence)
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for multiple sequences.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)
    parser.add_argument('-e', action="store", dest="e", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d
    e = args.e

    sequences = parse_fasta(fasta_file)
    print("Sequences parsed from FASTA:")
    print(len(sequences))
    # for i, sequence in enumerate(sequences):
    #     print(f"Sequence {i+1}: Length = {len(sequence)}")

    # with open(score_matrix_file) as f:
    #     s = json.loads(f.read().strip())

    # sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    msa_seqs = muscle(sequences)
    print(msa_seqs)
    

if __name__ == "__main__":
    main()




# EXAMPLE: 
# kmer_matrix = matrix_to_dict([ [0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39], [31, 34, 28, 0, 43], [23, 21, 39, 43, 0] ])
# EXAMPLE: 
# tree1 = {5: [0, 1], 6: [3, 4], 7: [5, 2], 8: [6, 7]}
# cluster_counts = [1, 1, 1, 1, 1, 2, 2, 3, 5]