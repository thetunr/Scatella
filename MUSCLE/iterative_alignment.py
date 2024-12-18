import numpy as np
import copy
from distance_matrix import get_distance_matrix
from profile_functions import profile_alignment, combine_profiles, get_profile_sequence
from print_matrix import matrix_to_dict, dict_to_matrix, print_dict_matrix, print_profiles, print_profile_matrix, print_profile_sequence

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
    msa: list of sequences after multiple sequence alignment
'''
def progressive_alignment(sequences, root, guide_tree, ordering, cluster_counts = []):
    ordering.append(root)
    # msa = [[]] * (len(ordering))
    msa = [[]] * (2 * len(sequences) - 1)
    # msa_sequences = [[] for _ in range(len(ordering))] # MSA SEQUENCES
    msa_sequences = [[] for _ in range((2 * len(sequences) - 1))] # MSA SEQUENCES

    # Aligning profiles up guide tree before root:
    for node in ordering:
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


''' Returns aligned sequences given MSA profiles. 
Arguments: 
    sequences: sequences to be aligned
    msa: profiles
Return:
    msa_sequences: aligned sequences from profiles
'''
def get_msa_sequences(sequences, msa):
    msa_sequences = []
    for i in range(len(sequences)):
        msa_sequences.append(get_profile_sequence(msa[i]))
    return msa_sequences

'''
'''
def merge_msa(msa1, msa2):
    merged_msa = [[]] * len(msa1)
    for i in range(len(msa1)):
        if msa1[i] == []:
            merged_msa[i] = msa2[i]
        else:
            merged_msa[i] = msa1[i]
    return merged_msa
'''
Input: sequences: list of sequences (not alligned)
Output: a multiple sequence allignment
'''
def muscle(sequences):
    # Muscle Step 1: Draft Progressive
    # 1.1 k-mer counting
    kmer_matrix = matrix_to_dict(get_distance_matrix(sequences, "kmers", 5))
    # EXAMPLE: 
    kmer_matrix = matrix_to_dict([ [0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39], [31, 34, 28, 0, 43], [23, 21, 39, 43, 0] ])
    # 1.2 UPGMA
    E, uD, root, cluster_counts = upgma(kmer_matrix)
    tree1 = assemble_tree(root, E)
    # EXAMPLE: 
    tree1 = {5: [0, 1], 6: [3, 4], 7: [5, 2], 8: [6, 7]}
    cluster_counts = [1, 1, 1, 1, 1, 2, 2, 3, 5]
    # 1.3 progressive alignment
    post_ordering = get_ordering(root, tree1)
    msa1 = progressive_alignment(sequences, root, tree1, post_ordering, cluster_counts)
    msa1_sequences = get_msa_sequences(sequences, msa1)
    
    # Muscle Step 2: Improved progressive
    # 2.1 kimura distance
    kimura_matrix = matrix_to_dict(get_distance_matrix(msa1_sequences, "kimura"))
    # 2.2 UPGMA
    E, uD, root, cluster_counts = upgma(kimura_matrix)
    tree2 = assemble_tree(root, E)
    # 2.3 progressive alignment
    post_ordering = get_ordering(root, tree2)
    msa2 = progressive_alignment(sequences, root, tree2, post_ordering, cluster_counts)
    msa2_sequences = get_msa_sequences(sequences, msa2)
    print(msa2_sequences)
    # CAN BE ITERATED

    # Muscle Step 3: Refinement
    tree2 = {key: tree2[key] for key in sorted(tree2)}
    # print(tree2)
    # print(post_ordering)
    for i, node in enumerate(tree2): # For each node / two edges
        # For each leaf of the node
        for leaf in tree2[node]:
            # Split tree into two subtrees
            tree2A = {}
            seqs2A = []
            visit = [leaf]
            while visit != []:
                if visit[0] in tree2:
                    tree2A[visit[0]] = tree2[visit[0]]
                    visit += tree2[visit[0]]
                seqs2A.append(visit[0])
                visit.remove(visit[0])
            tree2B = {key: [val for val in val_list if val not in seqs2A] for key, val_list in tree2.items() if key not in tree2A}
            seqs2B = [val for val in post_ordering if val not in seqs2A]            
            seqs2A.reverse()
            # Make profiles of each half of the tree
            msa2A = progressive_alignment(sequences, leaf, tree2A, seqs2A.copy())
            msa2B = progressive_alignment(sequences, node, tree2B, seqs2B.copy())
            msa2AB = merge_msa(msa2A, msa2B)
            msa2AB_seqs = get_msa_sequences(sequences, msa2AB)
            # Re-align profiles
            x_count = len([val for val in post_ordering if val in seqs2A])
            y_count = len([val for val in post_ordering if val in seqs2B])
            score, (p_x, p_y), x_gaps, y_gaps = profile_alignment(msa2A[leaf], msa2B[node], x_count, y_count, 60)
            for seq_node in seqs2A:
                if seq_node < len(sequences):
                    for x_i in x_gaps:
                        s = msa2AB_seqs[seq_node]
                        msa2AB_seqs[seq_node] = s[:x_i] + '-' + s[x_i:]
            for seq_node in seqs2B:
                if seq_node < len(sequences):
                    for y_i in y_gaps:
                        s = msa2AB_seqs[seq_node]
                        msa2AB_seqs[seq_node] = s[:y_i] + '-' + s[y_i:]
            print("node: ", node, ", leaf: ", leaf)
            print(msa2AB_seqs)

            # Accept or reject new alignment
            
            

def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # PRINT STATEMENTS
    # print(pairwise_distance(sequences[3], sequences[2]))
    # print(kmerdistance(5, sequences[0], sequences[1]))

    muscle(sequences)
    

if __name__ == "__main__":
    main()