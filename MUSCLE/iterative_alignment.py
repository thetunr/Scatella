import numpy as np
from distance_matrix import get_distance_matrix
from print_matrix import matrix_to_dict, dict_to_matrix, print_dict_matrix

''' Finds the post-order traversal of the tree given by tree_map
    rooted at root.

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
    E: A list storing the edges chosen from the UPGMA algorithm in the form of tuples: (index, index). 
        For example [(3,1),(3,2)] represents an rooted UPGMA tree of two edges, 
        3<-->1 and 3<-->2, where 1 & 2 are indexes of leaf nodes in the tree,
        and 3 is the index of the internal node you added.
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes),
        it's of the same format as D, storing all edge lengths of the UPGMA tree whose topology is specified by E.
        For example, {1: {1: 0.0, 2: 1.0, 3: 1.5}, 2: {1: 1.0, 2: 0.0, 3: 2.0}, 3: {1: 1.5, 2: 2.0, 3: 0.0}}
        will fully specify the edge lengths for the tree represented by the E example ([(3,1),(3,2)]):
        Length(3<-->1) = 1.5, Length(3<-->2) = 2.0.
    root: An integer that represents which node index represent the root.


    cluster_counts: A list of leaf counts under the pertaining node index.
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
        count_x = cluster_counts[x]; count_y = cluster_counts[y]

        for k in D:
            if k != x and k != y and k != z:
                d_zk = (D[x][k] * count_x + D[y][k] * count_y) / (count_x + count_y)
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
        cluster_counts.append(count_x + count_y)

        # print(len(cluster_counts))
        # print("D:")
        # print_dict_matrix(D)
        # print("uD:")
        # print_dict_matrix(uD)
        # print("")

        n = len(D)
        z += 1
        
    root = z - 1
    return E, uD, root

''' Helper function for defining a tree data structure.
    First finds the root node and find its children, and then generates 
    the whole binary tree based on.

Arguments:
    root: which node index represent the root
    E：A list storing the edges chosen from the NJ algorithm in the form of tuples: (index, index). 
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

''' Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t):
    a_x = ""
    a_y = ""

    i = len(t) - 1
    j = len(t[0]) - 1

    while (i > 0 or j > 0):
        match t[i, j]:
            case 0: # diagonal move
                i -= 1
                j -= 1
                a_x = x[-1] + a_x
                x = x[:-1]
                a_y = y[-1] + a_y
                y = y[:-1]
            case 1: # down move
                i -= 1
                a_x = x[-1] + a_x
                x = x[:-1]
                a_y = "-" + a_y
            case -1: # right move
                j -= 1
                a_x = "-" + a_x
                a_y = y[-1] + a_y
                y = y[:-1]
    return a_x, a_y

''' Computes the score and alignment of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def sequence_alignment(x, y, s, d):
    n = len(x)
    m = len(y)

    ''' Recurrence matrix '''
    F = np.zeros([n + 1, m + 1])
    ''' Traceback matrix '''
    t = np.zeros([n + 1, m + 1])

    # initializing F matrix
    for i in range(1, n + 1):
        F[i, 0] = F[i - 1, 0] - d
        t[i, 0] = 1 # representing down move as 1

    for j in range(1, m + 1):
        F[0, j] = F[0, j - 1] - d
        t[0, j] = -1 # representing right move as -1

    # filling in F matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            first_case = F[i - 1, j - 1] + s[x[i - 1]][y[j - 1]] # diagonal
            second_case = F[i - 1, j] - d # down move
            third_case = F[i, j - 1] - d # right move

            if (first_case > second_case and first_case > third_case):
                F[i, j] = first_case
                t[i, j] = 0 # representing diagonal move as 0
            elif (second_case > third_case):
                F[i, j] = second_case
                t[i, j] = 1 # representing down move as 1
            else:  
                F[i, j] = third_case
                t[i, j] = -1 # representing right move as -1
                
    score = F[n, m]
    a_x, a_y = traceback(x, y, t)
    return score, (a_x, a_y)

''' Aligns sequences from guide tree. 
    

'''
def profile_profile_alignment(sequences, tree, ordering):

    # default score matrix from BLASTZ
    score_matrix = {"A": {"A": 91, "C": -114, "T": -123, "G": -31}, "C": {"A": -114, "C": 100, "T": -31, "G": -125}, "T": {"A": -123, "C": -31, "T": 91, "G": -114}, "G": {"A": -31, "C": -125, "T": -114, "G": 100}}

    print(score_matrix)




    return 0


'''
Input: sequences: list of sequences (not alligned)
Output: a multiple sequence allignment
'''
def muscle(sequences):
    
    # Muscle Step 1: Draft Progressive
    # 1.1 k-mer counting
    kmer_matrix = matrix_to_dict(get_distance_matrix(sequences, "kmers", 5)) # k-mer distance matrix
    # 1.2 UPGMA

    # remove later
    kmer_matrix = matrix_to_dict([ [0, 17, 21, 31, 23], [17, 0, 30, 34, 21], [21, 30, 0, 28, 39], [31, 34, 28, 0, 43], [23, 21, 39, 43, 0] ])

    E, uD, root = upgma(kmer_matrix)
    tree1 = assemble_tree(root, E) # TREE1

    # 1.3 progressive alignment
    post_ordering = get_ordering(root, tree1)
    # print(post_ordering)
    msa1 = profile_profile_alignment(sequences, tree1, post_ordering) # MSA1

    # print(E)
    # print(dict_to_matrix(uD))
    # print(root)
    # print(tree1)

    # Muscle Step 2: Improved progressive
    # 2.1 kimura distance
    # use MSA1
    # 2.2 UPGMA
    # 2.3 progressive alignment
    # MSA2 is computed



    # make 2D matrix with kimura distance for each of the allignments
    # M = [[0 for k in len(sequences)] for i in len(sequences)]
    # for x in sequences:
    #     for y in sequences:
    #         if M[x][y] == 0:
    #             M[x][y] = kimura_distance(x,y)
    #             M[y][x] = M[x][y]

    #build a guide tree using M (same as the progressive method? or maybe w/ upgma?)
    #do progressive allignment on the guide tree
    

def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # print(pairwise_distance(sequences[3], sequences[2]))
    # print(kmerdistance(5, sequences[0], sequences[1]))

    muscle(sequences)
    

if __name__ == "__main__":
    main()