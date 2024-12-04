import numpy as np
from collections import Counter


'''
Computes the pairwise distance matrix.
Arguments: 
    sequences: the list of sequences to align
Returns:
    m: the score matrix of pairwise distances
'''
def get_matrix(sequences, func):
    m = np.zeros((len(sequences), len(sequences)))

    for s1 in range(len(sequences)):
        for s2 in range(s1 + 1, len(sequences)):
            score = func(5, sequences[s1], sequences[s2])
            m[s1][s2] = score
            m[s2][s1] = score
    return m


"""Generate all k-mers from the input string s."""
def get_kmers(s, k):
    return [s[i:i+k] for i in range(len(s) - k + 1)]


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
Returns the kimura distance between two sequences x and y
Assumes that x and y are of the same length
'''
def kimura_distance(x, y):
    """
    returns True if a and b represent a transition, False if they represent
    a transversion
    """
    def is_transition(a,b):
        if a == "A" or a == "G":
            if b == "A" or b == "G":
                return True
            else:
                return False
        else:
            if b == "C" or b == "T":
                return True
            else:
                return False
            
    n = len(x)
    transitions = 0
    transversions= 0
    for i in range(n):
        if x[i] != "-" and y[i] != "-":
            if is_transition(x[i],y[i]):
                transitions = transitions + 1
            else:
                transversions = transversions + 1
    # p is the proportion of transitions (against the full length)
    p = transitions/n
    # q is the proportion of transversions (against the full length)
    q = transversions/n

    return -0.5(np.log(1-2*p-q))


''' Creates a rooted tree using UPGMA.
Arguments:
    D: distance matrix
Returns:
    E : A list storing the edges chosen from the UPGMA algorithm in the form of tuples: (index, index). 
        For example [(3,1),(3,2)] represents an rooted UPGMA tree of two edges, 
        3<-->1 and 3<-->2, where 1 & 2 are indexes of leaf nodes in the tree,
        and 3 is the index of the internal node you added.
    uD: A dictionary of dictionary, defining distances between all nodes (leaves and internal nodes),
        it's of the same format as D, storing all edge lengths of the UPGMA tree whose topology is specified by E.
        For example, {1: {1: 0.0, 2: 1.0, 3: 1.5}, 2: {1: 1.0, 2: 0.0, 3: 2.0}, 3: {1: 1.5, 2: 2.0, 3: 0.0}}
        will fully specify the edge lengths for the tree represented by the E example ([(3,1),(3,2)]):
        Length(3<-->1) = 1.5, Length(3<-->2) = 2.0.
    root: which node index represent the root
'''
def upgma(D):
    # init
    E = []
    uD = {i: {j: val for j, val in D[i].items()} for i in D}
    n = len(D)
    z = len(D)
    cluster_counts = [1] * len(E)

    # until D is 2x2
    while n > 2:
        # score: 
        x = 0; y = 0
        min_score = float('inf')

        for i in D:
            for j in D[i]:
                if i != j:
                    score = D[i][j]
                    if score < min_score:
                        min_score = score; x = i; y = j

        # Update distance matrix:
        D[z] = {}; uD[z] = {}
        count_x = cluster_counts[x]; count_y = cluster_counts[y]
        for k in D:
            if k != x and k != y:
                d_zk = (D[x][k] * count_x + D[y][k] * count_y) / (count_x + count_y)
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
        cluster_counts.append(count_x + count_y)
        n = len(D)
        z += 1

    root = z - 1
    return E, uD, root


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
def assemble_tree(root, E):
    ''' Complete this function. '''
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


'''
'''
def profile_profile_alignment(tree):
    return


"""
Input: sequences: list of sequences (not alligned)
Output: a multiple sequence allignment
"""
def muscle(sequences):
    
    # Muscle Step 1: Draft Progressive
    # 1.1 k-mer counting
    kmer_matrix = get_matrix(sequences, kmerdistance) # k-mer distance matrix
    # 1.2 UPGMA
    E, uD, root = upgma(kmer_matrix)
    tree1 = assemble_tree(root, E) # TREE1
    # 1.3 progressive alignment
    msa1 = profile_profile_alignment(tree1) # MSA1

    # Muscle Step 2: Improved progressive
    # 2.1 kimura distance
    # use MSA1
    # 2.2 UPGMA
    # 2.3 progressive alignment
    # MSA2 is computed



    # make 2D matrix with kimura distance for each of the allignments
    M = [[0 for k in len(sequences)] for i in len(sequences)]
    for x in sequences:
        for y in sequences:
            if M[x][y] == 0:
                M[x][y] = kimura_distance(x,y)
                M[y][x] = M[x][y]

    #build a guide tree using M (same as the progressive method? or maybe w/ upgma?)
    #do progressive allignment on the guide tree
    

def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # print(pairwise_distance(sequences[3], sequences[2]))
    muscle(sequences)

if __name__ == "__main__":
    main()