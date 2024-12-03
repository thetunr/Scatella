import numpy as np

''' Computes the pairwise distance score of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
'''
def pairwise_distance(x, y):
    ''' Recurrence matrix, redefine/use as necessary. '''
    m = np.zeros((len(x) + 1, len(y) + 1))

    for i in range(len(x) + 1):
        for j in range(len(y) + 1):
            if i == 0:
                m[i][j] = j
            elif j == 0:
                m[i][j] = i
            else: 
                diag = m[i - 1][j - 1] + ((x[i - 1] != y[j - 1]) if 1 else 0)
                right = m[i][j - 1] + 1
                down = m[i - 1][j] + 1
                m[i][j] = min(diag, right, down)
    ''' Complete this function. '''
    # print(m)
    return m[len(x)][len(y)]

''' Computes the pairwise distance matrix
Arguments: 
    sequences: the list of sequences to align
Returns:
    m: the score matrix of pairwise distances
'''
def get_matrix(sequences):
    m = np.zeros((len(sequences), len(sequences)))

    for s1 in range(len(sequences)):
        for s2 in range(s1 + 1, len(sequences)):
            score = pairwise_distance(sequences[s1], sequences[s2])
            m[s1][s2] = score
            m[s2][s1] = score
    return m


''' Converts 2d matrix to dictionary
Arguments: 
    matrix: 2d matrix
Returns:
    D: dictionary
'''
def matrix_to_dict(matrix):
    D = {}
    for i, row in enumerate(matrix):
        D[i] = {j: val for j, val in enumerate(row)}
    return D


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
def neighbor_join(D, og):
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

    for i in range(0, len(E)):
        x, y = E[i]
        if x == og or y == og:
            del E[i]
            E.append((z, x)); E.append((z, y))
            dist = uD[x][y] / 2
            uD[z] = {}
            uD[z][x] = dist; uD[z][y] = dist
            fake_root = z
            break
    
    return E, uD, fake_root

''' Performs profile-profile alignment, which aligns alignments.
'''
def profprof_alignment():
    return


def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # print(pairwise_distance(sequences[3], sequences[2]))
    m = get_matrix(sequences)
    D = matrix_to_dict(m)
    print(m)
    print(D)
    E, uD, fake_root = neighbor_join(D, 2) # have to choose outgroup

    print(E) # edges
    # print(uD)
    print(fake_root)



if __name__ == "__main__":
    main()