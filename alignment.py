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


def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # print(pairwise_distance(sequences[3], sequences[2]))
    print(get_matrix(sequences))

if __name__ == "__main__":
    main()