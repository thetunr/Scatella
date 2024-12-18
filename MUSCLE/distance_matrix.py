import numpy as np
from collections import Counter

''' Generate all k-mers from the input string s. '''
def get_kmers(s, k):
    return [s[i:i+k] for i in range(len(s) - k + 1)]

''' Computes the k-mer distance score of two strings.
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

''' Returns True if a and b represent a transition, False if transversion
'''
def is_transition(a, b):
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

''' Returns the kimura distance between two sequences x and y
Assumes that x and y are of the same length
'''
def kimura_distance(x, y):
    n = len(x)
    transitions = 0
    transversions= 0
    for i in range(n):
        if x[i] != "-" and y[i] != "-":
            if is_transition(x[i], y[i]):
                transitions = transitions + 1
            else:
                transversions = transversions + 1
    # p is the proportion of transitions (against the full length)
    p = transitions / n
    # q is the proportion of transversions (against the full length)
    q = transversions / n

    return -0.5 * (np.log(1 - 2 * p - q))

''' Computes the pairwise distance matrix.
Arguments: 
    sequences: the list of sequences to align
    func: string of which function to use for distances ("kimura" or "kmers")
    k: value of k for kmers if func is "kmers"
Returns:
    m: the score matrix of pairwise distances
'''
def get_distance_matrix(sequences, func, k = 0):
    m = np.zeros((len(sequences), len(sequences)))

    for s1 in range(len(sequences)):
        for s2 in range(s1 + 1, len(sequences)):
            if func == "kmers":
                score = kmerdistance(k, sequences[s1], sequences[s2])
            elif func =="kimura":
                score = kimura_distance(sequences[s1], sequences[s2])
            else:
                print("invalid function name")
            m[s1][s2] = score
            m[s2][s1] = score
    return m
