import numpy as np
from collections import Counter


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
            score = kmerdistance(5, sequences[s1], sequences[s2])
            m[s1][s2] = score
            m[s2][s1] = score
    return m

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
    

def main():
    sequences = ["CAGGATTAG", "CAGGTTTAG", "CATTTTAG", "ACGTTAA", "ATGTTAA"]
    # print(pairwise_distance(sequences[3], sequences[2]))
    print(get_matrix(sequences))

if __name__ == "__main__":
    main()