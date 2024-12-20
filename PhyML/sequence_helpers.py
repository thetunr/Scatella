from collections import Counter

''' 
Reads data from ```filename``` in fasta format.

Arguments:
    filename: name of fasta file to read

Returns:
    sequences: dictionary of outputs (string (sequence id) -> sequence (string))
    size: length of each sequence
'''
def read_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        sequences = {}
        output = ""
        size = 0
        curr = ""
        for l in lines:
            if l[0] == ">": 
                if (len(output) != 0):
                    sequences[curr] = output.replace("N", "-")
                    size = len(output)
                    output = ""
                curr = l[1:].strip()
            else:
                output += l.strip()
        sequences[curr] = output.replace("N", "-")
    return sequences, size


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
Calculates the distance matrix of the sequence data using kmer distance

Arguments:
    sequences: sequence data (dict: name of sequence owner -> sequence)
    
Returns:
    D: A dictionary of dictionaries, defining distances between all species, every key is a species index,
        and the corresponding value is a dictionary containing all species indexes as keys. The values of
        these keys are the distance between species. For example {1: {1: 0.0, 2: 1.0}, 2: {1: 1.0, 2: 0.0}}
        defines the distance between two species, 1 and 2.
    mapping: A dictionary mapping indices to species names. For example, {1: 'Chimp', 2: 'Bonobo'}.
'''
def get_distances(sequences):
    D = {}
    mapping = {i: s for i, s in enumerate(sequences)}
    for i in mapping:
        D[i] = {}
        for j in mapping:
            D[i][j] = kmerdistance(4, sequences[mapping[i]], sequences[mapping[j]])
    return D, mapping
