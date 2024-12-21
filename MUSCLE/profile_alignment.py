import math
from print_helpers import print_dict_matrix, print_profiles, print_profile_matrix, print_profile_sequence

''' Default score matrix from BLASTZ. '''
score_matrix = {"A": {"A": 91, "C": -114, "T": -123, "G": -31}, "C": {"A": -114, "C": 100, "T": -31, "G": -125}, "T": {"A": -123, "C": -31, "T": 91, "G": -114}, "G": {"A": -31, "C": -125, "T": -114, "G": 100}}

''' Mapping from nucleotide to integer that represents it in profile. '''
nucleotide_mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3, '-': 4}

''' Mapping from integer to nucleotide that it represents in profile. '''
nucleotide_list = ['A', 'T', 'G', 'C', '-']


''' Computes log expectation score for two probability vectors. '''
def log_expectation_score(x_probs, y_probs, epsilon=1e-10):
    score = 0
    for x, y in zip(x_probs, y_probs):
        score += x * math.log((x + epsilon) / (y + epsilon))
    return score


''' Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first profile we're aligning
    y: the second profile we're aligning
    t: the traceback matrix
Returns:
    p_x: profile alignment of profile x
    p_y: profile alignment of profile y
    x_gaps: indices in x where gaps were inserted
    y_gaps: indices in y where gaps were inserted
'''
def profile_traceback(x, y, t):
    p_x = []
    p_y = []
    x_gaps = []
    y_gaps = []
    i = len(x)
    j = len(y)
    gap_probs = [0.0 for _ in range(len(x[0]) - 1)] + [1.0]
    while (i > 0 or j > 0):
        if t[i][j] == 0: # diagonal move
            p_x.append(x[i - 1])
            p_y.append(y[j - 1])
            i -= 1; j -= 1
        elif t[i][j] == 1:   
            p_x.append(x[i - 1])
            p_y.append(gap_probs)
            i -= 1
            y_gaps.append(i)
        else: # right move
            p_x.append(gap_probs)
            p_y.append(y[j - 1])
            j -= 1
            x_gaps.append(j)
    p_x.reverse()
    p_y.reverse()
    return p_x, p_y, x_gaps, y_gaps


''' Combines two profiles into one profile.
Arguments: 
    p_x: first profile
    p_y: second profile
    x_count: number of leaves present in profile x
    y_count: number of leaves present in profile y
Returns:
    combined_profile: combined profile
'''
def combine_profiles(p_x, p_y, x_count, y_count):
    combined_profile = [[0.0] * len(p_x[0]) for _ in range(len(p_x))]
    for i in range(len(p_x)): 
        for j in range(len(p_x[0])): 
            combined_profile[i][j] = (p_x[i][j] * x_count + p_y[i][j] * y_count) / (x_count + y_count)
    return combined_profile


''' Computes the score and profile alignment of two profile.
Arguments:
    x: the first profile we're aligning
    y: the second profile we're aligning    
    d: the gap opening/extension penalty
Returns:
    score: the score of the profile alignment
    profile: the aligned profile
'''
def profile_alignment(x, y, x_count, y_count, d):
    n = len(x)
    m = len(y)
    ''' Recurrence matrix '''
    F = [[0] * (m + 1) for _ in range(n + 1)]
    ''' Traceback matrix '''
    t = [[0] * (m + 1) for _ in range(n + 1)]

    # initializing F matrix
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] - d
        t[i][0] = 1 # representing down move as 1
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] - d
        t[0][j] = -1 # representing right move as -1
    
    # filling in F matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            if i % 500 == 0 and j == 1:
                print("    profile alignment: i -", i)

            score = log_expectation_score(x[i - 1], y[j - 1])
            # total_count = x_count + y_count
            # weighted_x_probs = [p * x_count for p in x[i - 1]]
            # total_sum_x = sum(weighted_x_probs)
            # normalized_x_probs = [p / total_sum_x for p in weighted_x_probs] if total_sum_x != 0 else [0 for _ in weighted_x_probs]
            # weighted_y_probs = [p * y_count for p in y[j - 1]]
            # total_sum_y = sum(weighted_y_probs)
            # normalized_y_probs = [p / total_sum_y for p in weighted_x_probs] if total_sum_y != 0 else [0 for _ in weighted_y_probs]
            # score = log_expectation_score(normalized_x_probs, normalized_y_probs)

            first_case = F[i - 1][j - 1] - score # diagonal
            second_case = F[i - 1][j] - d # down move
            third_case = F[i][j - 1] - d # right move

            # for traceback
            if (first_case >= second_case and first_case >= third_case):
                F[i][j] = first_case
                t[i][j] = 0 # representing diagonal move as 0
            elif (second_case >= third_case):
                F[i][j] = second_case
                t[i][j] = 1 # representing down move as 1
            else:  
                F[i][j] = third_case
                t[i][j] = -1 # representing right move as -1
    
    score = F[n][m]
    p_x, p_y, x_gaps, y_gaps = profile_traceback(x, y, t)

    return score, (p_x, p_y), x_gaps, y_gaps


''' Gets best sequence given a profile. '''
def get_profile_sequence(profile):
    best_sequence = []
    for row in profile:
        best_sequence.append(nucleotide_list[row.index(max(row))])
    return "".join(best_sequence)