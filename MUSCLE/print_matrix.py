''' Converts 2d array matrix to 2d dictionary matrix. '''
def matrix_to_dict(matrix):
    D = {}
    for i, row in enumerate(matrix):
        D[i] = {j: val for j, val in enumerate(row)}
    return D

''' Converts 2d dictionary matrix to 2d array matrix. '''
def dict_to_matrix(D):
    num_rows = max(D.keys()) + 1
    num_cols = max(max(row.keys()) for row in D.values()) + 1 if D else 0
    matrix = [[0 for _ in range(num_cols)] for _ in range(num_rows)]
    for i, row in D.items():
        for j, val in row.items():
            matrix[i][j] = val
    return matrix

''' Prints a dictionary-based 2D matrix as a properly formatted table, 
    with rows and columns aligned by their sorted keys.
'''
def print_dict_matrix(D):    
    all_keys = sorted(D.keys())
    print("    " + " ".join(f"{key:5}" for key in all_keys))
    for row in all_keys:
        row_values = [D[row].get(col, 0) for col in all_keys] 
        print(f"{row:3} " + " ".join(f"{value:5.1f}" for value in row_values))

''' Prints list of profiles, which are 2d array matrices. '''
def print_profiles(msa):
    for i, profile in enumerate(msa):
        if profile != []:
            print(i)
            print_profile_sequence(profile)
            print_profile_matrix(profile)
            print("\n")

''' Mapping from integer to nucleotide that it represents in profile. '''
nucleotide_list = ['A', 'T', 'G', 'C', '-']

''' Prints profile. '''
def print_profile_matrix(profile):
    for i, row in enumerate(profile):
        print(nucleotide_list[i], ": ", row)


''' Prints best sequence given a profile. '''
def print_profile_sequence(profile):
    best_sequence = []
    for i in range(len(profile[0])): 
        column = [profile[j][i] for j in range(len(profile))]
        best_index = column.index(max(column))
        best_nucleotide = nucleotide_list[best_index]
        best_sequence.append(best_nucleotide)
    print("".join(best_sequence))