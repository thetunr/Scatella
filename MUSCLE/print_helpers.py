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