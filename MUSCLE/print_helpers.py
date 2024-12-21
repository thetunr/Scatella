import copy
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

def transpose(l1, l2):
    l2 =[[row[i] for row in l1] for i in range(len(l1[0]))]
    return l2

''' Prints profile. '''
def print_profile_matrix(profile):
    mt = []
    mt = transpose(profile, mt)
    for i, row in enumerate(mt):
        print(nucleotide_list[i], ": ", row)

''' Prints best sequence given a profile. '''
def print_profile_sequence(profile):
    best_sequence = []
    for row in profile:
        best_sequence.append(nucleotide_list[row.index(max(row))])
    print("".join(best_sequence))
