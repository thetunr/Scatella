from profile_alignment import get_profile_sequence

''' Returns aligned sequences given MSA profiles. 
Arguments: 
    sequences: sequences to be aligned
    msa: list of profiles
Returns:
    msa_sequences: aligned sequences from profiles
'''
def get_msa_sequences(sequences, msa):
    msa_sequences = []
    for i in range(len(sequences)):
        msa_sequences.append(get_profile_sequence(msa[i]))
    return msa_sequences


''' Merges two MSAs into one.
Arguments:
    msa1: first MSA to merge
    msa2: second MSA to merge
Returns:
    merged_msa: merged MSA for the two MSAs given (MSAs are lists of profiles)
'''
def merge_msa(msa1, msa2):
    merged_msa = [[]] * len(msa1)
    for i in range(len(msa1)):
        if msa1[i] == []:
            merged_msa[i] = msa2[i]
        else:
            merged_msa[i] = msa1[i]
    return merged_msa
