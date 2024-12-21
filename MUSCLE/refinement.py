''' Partitions tree into two subsets given a leaf to disconnect. 
Arguments:
    leaf: node to disconnect from rest of tree
    tree: entire tree structure
    post_ordering: list of post-order traversal of tree
Returns
    tree2A: subtree including leaf and all of leaf's children
    tree2B: rest of tree not included in tree2A
    seqs2A: sequences of leaves in tree2A
    seqs2B: sequences of leaves in tree2B
'''
def tree_bipartition(leaf, tree, post_ordering):
    tree2A = {}
    seqs2A = []
    visit = [leaf]
    while visit != []:
        if visit[0] in tree:
            tree2A[visit[0]] = tree[visit[0]]
            visit += tree[visit[0]]
        seqs2A.append(visit[0])
        visit.remove(visit[0])
    tree2B = {key: [val for val in val_list if val not in seqs2A] for key, val_list in tree.items() if key not in tree2A}
    seqs2B = [val for val in post_ordering if val not in seqs2A]            
    seqs2A.reverse()
    return tree2A, tree2B, seqs2A, seqs2B


''' Inserts gaps into sequences.
Arguments: 
    msa_seqs: sequences to add gaps to
    sequences: which sequences of msa_seqs to affect
    gaps: list of new gap locations to be inserted
    num_leaves: how many leaves there are in msa_seqs
Returns:
    msa_seqs: sequences after adding gaps
'''
def insert_gaps_sequences(msa_seqs, sequences, gaps, num_leaves):
    inserted = 0
    gaps = sorted(gaps)
    for seq_node in sequences:
        if seq_node < num_leaves:
            for i in gaps:
                s = msa_seqs[seq_node]
                msa_seqs[seq_node] = s[:i + inserted] + '-' + s[i + inserted:]
    return msa_seqs


''' Substitution scoring function. '''
def subst(x, y):
    return 1 if x == y else -1

''' Gap cost function. '''
def gap_cost(interval):
    return interval['length']


''' Computes gap intervals for a sequence.
Arguments:
    seq: sequence
Return:
    lg: list of gap intervals for sequence
'''
def compute_gap_intervals(seq):
    lg = []
    ig = None
    for i in range(len(seq)):
        if seq[i] == '-' and ig is None:  # Start a new gap interval
            ig = {'start': i}
        elif seq[i] != '-' and ig is not None:  # End current gap interval
            ig['end'] = i - 1
            ig['length'] = ig['end'] - ig['start'] + 1
            lg.append(ig)
            ig = None
    # Handle terminal gap
    if ig is not None:
        ig['end'] = len(seq) - 1
        ig['length'] = ig['end'] - ig['start'] + 1
        lg.append(ig)
    return lg

''' Compute the SP score for a given set of aligned sequences.
Arguments: 
    msa_seqs: list of sequence strings to compute SP score of
    subst: substitution scoring function for two characters
    gap_cost: gap cost function
Returns:
    float: SP score for the alignment
'''
def compute_sp_score(msa_seqs, subst, gap_cost):
    sp_score = 0  # Initialize total SP score
    # Iterate over all sequence pairs
    for i in range(len(msa_seqs)):
        for j in range(i + 1, len(msa_seqs)):
            seq_i = msa_seqs[i]
            seq_j = msa_seqs[j]
            if len(seq_i) != len(seq_j):
                ValueError("Sequence ", i, " does not have the same length as sequence ", j)
            # Compute substitution costs
            for k in range(len(seq_i)):
                if seq_i[k] != '-' and seq_j[k] != '-':
                    sp_score += subst(seq_i[k], seq_j[k])
    # Add gap penalties
    for seq in msa_seqs:
        for gap_interval in compute_gap_intervals(seq):
            sp_score += gap_cost(gap_interval)
    return sp_score
