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