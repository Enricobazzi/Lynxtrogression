import sys
import numpy as np

def split_par_file(par_file: str) -> list:
    """
    Split the par file into its sections.    
    """
    with open(par_file, 'r') as file:
        file_content = file.read()
        split_data = file_content.split('//')
        return split_data

def get_n_mig_matrices(par_file: str) -> int:
    """
    Get the number of migration matrices from the par file.
    """
    return int(split_par_file(par_file)[5].strip().split('\n')[1])

def get_matrix(matrix_string: str):
    """
    Convert the migration matrix string into a numpy array.
    """
    return np.array([s.split(' ') for s in matrix_string.strip().split('\n')[1:]], dtype=float)

def get_matrix_list(par_file: str) -> list:
    """
    Get the migration matrices from the par file.
    """
    matrix_list = []
    for n in range(6, get_n_mig_matrices(par_file)+6):
        matrix_list.append(get_matrix(split_par_file(par_file)[n]))
    return matrix_list

def get_mean_mig(matrix: np.ndarray) -> float:
    """
    Get the mean migration rate from a migration matrix.
    """
    return np.mean(matrix)

def main():
    par_file = sys.argv[1]
    mean_mig_last = get_mean_mig(get_matrix_list(par_file)[-1])
    print(mean_mig_last)

if __name__ == "__main__":
    main()
