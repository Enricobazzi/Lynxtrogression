import os
import argparse

def get_observed_sfs(log_file):
    """
    Extracts the observed site frequency spectrum (SFS) from a log file.

    Parameters:
    log_file (str): The path to the log file.

    Returns:
    list: The observed SFS as a list of floats.
    """
    observed_sfs = []  # Initialize the variable with an empty list
    with open(log_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("REPLICATION FOUR_SFS"):
                observed_sfs = [float(bin.strip()) for bin in lines[i].split('REPLICATION FOUR_SFS')[1].strip().strip('[]').split(',')]
    if observed_sfs:
        return observed_sfs
    else:
        return ["NA", "NA", "NA", "NA", "NA", "NA", "NA"]

def get_simulated_sfs(log_file):
    """
    Retrieves the simulated site frequency spectrum (SFS) from a log file.

    Parameters:
    log_file (str): The path to the log file containing the simulated SFS.

    Returns:
    list: The simulated SFS as a list of floating-point values.
    """
    simulated_sfs = []  # Initialize the variable with an empty list
    with open(log_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("BEST SOLUTION"):
                simulated_sfs = [float(bin.strip()) for bin in lines[i-1].strip().strip('[]').split(',')]
    if simulated_sfs:
        return simulated_sfs
    else:
        return ["NA", "NA", "NA", "NA", "NA", "NA", "NA"]

def get_fitness(log_file):
    """
    Retrieves the fitness values from a log file.

    Parameters:
    log_file (str): The path to the log file containing the fitness values.

    Returns:
    list: The fitness values as a list of floating-point values.
    """
    fitness = None
    with open(log_file, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.startswith("BEST SOLUTION"):
                fitness = float(line.split(' ')[4])
    if fitness:
        return fitness
    else:
        return "NA"
    
def get_run_number(log_file):
    """
    Retrieves the run number from a log file.

    Parameters:
    log_file (str): The path to the log file containing the run number.

    Returns:
    int: The run number.
    """

    return int(log_file.split('.log')[0].split('_')[-1])

def print_header():
    """
    Prints the header for the data columns.

    The header includes the column names for the run number, sfs bins, 
    fitness value and data type.
    """
    print("run_number", end=",")
    print("0_1", end=",")
    print("0_2", end=",")
    print("1_0", end=",")
    print("1_1", end=",")
    print("1_2", end=",")
    print("2_0", end=",")
    print("2_1", end=",")
    print("fitness", end=",")
    print("type")

def print_row(log_file):
    """
    Print one row for observed data from a log file and one
    for simulated data.

    Parameters:
    log_file (str): The path to the log file containing the SFS data.

    Returns:
    None
    """
    observed_sfs = get_observed_sfs(log_file)
    simulated_sfs = get_simulated_sfs(log_file)
    fitness = get_fitness(log_file)
    run_number = get_run_number(log_file)

    print(run_number, end=",")
    for sfs in observed_sfs:
        print(sfs, end=",")
    print("0.0", end=",")
    print("observed")

    print(run_number, end=",")
    for sfs in simulated_sfs:
        print(sfs, end=",")
    print(fitness, end=",")
    print("simulated")

def get_log_files(dir):
    """
    Get a list of log files in the specified directory.

    Parameters:
    dir (str): The directory path to search for log files.

    Returns:
    list: A list of log file paths found in the directory.
    """
    log_files = []
    for file in os.listdir(dir):
        if file.endswith(".log"):
            log_files.append(f'{dir}/{file}')
    return log_files

def parse_arguments():
    """
    Parses command-line arguments.

    Parameters:
    None

    Returns:
    argparse.Namespace: The parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Prints a CSV file from log files in a directory.")
    parser.add_argument("dir", type=str, help="The directory containing the log files.")
    return parser.parse_args()

def main(dir):
    """
    Prints a CSV file from log files in a directory.

    Parameters:
    dir (str): The directory containing the log files.

    Returns:
    None
    """
    print_header()
    log_files = get_log_files(dir)
    for log_file in log_files:
        print_row(log_file)

if __name__ == "__main__":
    args = parse_arguments()
    main(args.dir)

