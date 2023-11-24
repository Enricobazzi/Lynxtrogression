"""
This script will be used to plot the fitness values of a run across its generations.

The fitness values will be extracted from the log file of the run.

Two options are available:
    - Plot the fitness values of a single run.
    - Plot the fitness values of multiple runs in the same folder.

The script can be run from the command line as follows:
    python plot_fitness_of_logs.py -l <log> -o <output_file> -x <x_min,x_max> -y <y_min,y_max> -g <n_gen>

and the script will either plot the fitness of the single run if log (-l) is a single file, 
or plot the fitness of all the runs in the folder if log (-l) is a folder.

The arguments -x, -y and -g are optional.

The output file will be a pdf file containing the plot.
"""

import os
import matplotlib.pyplot as plt
import argparse

def get_distances_list(log_file):
    distances_list = []
    with open(log_file, "r") as file:
        for line in file:
            if 'GENERATION' in line:
                distances_list.append(float(line.split(' ')[2]))
    return distances_list

def plot_distances(log_file, dist_list, x_lim = None, y_lim = None, n_gen = 60):
    # Generate x-axis values (assuming y-values are given in the list and x-values are indices)
    x_values = range(len(dist_list))
    leg = log_file.split('/')[-1].split('.')[0].replace('run_', '')
    # Create a line plot
    if n_gen == None:
        plt.plot(x_values, dist_list, label=leg)
    else:
        plt.plot(x_values[-n_gen:], dist_list[-n_gen:], label=leg)
    
    # Add labels and title
    plt.xlabel('generation')
    plt.ylabel('distance')
    # plt.title(log_file)
    # plt.legend()
    if x_lim:
        plt.xlim(x_lim[0], x_lim[1])
    if y_lim:
        plt.ylim(y_lim[0], y_lim[1])


def plot_distances_from_folder(log_folder, x_lim = None, y_lim = None, n_gen = 60):
    # Get all log files in the folder
    log_files = [os.path.join(log_folder, f) for f in os.listdir(log_folder) if f.endswith('.log')]
    # Get the distances list for each log file
    dist_lists = [get_distances_list(log_file) for log_file in log_files]
    # Plot the distances
    for log_file, dist_list in zip(log_files, dist_lists):
        plot_distances(log_file, dist_list, x_lim, y_lim, n_gen)
    # Add labels and title
    plt.xlabel('generation')
    plt.ylabel('distance')
    # plt.title(log_folder)
    # plt.legend(bbox_to_anchor=(0., 1.02, .5, .02), loc=0,
    #       ncol=1, mode="expand", borderaxespad=0.)
    # plt.legend()
    if x_lim:
        plt.xlim(x_lim[0], x_lim[1])
    if y_lim:
        plt.ylim(y_lim[0], y_lim[1])


def parse_arguments():
    parser = argparse.ArgumentParser(description='Plot the fitness values of a run or multiple runs.')
    parser.add_argument('-l', '--log', required=True, type=str, help='The log file or folder.')
    parser.add_argument('-o', '--output_file', required=False, type=str, help='The output file.')
    parser.add_argument('-x', '--x_lim', required=False, type=str, help='The x-axis limits.')
    parser.add_argument('-y', '--y_lim', required=False, type=str, help='The y-axis limits.')
    parser.add_argument('-g', '--n_gen', required=False, type=int, help='The number of generations to plot.')
    return parser.parse_args()

def main():
    args = parse_arguments()
    log = args.log
    output_file = args.output_file
    x_lim = args.x_lim
    y_lim = args.y_lim
    n_gen = args.n_gen

    if x_lim:
        x_lim = [float(x) for x in x_lim.split(',')]
    if y_lim:
        y_lim = [float(y) for y in y_lim.split(',')]

    if not n_gen:
        n_gen = 60

    if os.path.isdir(log):
        plot_distances_from_folder(log, x_lim, y_lim, n_gen)
    else:
        dist_list = get_distances_list(log)
        plot_distances(log, dist_list, x_lim, y_lim, n_gen)

    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

if __name__ == "__main__":
    main()

