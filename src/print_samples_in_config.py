"""
This script prints sample names in the config file.
The config file is a yaml file with the following structure:
    sample_dict:
      sample1:
        fastq_pair_1:
          ['path/to/pair1', 'pair1_r1.fastq.gz', 'pair1_r2.fastq.gz']
        fastq_pair_2:
          ['path/to/pair2', 'pair2_r1.fastq.gz', 'pair2_r2.fastq.gz']
        sample2:
          fastq_pair_1:
            ['path/to/pair3', 'pair3_r1.fastq.gz', 'pair3_r2.fastq.gz']
Usage:
    python print_samples_in_config.py <config_file>
Example:
    python print_samples_in_config.py config/all_fastp_alignment.yaml
"""

import sys
import yaml

def main():
    # Load the YAML file
    with open(sys.argv[1], 'r') as f:
        data = yaml.safe_load(f)
    
    # Loop through each sample in the YAML file
    for sample in data['sample_dict']:
        print(sample)

if __name__ == '__main__':
    main()

