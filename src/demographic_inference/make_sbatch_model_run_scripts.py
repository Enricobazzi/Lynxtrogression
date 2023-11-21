"""
This script generates a bash script for a specific model and run number.

The bash script will run the model on the cluster.

It takes the run number and the model's configuration file as first two arguments.

The third argument is the path to where the script should be saved.

The configuration file has the following format:
    working_folder=PATH_TO_WORKING_FOLDER
    bed_file=BED_FILE_PREFIX
    individuals_training=SAMPLE1,SAMPLE2
    individuals_replication=SAMPLE3,SAMPLE4
    fsc_name=FSC_NAME
    fsc_folder=FSC_FOLDER
    model_to_run=MODEL_NUMBER
    model_name=MODEL_NAME
    scripts_folder=PATH_TO_SCRIPTS_FOLDER
    logs_folder=PATH_TO_LOGS_FOLDER
    jar_file=PATH_TO_JAR_FILE

Example usage:
python make_sbatch_model_run_scripts.py 1 config/model_1.config.properties

This will generate a bash script called run_model_1_run_1.sh

The bash script can be run on the cluster with:
sbatch run_model_1_run_1.sh
"""

import sys
import os
import textwrap


def main():
    run_number = int(sys.argv[1])
    config_file = sys.argv[2]
    
    # Read the configuration file
    config = {}
    with open(config_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            key, value = line.strip().split('=')
            config[key] = value

    # Get the model name
    model_name = config['model_name']

    # Get the working folder
    working_folder = config['working_folder']

    # Get the jar file
    jar_file = config['jar_file']

    # Get the logs folder
    logs_folder = config['logs_folder']

    script  = textwrap.dedent(
        f'''\
            #!/bin/bash
            #SBATCH --job-name={model_name}_run_{run_number}
            #SBATCH --output={logs_folder}/run_{model_name}_{run_number}.out
            #SBATCH --error={logs_folder}/run_{model_name}_{run_number}.log
            #SBATCH --time=3-00:00:00
            #SBATCH --mem=20G
            #SBATCH --cpus-per-task=4
            
            module load cesga/2020 jdk/17.0.2
            
            cd {working_folder}
            
            # copy the template run folder to a new folder for this run
            cp -r fastSimcoal2_template fastSimcoal2_{run_number}
            
            # run lynx_ea_abc
            java -jar {jar_file} {run_number} ../../../{config_file}
            '''
    )

    # Write the script to a file in the scripts folder
    scripts_folder = sys.argv[3]
    script_file = os.path.join(scripts_folder, f'run_{model_name}_{run_number}.sh')

    with open(script_file, 'w') as f:
        f.write(script)

    
if __name__ == '__main__':
    main()
