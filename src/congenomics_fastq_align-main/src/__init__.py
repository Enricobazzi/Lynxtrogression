import yaml
import os
from yaml.loader import FullLoader

class RunConfig:
    """
    A class to hold the configuration parameters for a run.
    
    Attributes:
    - sample_dict (dict): A dictionary mapping sample names to dictionaries mapping unique fastq_ids to the file names of the paired-end fastq files (r1 and r2).
    - samples (list): A list of sample names in the order they appear in the configuration file (within the dictionary sample_dict).
    - fastq_folder (str): The folder path where the fastq files are located.
    - reference_fasta (str): The path to the reference fasta file used for alignment.
    - output_folder (str): The folder path where the output files will be written.
    - threads (int): The number of threads to use for parallelized tools.
    - alignment_name (str): A label used to identify the alignment in the output file names.
    - call_bwa (str): The command or path to the executable for the BWA aligner.
    - call_samtools (str): The command or path to the executable for the SAMtools toolkit.
    - call_picard (str): The command or path to the executable for the Picard toolkit.
    - call_gatk (str): The command or path to the executable for the GATK toolkit.
    
    Methods:
    - __init__(self, config_file): Initializes a new RunConfig object from a YAML configuration file.
    - verify_arguments(self): Checks if the various elements of the configuration file exist
    """
    
    def __init__(self, config_file):
        """
        Initializes a new RunConfig object from a YAML configuration file.
        
        Parameters:
        - config_file (str): The path to the YAML configuration file.
        """
        with open(config_file) as file:
            data = yaml.load(file, Loader=FullLoader)
        
        self.sample_dict = data['sample_dict']
        self.samples = [sample for sample in list(self.sample_dict.keys())]
        self.reference_fasta = data['reference_fasta']
        self.output_folder = data['output_folder']
        self.threads = data['threads']
        self.alignment_name = data['alignment_name']
        self.call_bwa = data['call_bwa']
        self.call_samtools = data['call_samtools']
        self.call_picard = data['call_picard']
        self.call_gatk = data['call_gatk']        

    def verify_arguments(self):
        """
        Checks if the various elements of the configuration file exist
        """
        # check if path to fastqs exists
        if os.path.exists(self.reference_fasta) is False:
            print("reference genome {self.reference_fasta} not found!")
        else:
            print(f"reference genome found!")

            
def map_reads(fastq_id: str, sample: str, run: RunConfig) -> list:
    """
    Map reads from a pair of FASTQ files to a reference genome using BWA-MEM and SAMtools.
    
    Parameters:
    - fastq_id (str): The ID of the r1 and r2 FASTQ files pair to map reads from.
    - sample (str): The name of the sample to map reads for.
    - run (RunConfig): An object containing the parameters and file paths for the mapping run - created with the RunConfig() class.
    
    Returns:
    - list: A list of command-line arguments for running BWA and SAMtools to map reads.
    
    Example:
    >>> run = RunConfig(config/run.yml)
    >>> map_reads('id1', 'sample1', run)
    ['/path/to/bwa mem /path/to/reference.fasta /path/to/id1.fastq.gz /path/to/id1_2.fastq.gz -t 8 | /path/to/samtools view -hbS -@ 8 - -o /path/to/sample1_id1.bwa.bam']
    """
    
    # build BWA MEM + SAMTOOLS VIEW command
    bwa_command = [
        run.call_bwa, "mem", run.reference_fasta,
        f"{run.sample_dict[sample][fastq_id][0]}/{run.sample_dict[sample][fastq_id][1]}",
        f"{run.sample_dict[sample][fastq_id][0]}/{run.sample_dict[sample][fastq_id][2]}",
        "-t", str(run.threads),
        "|",
        run.call_samtools, "view", "-hbS", "-@", str(run.threads), "-",
        "-o", f"{run.output_folder}/{sample}_{fastq_id}_{run.alignment_name}.bam"
    ]
    
    return bwa_command


def sort_bam(bam_file: str, call_samtools: str, threads: int = 1) -> list:
    """
    Sort a BAM file using SAMtools sort command.
    
    Args:
    - bam_file (str): Path to the BAM file to be sorted.
    - call_samtools (str): Path to the SAMtools executable.
    - threads (int): Number of threads to use for sorting the BAM file (default=1).
    
    Returns:
    - list: A list of command-line arguments for running the SAMtools sort command to sort the input BAM file.
    
    Raises:
    - TypeError: If the BAM file does not end with '.bam'.

    Example:
    >>> sort_bam('/path/to/input.bam', '/path/to/samtools', 4)
    ['/path/to/samtools sort -@ 4 /path/to/input.bam -o /path/to/input_sorted.bam']
    """ 
    
    # Check that bam_file ends with '.bam'
    if not bam_file.endswith('.bam'):
        raise TypeError("bam file should end in .bam !")
    else:
        out_bam = bam_file[:-4] + '_sorted.bam'
    
    # Build SAMTOOLS SORT command
    samtools_sort = [
        call_samtools, "sort",
        "-@", str(threads),
        bam_file,
        "-o", out_bam
    ]
    return samtools_sort


def add_rg(bam_file: str, sample: str, fastq_id: str, call_picard: str) -> list:
    """
    Add read groups to a BAM file using PicardTools AddOrReplaceReadGroups.
    
    Args:
    - bam_file (str): Path to the input BAM file to add read groups to.
    - sample (str): Name of the sample to add read groups for.
    - fastq_id (str): ID of the FASTQ file to add read groups for.
    - call_picard (str): Path to the PicardTools executable to use.
    
    Returns:
    - list: A list of command-line arguments for running PicardTools AddOrReplaceReadGroups.

    Raises:
    - TypeError: If the BAM file does not end with '.bam'.
    
    Example:
    >>> call_picard = '/path/to/picard.jar'
    >>> bam_file = '/path/to/input.bam'
    >>> sample = 'sample1'
    >>> fastq_id = 'fastq1'
    >>> add_rg(bam_file, sample, fastq_id, call_picard)
    ['/path/to/picard.jar', 'AddOrReplaceReadGroups', 'I=/path/to/input.bam', 'O=/path/to/input_rg.bam', 'RGID=fastq1', 'RGLB=sample1_lib', 'RGPL=Illumina', 'RGPU=fastq1', 'RGSM=sample1', 'VALIDATION_STRINGENCY=SILENT']
    """
    
    # Check that bam_file ends with '.bam'
    if not bam_file.endswith('.bam'):
        raise TypeError("bam file should end in .bam !")
    else:
        out_bam = bam_file[:-4] + '_rg.bam'
    
    # Build PICARDTOOLS AddOrReplaceReadGroups command
    picard_add_rg = [
        call_picard, "AddOrReplaceReadGroups",
        f"I={bam_file}",
        f"O={out_bam}",
        f"RGID={fastq_id}",
        f"RGLB={sample}_lib",
        f"RGPL=Illumina",
        f"RGPU={fastq_id}",
        f"RGSM={sample}",
        f"VALIDATION_STRINGENCY=SILENT"
    ]
    return picard_add_rg


def make_rg_bamlist(sample: str, output_folder: str) -> list:
    """
    Constructs a command to generate a list of BAM files with read group (RG) information
    for a given sample and output folder.
    
    Args:
    - sample (str): The name of the sample to generate the BAM file list for.
    - output_folder (str): The folder where the BAM files are stored.
    
    Returns:
    - list: A list containing the ls command to generate the BAM file list.
    
    Examples:
    >>> make_rg_bamlist('sample1', '/path/to/output/folder')
    ['ls', '/path/to/output/folder/sample1_*_sorted_rg.bam', '>', '/path/to/output/folder/sample1.bam.list']
    """

    # Build ls command
    ls_bamlist = [
        "ls",
        f"{output_folder}/{sample}_*_sorted_rg.bam",
        ">",
        f"{output_folder}/{sample}.bam.list"
    ]

    return ls_bamlist


def merge_bams(bam_list: str, sample: str, call_samtools: str, alignment_name: str, output_folder: str, threads: int = 1) -> list:
    """
    Constructs a command to merge multiple BAM files into a single BAM file for a given sample.
    
    Args:
    - bam_list (str): The path to a file containing a list of BAM files to merge.
    - sample (str): The name of the sample to merge the BAM files for.
    - call_samtools (str): The command to call SAMTOOLS.
    - alignment_name (str): The name of the alignment
    - threads (int): The number of threads to use for the merging process (default=1).
    - output_folder (str): The folder where the merged BAM file will be stored.
    
    Returns:
    - list: A list containing the SAMTOOLS MERGE command to merge the BAM files.

    Examples:
        
    >>> merge_bams("/path/to/bam.list", "sample1", "samtools", 8, "/path/to/output/folder")
    ['samtools', 'merge', '-@', '8', '-r', '/path/to/output/folder/sample1_merged.bam', '-b', '/path/to/bam.list']
    """

    # Build SAMTOOLS MERGE command
    samtools_merge = [
        call_samtools, "merge",
        "-@", str(threads),
        "-r", f"{output_folder}/{sample}_{alignment_name}_sorted_rg_merged.bam",
        "-b", bam_list
    ]

    return samtools_merge

    
def rename_to_merged_bam(rg_bam: str, sample: str, alignment_name: str, output_folder: str) -> list:
    """
    Renames the given RG BAM file to a merged and sorted BAM file name as if it was the result of merging and resorting of a bamlist.

    Args:
    - rg_bam (str): The path and name of the input RG BAM file.
    - sample (str): The name of the sample.
    - alignment_name (str): The name of the alignment.
    - output_folder (str): The path of the output folder.

    Returns:
    - list: The command to rename the file.

    Examples:
    >>> rename_to_merged_bam("/path/to/sample1_alignment1_sorted_rg.bam", "sample1", "alignment1", "/path/to/output_folder")
    ["mv", "/path/to/sample1_alignment1_sorted_rg.bam", "/path/to/output_folder/sample1_alignment1_merged_sorted.bam"]
    """

    # Build mv command to rename the file
    mv_command = [
        "mv",
        rg_bam,
        f"{output_folder}/{sample}_{alignment_name}_sorted_rg_merged_sorted.bam"
    ]
    return mv_command


def mark_dups(bam_file: str, sample: str, call_picard: str, output_folder: str) -> list:
    """
    Marks duplicates in the input BAM file using Picard MarkDuplicates.

    Args:
    - bam_file (str): The path and name of the input BAM file.
    - sample (str): The name of the sample.
    - call_picard (str): The command to call Picard.
    - output_folder (str): The path of the output folder.

    Returns:
    - list: The command to run Picard MarkDuplicates.

    Raises:
    - TypeError: If the BAM file does not end with '.bam'.

    Examples:
        
    >>> mark_dups("/path/to/sample1.bam", "sample1", "/path/to/picard.jar", "/path/to/output_folder")
    ["/path/to/picard.jar", "MarkDuplicates", "METRICS_FILE=/path/to/output_folder/sample1_rmdup.txt", "I=/path/to/sample1.bam", "O=/path/to/output_folder/sample1_rmdup.bam", "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800"]
    """

    # Check that bam_file ends with '.bam'
    if not bam_file.endswith('.bam'):
        raise TypeError("bam file should end in .bam!")
    else:
        out_bam = bam_file[:-4] + '_rmdup.bam'

    # Build Picard MarkDuplicates command
    picard_command = [
        call_picard, "MarkDuplicates",
        f"METRICS_FILE={output_folder}/{sample}_rmdup.txt",
        f"I={bam_file}",
        f"O={out_bam}",
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800"
    ]

    return picard_command


def index_bam(bam_file: str, call_samtools: str) -> list:
    """
    Indexes the input BAM file using SAMtools index.

    Args:
    - bam_file (str): The path and name of the input BAM file.
    - call_samtools (str): The command to call SAMtools.

    Returns:
    - list: The command to run SAMtools index.

    Examples:
    >>> index_bam("/path/to/sample1.bam", "samtools")
    ["samtools", "index", "/path/to/sample1.bam"]
    """

    # Build SAMtools index command
    samtools_index = [
        call_samtools, "index",
        bam_file
    ]

    return samtools_index


def realigner_target_creator(bam_file: str, sample: str, call_gatk: str, reference_fasta: str, output_folder: str, threads: int = 1) -> list:
    """
    Creates a GATK RealignerTargetCreator command to generate a list of intervals to be used for local realignment of
    reads.

    Args:
    - bam_file (str): Path to input BAM file
    - sample (str): Sample name
    - call_gatk (str): Path to GATK executable
    - reference_fasta (str): Path to reference FASTA file
    - output_folder (str): Path to output folder
    - threads (int): Number of threads to use (default=1)
        
    Returns:
    - list: A list containing the GATK RealignerTargetCreator command and its arguments.

    Example:
    >>> bam_file = "/path/to/input.bam"
    >>> sample = "my_sample"
    >>> output_folder = "/path/to/output/folder"
    >>> threads = 4
    >>> reference_fasta = "/path/to/reference.fasta"
    >>> gatk_path = "/path/to/gatk.jar"
    >>> gatk_command = realigner_target_creator(bam_file, sample, output_folder, threads, reference_fasta, gatk_path)
    >>> print(gatk_command)
    ['java', '-jar', '/path/to/gatk.jar', '-T', 'RealignerTargetCreator', '-nt', '4', '-R', '/path/to/reference.fasta', '-I', '/path/to/input.bam', '-o', '/path/to/output/folder/my_sample_realignertargetcreator.intervals']
    """

    # build GATK RealignerTargetCreator command
    gatk_command = [
        call_gatk, "-T", "RealignerTargetCreator",
        "-nt", str(threads),
        "-R", reference_fasta,
        "-I", bam_file,
        "-o", f"{output_folder}/{sample}_realignertargetcreator.intervals"
    ]

    return gatk_command


def indel_realigner(bam_file: str, intervals: str, call_gatk: str, reference_fasta: str) -> list:
    """
    Realigns indels in the specified BAM file using GATK's IndelRealigner.

    Args:
    - bam_file (str): The BAM file to realign.
    - intervals (str): The intervals file produced by GATK's RealignerTargetCreator.
    - call_gatk (str): The path to the GATK executable.
    - reference_fasta (str): The path to the reference FASTA file.

    Returns:
    - list: A list representing the GATK IndelRealigner command.

    Example:
    >>> bam_file = "/path/to/input.bam"
    >>> intervals_file = "/path/to/intervals.intervals"
    >>> call_gatk = "/path/to/gatk"
    >>> reference_fasta = "/path/to/reference.fasta"
    >>> indel_realigner(bam_file, intervals_file, call_gatk, reference_fasta)
    ["/path/to/gatk", "-T", "IndelRealigner", "-R", "/path/to/reference.fasta", 
    "-targetIntervals", "/path/to/intervals.intervals", "-I", "/path/to/input.bam", 
    "-o", "/path/to/input_indelrealigner.bam"]
    """
    # Check that bam_file ends with '.bam'
    if not bam_file.endswith('.bam'):
        raise TypeError("bam file should end in .bam!")
    else:
        out_bam = bam_file[:-4] + '_indelrealigner.bam'

    # Build GATK IndelRealigner command
    gatk_command = [
        call_gatk, "-T", "IndelRealigner",
        "-R", reference_fasta,
        "-targetIntervals", intervals,
        "-I", bam_file,
        "-o", out_bam
    ]

    return gatk_command


def print_sample_pipeline(sample: str, run: object) -> list:
    """
    Print a pipeline of commands to align, sort, add read groups, merge bams, mark duplicates,
    create realignment targets, and realign indels for a given sample.
    
    Args:
    - sample (str): Name of the sample to be processed.
    - run (object): Object containing various pipeline parameters.
    
    Returns:
    - None
    
    Example:
    print_sample_pipeline(sample="sample1", run=run)
    """
    bash_script = ''
    
    bash_script = bash_script + f"# Align reads of {sample}\n"
    
    fastq_ids = list(run.sample_dict[sample].keys())
    
    for fastq_id in fastq_ids:
        
        bash_script = bash_script + f"\n# align reads from files '{run.sample_dict[sample][fastq_id][1]}' and '{run.sample_dict[sample][fastq_id][2]}':\n"
        map_reads_command = map_reads(sample = sample, fastq_id = fastq_id, run = run)
        bash_script = bash_script + ' '.join(map_reads_command)
        first_bam = map_reads_command[-1]
        bash_script = bash_script + f"\n# -> {first_bam}\n"
        
        bash_script = bash_script + f"\n# sort the resulting bam file:\n"
        sort_first_bam = sort_bam(bam_file = first_bam, call_samtools = run.call_samtools, threads = run.threads)
        bash_script = bash_script + ' '.join(sort_first_bam)
        sorted_bam = sort_first_bam[-1]
        bash_script = bash_script + f"\n# -> {sorted_bam}\n"

        bash_script = bash_script + f"\n# add read groups to the sorted bam file:\n"
        add_rgs_to_sorted = add_rg(bam_file = sorted_bam, sample = sample, fastq_id = fastq_id, call_picard = run.call_picard)
        bash_script = bash_script + ' '.join(add_rgs_to_sorted)
        rg_bam = add_rgs_to_sorted[3][2:]
        bash_script = bash_script + f"\n# -> {rg_bam}\n"
    
    if len(fastq_ids) > 1:
        bash_script = bash_script + f"\n# {sample}'s bams should be merged!\n"
        
        bash_script = bash_script + f"\n# building {sample}'s list of bams:\n"
        make_bamlist = make_rg_bamlist(sample = sample, output_folder = run.output_folder)
        bash_script = bash_script + ' '.join(make_bamlist)
        bam_list = make_bamlist[-1]
        bash_script = bash_script + f"\n# -> {bam_list}\n"
        
        bash_script = bash_script + f"\n# merging {sample}'s bams:\n"
        merge_sample_bams = merge_bams(bam_list = bam_list, sample = sample, call_samtools = run.call_samtools,
                                   alignment_name = run.alignment_name, threads = run.threads, output_folder = run.output_folder)
        bash_script = bash_script + ' '.join(merge_sample_bams)
        merged_bam = merge_sample_bams[5]
        bash_script = bash_script + f"\n# -> {merged_bam}\n"
        
        bash_script = bash_script + f"\n# sorting {sample}'s merged bam:\n"
        sort_merged_bam = sort_bam(bam_file = merged_bam, call_samtools = run.call_samtools, threads = run.threads)
        bash_script = bash_script + ' '.join(sort_merged_bam)
        merged_bam = sort_merged_bam[-1]
        bash_script = bash_script + f"\n# -> {merged_bam}\n"
        
    else:
        bash_script = bash_script + f"\n# {sample} only has one fastqid!\n"
        bash_script = bash_script + f"# renaming its only sorted_rg BAM to merged_sorted:\n"
        rename = rename_to_merged_bam(rg_bam = rg_bam, sample = sample, alignment_name = run.alignment_name, output_folder = run.output_folder)
        bash_script = bash_script + ' '.join(rename)
        merged_bam = rename[-1]
        bash_script = bash_script + f"\n# -> {merged_bam}\n"
    
    bash_script = bash_script + f"\n# marking duplicates of merged bam:\n"
    md_command = mark_dups(bam_file = merged_bam, sample = sample, call_picard = run.call_picard, output_folder = run.output_folder)
    bash_script = bash_script + ' '.join(md_command)
    md_bam = md_command[4][2:]
    bash_script = bash_script + f"\n# -> {md_bam}\n"
    
    bash_script = bash_script + f"\n# index the duplicate marked bam for gatk:\n"
    index_cmd = index_bam(bam_file = md_bam, call_samtools = run.call_samtools)
    bash_script = bash_script + ' '.join(index_cmd)
    bash_script = bash_script + f"\n# -> {md_bam}.bai\n"
    
    bash_script = bash_script + f"\n# identify indel realignment targets:\n"
    real_tar = realigner_target_creator(bam_file = md_bam, sample = sample, call_gatk = run.call_gatk, reference_fasta = run.reference_fasta, output_folder = run.output_folder, threads = run.threads)
    bash_script = bash_script + ' '.join(real_tar)
    target_file = real_tar[-1]
    bash_script = bash_script + f"\n# -> {target_file}\n"
    
    bash_script = bash_script + f"\n# realign indels:\n"
    ind_real = indel_realigner(bam_file = md_bam, intervals = target_file, call_gatk = run.call_gatk, reference_fasta = run.reference_fasta)
    bash_script = bash_script + ' '.join(ind_real)
    real_bam = ind_real[-1]
    bash_script = bash_script + f"\n# -> {real_bam}\n"
    
    bash_script = bash_script + f"\n# index the indel realigned bam for future analyses:\n"
    index_cmd_2 = index_bam(bam_file = real_bam, call_samtools = run.call_samtools)
    bash_script = bash_script + ' '.join(index_cmd_2)
    bash_script = bash_script + f"\n# -> {md_bam}.bai\n"
        
    return bash_script