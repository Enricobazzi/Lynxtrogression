
"""
In this script, we will make a plink range file from the masked regions file.

The masked regions file is a txt file with the following format:
chr1 1000 2000 : 1100 1200, 1400 1500, 1900 2000
chr1 3000 4000 : 3300 3650, 3900 3950

The output plink range file is a txt file with the following format:
chr1 1100 1200 1000-2000
chr1 1400 1500 1000-2000
chr1 1900 2000 1000-2000
chr1 3300 3650 3000-4000
chr1 3900 3950 3000-4000
"""

import argparse

# chromosome dictionary
chr_to_scafold_conversion = {}
chr_to_scafold_conversion["1"] = "Super_Scaffold_1"
chr_to_scafold_conversion["2"] = "Super_Scaffold_11"
chr_to_scafold_conversion["3"] = "Super_Scaffold_12"
chr_to_scafold_conversion["4"] = "Super_Scaffold_13"
chr_to_scafold_conversion["5"] = "Super_Scaffold_14"
chr_to_scafold_conversion["6"] = "Super_Scaffold_2"
chr_to_scafold_conversion["7"] = "Super_Scaffold_3"
chr_to_scafold_conversion["8"] = "Super_Scaffold_4"
chr_to_scafold_conversion["9"] = "Super_Scaffold_5"
chr_to_scafold_conversion["10"] = "Super_Scaffold_6"
chr_to_scafold_conversion["11"] = "Super_Scaffold_7"
chr_to_scafold_conversion["12"] = "Super_Scaffold_8"
chr_to_scafold_conversion["13"] = "Super_Scaffold_9"
chr_to_scafold_conversion["14"] = "scaffold_11_arrow_ctg1"
chr_to_scafold_conversion["15"] = "scaffold_17_arrow_ctg1"
chr_to_scafold_conversion["16"] = "scaffold_18_arrow_ctg1"
chr_to_scafold_conversion["17"] = "scaffold_21_arrow_ctg1"
chr_to_scafold_conversion["18"] = "scaffold_2_arrow_ctg1"

def main():
    # read the arguments
    parser = argparse.ArgumentParser(description="Make plink range file from masked regions file.")
    parser.add_argument("--masked_regions_file", help="Path to the masked regions file.")

    args = parser.parse_args()

    # read the masked regions file
    with open(args.masked_regions_file, "r") as f:
        masked_regions = f.readlines()

    # plink range file (remove .txt extension and add .plink_range extension)
    plink_range_file = args.masked_regions_file[:-4] + ".plink_range"

    # write the plink range file
    with open(plink_range_file, "w") as f:
        for line in masked_regions:
            line = line.strip().split()
            chrom = chr_to_scafold_conversion[line[0]]
            start = line[1]
            end = line[2]
            for i in range(4, len(line), 2):
                masked_region_start = line[i]
                masked_region_end = line[i+1].replace(",", "")
                f.write(chrom + " " + masked_region_start + " " + masked_region_end + " " + start + "-" + end + "\n")
    
if __name__ == "__main__":
    main()




