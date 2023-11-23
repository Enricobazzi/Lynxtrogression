
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
            chrom = line[0]
            start = line[1]
            end = line[2]
            for i in range(4, len(line), 2):
                masked_region_start = line[i]
                masked_region_end = line[i+1].replace(",", "")
                f.write(chrom + " " + masked_region_start + " " + masked_region_end + " " + start + "-" + end + "\n")
    
if __name__ == "__main__":
    main()




