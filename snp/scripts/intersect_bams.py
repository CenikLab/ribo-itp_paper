import subprocess
import sys
from glob import glob
import os

# This script create index for the bam files provided in the
# bam file folder.
# 
# Usage:
#         python intersect_bams.py vcf_file /path/to/bamfolder output_folder

def main():
    if len(sys.argv) < 3:
        print("Usage: python intersect_bams.py vcf_file /path/to/bamfolder output_folder")
        exit(1)

    vcf_file      = sys.argv[1]
    bam_folder    = sys.argv[2]
    output_folder = sys.argv[3]

    os.makedirs(output_folder, exist_ok=True)

    bam_files = glob( bam_folder + "/*bam" )

    for this_file in bam_files:
        print("Intersecting", this_file)

        file_name = os.path.basename(this_file)
        output_file = os.path.join(output_folder, file_name)

        output_stream = open(output_file, "wb")

        p1 = subprocess.run(["bedtools", "intersect", 
                        "-a",  this_file, "-b", vcf_file, 
                        "-wa"], stdout=output_stream)

        output_stream.close()



if __name__ == "__main__":
    main()