import subprocess
import sys
from glob import glob

# This script create index for the bam files provided in the
# bam file folder.
# 
# Usage:
#         python index_bams.py /path/to/bams

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py /path/to/bam/files")
        exit(1)

    bam_folder = sys.argv[1]

    print("Indexing the files in the folder: {}".format( bam_folder ) )

    bam_files = glob( bam_folder + "/*bam" )

    for this_file in bam_files:
        print("Indexing", this_file)
        subprocess.run(["samtools", "index", this_file])

if __name__ == "__main__":
    main()