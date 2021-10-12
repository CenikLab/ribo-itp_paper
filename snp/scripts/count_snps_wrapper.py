import subprocess
import sys
from glob import glob
import argparse
import os

def get_parameters():
    parser = argparse.ArgumentParser(
                  description = "Wrapper for count_snps.py")

    parser.add_argument("-o", "--out",
                        type     = str,
                        required = True,
                        help     = "Output folder")

    parser.add_argument("-i", "--input",
                        type     = str,
                        required = True,
                        help     = "Input folder containing bam files")

    parser.add_argument("--vcf",
                        type     = str,
                        required = True,
                        help     = "vcf file containing the snps")


 
    args = parser.parse_args()
    return args

####################################################################

def main():

    arguments  = get_parameters()

    out_folder = arguments.out

    bam_files  = glob( arguments.input + "/*bam" )

    script_folder = os.path.dirname(os.path.realpath(__file__))

    script = os.path.join(script_folder, "count_snps.py")

    os.makedirs(out_folder, exist_ok=True)

    for this_file in bam_files:
        print("Processing", this_file)
        experiment_name = os.path.basename(this_file).split(".")[0]
        output_file = os.path.join(out_folder, experiment_name + ".tsv.gz")
        run_tuple = ["python", script, "-o", output_file, "--bam", this_file, "--vcf", arguments.vcf]
        subprocess.run(run_tuple)

if __name__ == "__main__":
    main()