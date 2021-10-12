
import gzip
import argparse
import re

import sys
import os

import numpy as np


script_folder = os.path.dirname(os.path.realpath(__file__))


sys.path.insert(0, os.path.dirname(script_folder) )

import ref_lib
from ref_lib.Vcf import VcfFile

from collections import defaultdict, OrderedDict

import pysam

##################################################################

def get_parameters():
    parser = argparse.ArgumentParser(
                  description = "Reports the count of each nucleotide at each "
                                "nucleotide position. "
                                "Reads with INDELS are discarded.")

    parser.add_argument("-o", "--out",
                        type     = str,
                        required = True,
                        help     = "Output file")

    parser.add_argument("--bam",
                        type     = str,
                        required = True,
                        help     = "Bam file coming from transcriptomic alignments")

    parser.add_argument("--vcf",
                        type     = str,
                        required = True,
                        help     = "vcf file containing the snps")


 
    args = parser.parse_args()
    return args

##############################################################################

def read_vcf(vcf_file):
    vcf_reader = VcfFile(vcf_file)

    contents_dict = OrderedDict()

    for entry in vcf_reader:

        transcript_id = entry.fields["CHROM"]
        position      = int(entry.fields["POS"])
        ref           = entry.fields["REF"]
        alt           = entry.fields["ALT"]

        if not contents_dict.get(transcript_id):
            contents_dict[transcript_id] = OrderedDict()

        contents_dict[transcript_id][position] = { "REF": ref, "ALT": alt,
                                                   "A" : 0, "C" : 0,
                                                   "G" : 0, "T" : 0 , "N" : 0}

    return contents_dict


##############################################################################
        
def count_nucleotides(vcf_contents, bam_file):
    """
    Note that psam coordinates are 0-based and vcf coordinates are 1-based
    So we need to make a conversion
    """
    pass
    
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    for transcript_header, pos_dict in vcf_contents.items():
        iter = samfile.fetch(transcript_header)
        sorted_positions = sorted(list(pos_dict.keys() ) )

        for alignment in iter:
            #print(alignment)

            #print("###")
            
            # Coversion from 0-based to 1-based
            read_start_position = alignment.reference_start + 1
            #print(read_start_position)
            #print("###")
            # Coversion from 0-based to 1-based
            # Note that read end position is EXCLUDED
            # It doesn't belong to the read
            # Read end position is the position 1 after the last nuc of the read
            # This comes from pysam convention
            # also simplifies our search
            read_end_position = alignment.reference_end + 1
            #print(read_end_position)

            snp_start_index = np.searchsorted(sorted_positions, read_start_position)
            snp_stop_index  = np.searchsorted(sorted_positions, read_end_position)

            # If there are any insertions or deletions skip

            skip_read = False

            for code, count in alignment.cigartuples:
                if code != 0:
                    skip_read = True
                    # print(alignment.cigartuples)

            if skip_read:
                continue


            for i in range( snp_start_index, snp_stop_index ):
                #print( sorted_positions[i] )

                nucleotide_in_read = alignment.query_sequence[ sorted_positions[i] - read_start_position ]
            
                pos_dict[ sorted_positions[i] ][nucleotide_in_read] += 1


            # print("------------------------")

    return vcf_contents

###########################################################################

def write_nuc_counts_dict(count_dict, output_file):

    header_array = ["transcript", "position", "REF", "ALT", "A", "C", "G", "T"]
     
    if output_file.lower().endswith("gz"):
        myopen = gzip.open
    else:
        myopen = open

    with myopen(output_file, "wt") as output_stream:
        print("\t".join(header_array), file = output_stream)

        for transcript_header, pos_dict in count_dict.items():
            for position in sorted(pos_dict.keys()):
                row_array = [transcript_header, str(position)] + \
                            [ pos_dict[position][x] for x in header_array[2:] ]

                row_array = map(str, row_array )            

                print( "\t".join(row_array), file = output_stream)



###########################################################################

def main():
    arguments = get_parameters()

    vcf_contents = read_vcf(arguments.vcf)

    nucleotide_counts_dict = count_nucleotides(vcf_contents, arguments.bam)

    write_nuc_counts_dict(count_dict  = nucleotide_counts_dict, 
                          output_file = arguments.out)


if __name__ == "__main__":
    main()
