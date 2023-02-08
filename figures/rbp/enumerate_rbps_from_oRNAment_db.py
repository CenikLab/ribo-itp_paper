#!/usr/bin/env python3
"""Runs enumeration of RBPs from oRNAment db files."""

import argparse
import collections
import gzip
import os
import sys

DB_SEP = ","

DB_HEADER = ("ensembl_gene_id", "ensembl_transcript_id", "gene_biotype", "transcript_biotype",
             "transcript_position", "RBP", "score", "unpaired_probability",
             "chromosome", "region", "exon_start", "exon_end")

DB_INDEX = {e: i for (i, e) in enumerate(DB_HEADER)}

IDS_HEADER = ("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id_INT",
              "ensembl_gene_id_INT")

IDS_INDEX = {e: i for (i, e) in enumerate(IDS_HEADER)}

RBP_IDS_HEADER = ("RBP", "RBP_name")
RBP_IDS_INDEX = {e: i for (i, e) in enumerate(RBP_IDS_HEADER)}

ANNOT_TUPLE = collections.namedtuple("ANNOT_TUPLE", "trx_id, gene_id, gene_name")

DEFAULT_SCORE = 1.0
DEFAULT_OUTFILE = "oRNAment_db_res.txt"
UTR5_REGIONS = {"5;5", "5;C"}
UTR3_REGIONS = {"3;3", "C;3"}
CDS_REGIONS = {"C;C"}
UTR5_IDX = 0
CDS_IDX = 1
UTR3_IDX = 2
OUT_HEADER = ("Transcript_ID", "RBP", "UTR5_counts", "CDS_counts", "UTR3_counts")


def parse_commandline_params(args):
    """Parses command line parameters.

    :param list args: command line arguments, no script name
    :return argparse.Namespace: namespace object with dict-like access
    """

    parser = argparse.ArgumentParser(description="%s arguments" % __name__)

    parser.add_argument("-i", "--trx_int_ids", type=str, required=True,
                        help='oRNAment integer transcript IDs to select for enumeration, one per line.')

    parser.add_argument("-d", "--oRNAment_db_gz", type=str, required=True,
                        help='Gzipped oRNAment database file.')

    parser.add_argument("-t", "--oRNAment_trx_id_encoding_gz", type=str, required=True,
                        help='Gzipped oRNAment transcript ID encodings.')

    parser.add_argument("-r", "--oRNAment_rbp_id_encoding", type=str, required=True,
                        help='oRNAment RBP ID encodings.')

    parser.add_argument("-s", "--score", type=float, default=DEFAULT_SCORE, 
                        help='Score (matrix similarity score) for filtering RBP motifs.')

    parser.add_argument("-o", "--output_dir", type=str, required=False, default=".",
                        help='Optional output directory. Default current working directory.')

    parser.add_argument("-f", "--filename", type=str, required=False, default=DEFAULT_OUTFILE,
                        help='Optional filename to write to output_dir. Default %s.' % DEFAULT_OUTFILE)

    parsed_args = vars(parser.parse_args(args))
    return parsed_args


def main():
    """Runs the workflow when called from command line."""

    parsed_args = parse_commandline_params(sys.argv[1:])

    outdir = parsed_args["output_dir"]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Read in the ID encodings
    with gzip.open(parsed_args["oRNAment_trx_id_encoding_gz"], "rt") as trx_encoding_fh:
        trx_id_encoding = {}
        for line in trx_encoding_fh:
            line_split = line.strip().split(DB_SEP)
            trx_id_int = line_split[IDS_INDEX["ensembl_transcript_id_INT"]]
            trx_id = line_split[IDS_INDEX["ensembl_transcript_id"]]
            gene_id = line_split[IDS_INDEX["ensembl_gene_id"]]
            gene_name = line_split[IDS_INDEX["external_gene_name"]]
            trx_id_encoding[trx_id_int] = ANNOT_TUPLE(trx_id, gene_id, gene_name)

    with open(parsed_args["oRNAment_rbp_id_encoding"], "r") as trx_encoding_fh:
        rbp_id_encoding = {}
        for line in trx_encoding_fh:
            line_split = line.split(DB_SEP)
            rbp_id_int = int(line_split[RBP_IDS_INDEX["RBP"]])
            rbp = line_split[RBP_IDS_INDEX["RBP_name"]].lstrip(' "').rstrip('"\n')
            rbp_id_encoding[rbp_id_int] = rbp

    # Get the transcripts of interest
    with open(parsed_args["trx_int_ids"], "r") as trx_int_ids_fh:
        trx_int_ids = {line.strip() for line in trx_int_ids_fh}

    # Now enumerate RBPs using a dict keyed by trx_id_int, rbp
    rbp_counts = collections.OrderedDict()
    with gzip.open(parsed_args["oRNAment_db_gz"], "rt") as ornament_db_fh:

        for line in ornament_db_fh:

            line_split = line.strip().split(DB_SEP)
            trx_id_int = line_split[DB_INDEX["ensembl_transcript_id"]].strip(".")

            # For some reason versions are include in the oRNAment db trx INT ids
            # gzcat Mus_musculus_cDNA_oRNAment.csv.gz | cut -d ',' -f2 | fgrep -c . -

            # But they are not in the encodings file:
            # gzcat Mus_musculus_string_to_int_ID_conversion.csv.gz | cut -d ',' -f4 | head
            trx_id_int_noversion = trx_id_int.split(".")[0]

            if trx_id_int_noversion not in trx_int_ids:
                continue

            rbp_id_int = line_split[DB_INDEX["RBP"]]

            score = float(line_split[DB_INDEX["score"]])
            if score < parsed_args["score"]:
                continue

            key = (trx_id_int_noversion, int(rbp_id_int))
            if key not in rbp_counts:
                rbp_counts[key] = [0, 0, 0]

            region = line_split[DB_INDEX["region"]]
            if region in UTR5_REGIONS:
                rbp_counts[key][UTR5_IDX] += 1
            elif region in UTR3_REGIONS:
                rbp_counts[key][UTR3_IDX] += 1
            else:
                rbp_counts[key][CDS_IDX] += 1

    # Write the results
    outfile = os.path.join(outdir, parsed_args["filename"])
    with open(outfile, "w") as out_fh:
        out_fh.write("\t".join(OUT_HEADER) + "\n")
        for k, v in rbp_counts.items():
            trx_id_int, rbp_id_int = k
            rbp_name = rbp_id_encoding[rbp_id_int]
            trx_id = trx_id_encoding[trx_id_int].trx_id
            utr5_count = v[UTR5_IDX]
            cds_count = v[CDS_IDX]
            utr3_count = v[UTR3_IDX]
            out_fh.write("\t".join(list(map(str, (trx_id, rbp_name, utr5_count, cds_count, utr3_count)))) + "\n")


if __name__ == "__main__":
    main()
