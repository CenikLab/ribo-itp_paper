# Ribo-ITP: Scripts, Pipelines, Notebooks and References

This repository contains the pipelines & their parameters, downstream processing & analysis scripts 
as well as the transcriptome reference and annotation files used for the manuscript **(link needed)**.

## Contents
  * [RiboFlow pipeline and parameter files](https://github.com/CenikLab/ribo-itp_paper/tree/main/Riboflow)
  * [RNA-Seq Pipeline](https://github.com/CenikLab/ribo-itp_paper/tree/main/rnaseq/pipeline)
  * [SNP calling](https://github.com/CenikLab/ribo-itp_paper/tree/main/snp)
  * [Analyses and figures](https://github.com/CenikLab/ribo-itp_paper/tree/main/figures)
  * [References (in another repository)](https://github.com/CenikLab/ribo-itp_paper_references) 

## Requirements

  1) **.ribo files:** To run the analysis scripts, ribo files, for human and mouse samples  **(link to GEO)**, need to be downloaded by the users.  

  2) **Sequencing Files:** To run the ribosome profiling pipeline, RiboFlow, and the RNA-Seq data, fastq files, available in GEO **(link needed)**, are required.
  
  
## [RiboFlow pipeline and parameter files](https://github.com/CenikLab/ribo-itp_paper/tree/main/Riboflow)

Ribosome profiling data were processed using a modified version of [RiboFlow](https://github.com/ribosomeprofiling/riboflow), which is available [here](https://github.com/CenikLab/ribo-itp_paper/tree/main/Riboflow).
This modified version of RiboFlow can handle ribosome profiling libraries with unique molecular identifiers (UMIs). 
We also provide parameter files for each RiboFlow run.

RNA-Seq data were **NOT** processed via RiboFlow. We used a custom pipeline to process RNA-Seq files whcih is available in this repository (see below).

We recommend running RiboFlow in a [conda](https://conda.io/projects/conda/en/latest/) environment. 
[Umi-tools](https://umi-tools.readthedocs.io/en/latest/) is needed to process ribosome profiling data.


## [RNA-Seq Pipeline](https://github.com/CenikLab/ribo-itp_paper/tree/main/rnaseq/pipeline)

We processed the RNA-seq data using our custom pipeline. Reads were mapped, deduplicated and quantified by [this](https://github.com/CenikLab/ribo-itp_paper/blob/main/rnaseq/pipeline/Snakefile) Snakemake pipeline. 


## [SNP calling](https://github.com/CenikLab/ribo-itp_paper/tree/main/snp)

Scripts and python notebooks, for calling single nucleotide polymorphisms (SNPs), are provided in this folder. [The list of SNPs](https://github.com/CenikLab/ribo-itp_paper/blob/main/snp/reference_files/cds_of_transcriptomic_variants.vcf.gz), 
on the coding sequence of the transcripts, in VCF format, is also provided.

We also provide a detailed list of counts for each SNP [in this file](https://github.com/CenikLab/ribo-itp_paper/blob/main/snp/notebooks/snp_dataframes/riboseq_detailed_snps.csv.gz) for ribosome profiling experiments. For RNA-Seq experiments, a [similar list](https://github.com/CenikLab/ribo-itp_paper/blob/main/snp/notebooks/snp_dataframes/rnaseq_experimentwise_snp_counts.csv) is provided.

## [Analyses and figures](https://github.com/CenikLab/ribo-itp_paper/tree/main/figures)

R scripts, used for differential ribosome occupancy, reproducibility and other figures and analyses are provided in this folder. Note that additional packages (for example [Enhanced Volcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html), [Seurat](https://satijalab.org/seurat/), etc.) are needed to run these R scripts. The paths to the .ribo files and count files are given relative to the scripts. The users might need to adjust these paths based on their active directory or file organization. 

These scripts heavily use the number of reads on the CDS of the transcripts. For RNA-Seq samples, this is given in a [csv file](https://github.com/CenikLab/ribo-itp_paper/blob/main/cds_counts.csv.gz). For ribosome profiling samples, .ribo files are used to extract the counts. Also, for convenience, we combined ribosome profiling and RNA-Seq counts into [one file](https://github.com/CenikLab/ribo-itp_paper/blob/main/ribo_and_rna_cds_counts.csv.gz).

## [References (in another repository)](https://github.com/CenikLab/ribo-itp_paper_references) 

[In a separate repository](https://github.com/CenikLab/ribo-itp_paper_references), we provide the reference files for the mouse and human transcriptomes used in this study. These files are used for running RiboFlow (for ribosome profiling experiments) and our RNA-Seq pipelines. 

In the mouse transcriptome, we masked the nucleotides, overlapping the SNPs, with Ns. 

Also, in human and mouse  rtRNA filters, we added some extra sequences for non-coding RNAs.
