# Snp-itp

List of scripts for the SNP calling on mouse ITP data.

Note that many of our scripts are dependent on the conventions of the VCF, Fasta and GTF files being used. For example, the fasta file header condtains transcript id, CDS positions etc.

## VCF Coordinate Conversion

Our original VCF file `CAST.SNPs.validated.vcf.gz`
([taken from SmartSeq3 repo](https://github.com/sandberg-lab/Smart-seq3/blob/master/allele_level_expression/CAST.SNPs.validated.vcf.gz)) is in **genomic** coordinates. The vcf file containing the transcriptomic coordinates is `transcriptomic_variants.vcf.gz`.

We transformed the coordinates in the original VCF file to transcriptomic coordinates. We used our GTF file, `gencode.vM25.annotation.gtf.gz`, for this conversion. 

The code used to convert the coordinates is in the notebook:
`translate_genomic_to_transcriptomic_coord.ipynb`.

We also masked the SNPs in the reference transcriptome  fasta file and obtained the masked fsata file 
`variant_masked_mouse_transcriptome.fa.gz`.
We used this file to map the reads against to call the SNPs.


