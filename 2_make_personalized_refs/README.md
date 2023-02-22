# Step 2: Make personalized reference and annotation files

## Download the reference genome and annotations
You will need to download the original reference genome and GTF files to [../example_data/GRCh38_refs](../example_data/GRCh38_refs) (files are too large for GitHub).

* Reference genome (e.g. GRCh38.primary_assembly.genome.fa): available [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)
* Gene annotation file (e.g. gencode.v38.annotation.gtf): available [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)

We recommend filtering the starting GTF file to genes types of interest (e.g. protein_coding, lncRNA, TCR, Ig) using the steps in [steps_to_filter_original_gtf.sh](steps_to_filter_original_gtf.sh).

## Tutorial
The [tutorial](2_make_personalized_refs/tutorial_make_pers_and_mask_GRCh38.ipynb) demonstrates how to generate personalized contigs (FASTA) and annotations (GTF) files (that will be combined with the masked reference) and how to mask the reference files.