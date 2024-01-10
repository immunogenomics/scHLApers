# scHLApers
Code to run the scHLApers pipeline for quantifying **s**ingle-**c**ell **HLA** expression using **pers**onalized reference genomes (Kang et al., Nat Genetics 2023).
![Overview](images/overview.png)

## Requirements
R program requires (listed version or higher):
* R=4.0.5
* Biostrings=2.58.0
* purrr=0.3.4
* readr=2.1.2
* stringi=1.7.8
* stringr=1.4.0
* tidyverse=1.3.1
* rtracklayer=1.50.0

Other software:
* STAR=2.7.10a [https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)
* samtools=1.4.1 [http://www.htslib.org/download/](http://www.htslib.org/download/)

Data:
* Reference genome (e.g. GRCh38.primary_assembly.genome.fa): available [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)
* Gene annotation file (e.g. gencode.v38.annotation.gtf): available [here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)
* Cell barcode whitelist: more info [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)

## Pipeline and example data
Each step has its own directory with necessary scripts and a tutorial walking through the steps. The [example_data](example_data) and [example_output](example_outputs) directories contain example input and output files for 2 samples. The raw scRNA-seq data for the example was obtained from [Yazar et al. Science 2022](https://pubmed.ncbi.nlm.nih.gov/35389779/) study, publicly available on GEO (GSE196830).

### Input
The inputs to scHLApers are:
* Raw scRNA-seq data (either FASTQ or BAM format)
* HLA allele calls (in CSV format, labeled as "SampleX_alleles.csv", see `example_data/inputs/alleles` for format)

See the [HLA analyses tutorial](https://github.com/immunogenomics/HLA_analyses_tutorial) from [Sakaue et al.](https://www.biorxiv.org/content/10.1101/2022.08.24.504550v1) for protocol for imputing HLA alleles from genotype array data.

### Step 1: Prepare HLA allelic sequence database
We provide a [pre-prepared database](1_make_HLA_database/IMGTHLA_all_alleles_FINAL.fa) generated from IPD-IMGT/HLA version 3.47 that can be directly used in Step 2. Alternatively, you can prepare your own database using the latest IPD-IMGT/HLA verison following the [tutorial](1_make_HLA_database/tutorial_make_database.ipynb).

### Step 2: Make personalized reference and annotation files
The [tutorial](2_make_personalized_refs/tutorial_make_pers_and_mask_GRCh38.ipynb) demonstrates how to generate personalized contigs (FASTA) and annotations (GTF) files (that will be combined with the masked reference) and how to mask the reference.

### Step 3: Quantify single-cell expression with STARsolo
Example scripts for how to run [STARsolo](https://github.com/alexdobin/STAR) for read alignment and expression quantification in single-cell data. Script will need to be modified based on the specifics of your dataset (e.g. UMI length, input format, barcode whitelist path, STAR executable). Please see the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for all options.

### Outputs
The output of scHLApers is a genes by cells expression matrix, with improved classical HLA expression estimates. In the example output, we have filtered the raw STARsolo counts matrix (to remove empty droplets) using a provided list of cell barcodes (see `example_data/cell_meta_example.csv`).

The raw counts matrix output by the pipeline for example `Sample_1006_1007` can be found [here](../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/Sample_1006_1007_scHLApers_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx):
`../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/Sample_1006_1007_scHLApers_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx`

A filtered version is located [here](../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/exp_EM.rds) (read into R using `readRDS`):
`../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/exp_EM.rds`

Note: The classical HLA genes are named `IMGT_A`, `IMGT_C`, `IMGT_B`, `IMGT_DRB1`, `IMGT_DQA1`, `IMGT_DQB1`, `IMGT_DPA1`, `IMGT_DPB1`.

## Support
For questions and assistance not answered in tutorials, you can contact Joyce Kang (joyce_kang AT hms.harvard DOT edu).

## Reproducing results from the manuscript
Code to reproduce the figures and analyses from Kang et al. will become available at [https://github.com/immunogenomics/hla2023](https://github.com/immunogenomics/hla2023).
