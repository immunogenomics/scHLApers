# Step 3: Quantify single-cell expression with STARsolo

See the [script](Run_STARsolo_scHLApers_BAMinput.sh), which runs the personalized alignment for one sample. 
Arguments:
1. Sample name
2. Path to STAR executable

## STARsolo parameters

Modify the script based on your dataset:
- `UMIlen` (12 for 10x v3, 10 for 10x v2)
- `whitelist` (barcode whitelist file)
- `soloType` (CB_UMI_Simple for droplet-based data)
- `cell_meta` (.csv file containing cell metadata, used by the `starsolo_to_genesXcells.R` to filter the resulting raw output expression matrix to remove empty droplets. If you do not have a list of desired cells, see `--soloCellFilter` option in the [STARsolo manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- `numThreads` (for parallelization)
- `readsDir` (directory with input bam files named "XXX.bam", where "XXX" is sample name)
- Change filepaths for `dir` (output directory), masked GRCh38 genome and annotation files, and personalized reference files, as needed.

Do not modify (required for scHLApers):
- `soloFeatures` (GeneFull_Ex50pAS required for intronic and exonic counting)
- `soloMultiMappers` (EM required for distributing multimapping reads)

## How to change the script to use FASTQ format for input

The example input data is in BAM format. If your input files are FASTQ format, input the FASTQ files in the order of:
`--readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz`

Also, replace the arguments:
`--readFilesType SAM SE --readFilesCommand samtools view -F 0x100` with `--readFilesCommand zcat`.

## Output files

The raw counts matrix output by the pipeline for example `Sample_1006_1007` can be found [here](../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/Sample_1006_1007_scHLApers_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx):
`../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/Sample_1006_1007_scHLApers_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx`

A filtered version is located [here](../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/exp_EM.rds) (read into R using `readRDS`):
`../example_outputs/STARsolo_results/Sample_1006_1007_scHLApers/exp_EM.rds`

Note: The classical HLA genes are named `IMGT_A`, `IMGT_C`, `IMGT_B`, `IMGT_DRB1`, `IMGT_DQA1`, `IMGT_DQB1`, `IMGT_DPA1`, `IMGT_DPB1`.