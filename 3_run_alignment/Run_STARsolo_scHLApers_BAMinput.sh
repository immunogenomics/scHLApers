#!/bin/bash
## Amber Shen and Joyce Kang
## Last updated Feb 21, 2023

## Script that runs the personalized (scHLApers) pipeline for a single sample,
## Input is the sample name. Script outputs all intermediary files
## and will automatically rerun any job that fails to generate the appropriate output files.

# Load samtools
module load samtools # can replace with your installation

# Set up input files and parameters
sample=$1 # sample name
STAR_executable=$2 # e.g. /PHShome/jbk37/tools/STAR-2.7.10a/source/STAR
numThreads=4
UMIlen=10 # 12 for 10x v3, 10 for 10x v2
soloType="CB_UMI_Simple" # type of CB/UMI used
whitelist="../example_data/barcode_whitelists/737K-august-2016.txt" # 10x v2 whitelist
# For 10x data, whitelists are packaged with Cell Ranger:
# path/to/cellranger/1.3.1/cellranger-cs/1.3.1/lib/python/cellranger/barcodes/737K-august-2016.txt

readsDir="../example_data/inputs/"
samtools collate ${readsDir}${sample}.bam ${readsDir}${sample}_shuffled # shuffle the input bam file
input_bam=${readsDir}${sample}_shuffled.bam

# Path to script to make genesXcells using a list of cell barcodes passing QC
make_exp_script="starsolo_to_genesXcells.R"
cell_meta="../example_data/cell_meta_example.csv"

##########################
# Run scHLApers Pipeline
##########################
echo 'Run scHLApers'

label="scHLApers"
dir="../example_outputs/"
mkdir -p ${dir}STARsolo_results # create output directory if does not exist
out="${dir}STARsolo_results/${sample}_${label}/${sample}_${label}_" # directory to store results
genomeDir="${dir}${sample}_${label}_index/" # store index

finished=false
while ! $finished; do

    ##### Run Personalized Pipeline #####
    rm -rf ${dir}STARsolo_results/${sample}_${label} # remove the output directory if it exists
    mkdir -p ${dir}STARsolo_results/${sample}_${label}
    
    rm -rf $genomeDir # remove the genome directory if it exists
    mkdir $genomeDir

    # Combine personalized genome with masked reference (GRCh38)
    GRCh38_genome="../example_data/GRCh38_refs/GRCh38.primary_assembly.genome_maskedHLAClassical.fa"
    GRCh38_annot="../example_data/GRCh38_refs/gencode.v38.annotation.filtered.masked.gtf"
    perGenome="../example_outputs/personalized_references/genomes/${sample}_genome.fa"
    perAnnot="../example_outputs/personalized_references/annotations/${sample}_annotation.gtf"
    combPerGenome="${dir}STARsolo_results/${sample}_${label}_genome.fa"
    combPerAnnot="${dir}STARsolo_results/${sample}_${label}_annot.gtf"
    cat $GRCh38_genome $perGenome > $combPerGenome
    cat $GRCh38_annot $perAnnot > $combPerAnnot

    # Generate Genome Index
    CMD="$STAR_executable --runThreadN $numThreads \
        --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $combPerGenome \
        --sjdbGTFfile $combPerAnnot --outTmpDir ${out}_STARtmp_genome_generate"
    $CMD

    # Align Reads
    CMD="$STAR_executable --genomeDir $genomeDir --readFilesIn $input_bam \
          --runThreadN $numThreads --runDirPerm All_RWX --soloUMIlen $UMIlen \
          --outFileNamePrefix $out --soloType $soloType --soloCBwhitelist $whitelist \
          --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR \
          --soloFeatures GeneFull_Ex50pAS \
          --soloMultiMappers EM --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
          --readFilesType SAM SE --readFilesCommand samtools view -F 0x100 \
          --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY \
          --readFilesSAMattrKeep None \
          --outSAMattributes NH HI nM AS CR UR GX GN CB UB \
          --outTmpDir ${out}_STARtmp"
    echo $CMD
    $CMD

    ## Check whether the run generated output
    if [ -f "${out}Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx" ] && [ -f "${out}Aligned.sortedByCoord.out.bam" ]; then
        echo "File exists, setting finished=true"
        finished=true
    else
        echo "Error: File does not exist, rerunning STARsolo personalized pipeline"
    fi
done

# Optional: Sort the output bam file and save HLA region and unmapped reads
bam="${out}Aligned.sortedByCoord.out.bam"
CMD="samtools index -@ 6 ${bam}"
echo $CMD
$CMD

pers_chrs=$(samtools idxstats ${bam} | cut -f 1 | grep "IMGT_") # extract names from bam file
pers_chrs=$(echo $pers_chrs | sed -e 's/ /: /g') # add ":" after every space
pers_chrs="${pers_chrs}:" # add ":" after final chromosome name

bam_HLA="${out}HLA.bam"
bam_unmapped="${out}unmapped.bam"
   
CMD="samtools view -b ${bam} chr6:28000000-34000000 ${pers_chrs} -o ${bam_HLA}"
$CMD
CMD="samtools view -b -f 4 ${bam} -o ${bam_unmapped}"
$CMD
CMD="samtools index -@ 6 ${bam_HLA}"
$CMD

# Optional: Remove original output bam file to save space
rm $bam
rm "${bam}.bai"

# Optional: Delete intermediate files to save space
rm $combPerGenome
rm $combPerAnnot
rm -R $genomeDir

# Subset to only the cells of interest
CMD="Rscript ${make_exp_script} $sample ${dir}STARsolo_results/ $label $cell_meta"
echo $CMD
$CMD

echo "Finished!"
