#!/usr/bin/env Rscript

# Converts STARsolo output for 1 sample into a genesXcells matrix, subset the full matrix
# output to only include high-quality cells.

# Input file format requirements:
#  - sample: sample name
#  - STARsolo_sample_dir: directory of STARsolo output for the sample
#  - label: label to add to file names (must match the "label" used in Run_STARsolo_scHLApers)
#  - cell_metadata_path: file must be a .csv file, including columns named "Cell" and "Sample" for 
#    cell barcodes and sample names (used to subset the full matrix)

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(stats)
    library(tidyverse)
    library(Matrix)
    library(tidyr)
})

args = commandArgs(trailingOnly=TRUE)
sample = args[1] # sample name
STARsolo_results_dir = args[2] 
label=args[3] # label to add to file names
cell_meta_path = args[4]  # path to cell metadata (.csv) file
STARsolo_sample_dir = paste0(STARsolo_results_dir, sample, '_', label, '/')
cell_meta = read.csv(cell_meta_path)

# get cell and gene names
result_cells = readLines(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/barcodes.tsv'))
result_genes = readLines(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/features.tsv'))

get_gene_name = function(gene) { return(str_split(gene, '\t')[[1]][2]) }
result_gene_names = unlist(lapply(result_genes, get_gene_name))

# build starting matrices
result_EM_mat = as(readMM(paste0(STARsolo_sample_dir, sample, '_', label, '_Solo.out/GeneFull_Ex50pAS/raw/UniqueAndMult-EM.mtx')), 'dgCMatrix') # EM 
dimnames(result_EM_mat) = list(result_gene_names, result_cells)
colnames(result_EM_mat) = paste0(sample, '_', colnames(result_EM_mat))

cells_sample = cell_meta$Cell[which(cell_meta$Sample == sample)] # Get list of cells for that sample

# subset matrix by QC'd cells
result_EM_mat = result_EM_mat[, which(colnames(result_EM_mat) %in% cells_sample)]

# save
saveRDS(result_EM_mat, paste0(STARsolo_sample_dir, 'exp_EM.rds'))