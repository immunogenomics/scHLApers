#!/usr/bin/env Rscript

# Processes HLA alignment data in *.nuc files and return a data frame, 
# filling in missing bases for any allele with the bases from the closest allele.

suppressPackageStartupMessages({
    library(Biostrings)
    library(dplyr)
    library(purrr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(rtracklayer)
})

source('hlaseqlib_code.R') # contains the hla_compile_index function definition

if (! dir.exists('IMGTHLA/alignments_nuc_aligned_filled')) {
    dir.create('IMGTHLA/alignments_nuc_aligned_filled') # make output directory
}

for (g in c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1')) {
    hladb <- tibble(locus = g) %>%
        mutate(data = map(locus, ~hla_compile_index(., 'IMGTHLA', imgtfile = 'nuc'))) %>%
        filter(!is.na(data)) %>%
        unnest(data) %>%
        filter(!grepl("N$", allele)) %>%
        select(-locus) %>%
        mutate(allele = paste0("IMGT_", allele)) %>%
        split(.$allele) %>%
        map_chr("cds") %>%
        DNAStringSet() %>%
        writeXStringSet(paste0('IMGTHLA/alignments_nuc_aligned_filled/', g, '_nuc_aligned_filled.fa'))   
}