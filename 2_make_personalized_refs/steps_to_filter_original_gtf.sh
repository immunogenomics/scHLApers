#!/bin/bash

gtf_in="gencode.v38.annotation.gtf"

BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""

cat "$gtf_in" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "gene_allowlist"

# Filter the GTF file based on the gene allowlist
gtf_filtered="gencode.v38.annotation.filtered.gtf"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_in" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "gene_allowlist" "$gtf_in" \
    >> "$gtf_filtered"


