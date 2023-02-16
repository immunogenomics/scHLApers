## Useful R functions for creating the personalized contigs and GTF file
## By Amber Shen


# Formats the name of an HLA allele from HLA_X*ii:jj:kk to X-ii-jj
format_allele = function(allele) {
    gene = str_split(str_split(allele, '_')[[1]][2], '\\*')[[1]][1]
    two_digit = str_split(str_split(allele, '\\*')[[1]][2], ':')[[1]][1]
    four_digit = str_split(str_split(allele, '\\*')[[1]][2], ':')[[1]][2]
    return(paste(gene, two_digit, four_digit, sep = '-'))
}

# Makes the personalized .fa file for extra contigs to be added to genome
# Inputs:
# - personalized_alleles: .csv file containing personalized HLA alleles imputed for a sample
#    Assumes the file is labeled "sample_alleles.csv", where "sample" is a unique sample name.
# - genome_out: Directory for output
# Outputs:
# - Writes extra contigs to genome_out/sample_genome.fa
# - Appends any alleles missing from the database to genome_out/missing_alleles.csv
make_genome = function(personalized_alleles, genome_out) {
    idxs = c()
    warnings = c()
    sample = str_split(tail(str_split(personalized_alleles, '/')[[1]], n = 1), '_alleles')[[1]][1]
    alleles = read.csv(personalized_alleles)
    for (i in 1:nrow(alleles)) {
        
        to_match = format_allele(alleles[i, 'ID']) # HLA_X*ii:jj:kk to X-ii-jj
        idx = which(formatted_database_names == to_match) # get indices of matching sequences
        
        # If no matches, skip allele
        if (length(idx) == 0) {
            warnings = c(warnings, paste0('WARNING: ', to_match, ' not found for sample ', sample))
            next
        }
        idxs = c(idxs, idx[1])   
    }
    Biostrings::writeXStringSet(database[idxs], paste0(genome_out, sample, '_genome.fa')) # save file
    lapply(warnings, cat, '\n', file=paste0(genome_out, 'missing_alleles.csv'), append=TRUE) # saves missing alleles
}

# Makes the personalized annotation .gtf file to be added to Gencode gtf file
# Inputs:
# - sample: sample name
# - annot_out: Directory for output
# - unique: whether to label alleles uniquely in the GTF file (TRUE) or count them towards the overall gene (FALSE)
#          (default is FALSE)
# Outputs:
# - Writes personalized .gtf to annot_out file
make_annotation = function(sample, annot_out, unique = FALSE) {
    genome = Biostrings::readDNAStringSet(paste0(out, 'genomes/', sample, '_genome.fa'))
    col_names = c('seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')
    annot = data.frame(matrix(ncol=9, nrow=0, dimnames=list(NULL, col_names)))
    
    for (allele in names(genome)) {
        end = width(genome[allele])
        if (unique) {
            name = allele
        } else {
            name = str_split(allele, '\\*')[[1]][1] # nonunique case
        }
        attribute = paste0('transcript_id "', allele, '"; gene_id "', name, '"; gene_name "', name, '";')
        annot[nrow(annot)+1,] = c(allele, 'IMGTHLA', 'exon', 1, end, '.', '+', '.', attribute) 
    }
    write.table(annot, file=annot_out, sep='\t', quote = FALSE, col.names=FALSE, row.names=FALSE)
}

# Inputs:
# - personalized_alleles_path: directory containing Sample_alleles.csv files
#     Important: Make sure no other files are in the directory
# - sequence_database_path: database of HLA allele sequences (generated in 1_make_HLA_database
# - out: Output prefix
# Outputs:
# - Personalized FASTA (in out/genomes directory) and GTF files (in out/annotations directory)
make_references = function(personalized_alleles_path, sequence_database_path, out) {
    dir.create(out, showWarnings = FALSE)
    dir.create(paste0(out, 'genomes/'), showWarnings = FALSE)
    dir.create(paste0(out, 'annotations/'), showWarnings = FALSE)
    
    # Get list of personalized allele files
    personalized_alleles = list.files(personalized_alleles_path, full.names=TRUE)
    
    # make genomes
    for (i in 1:length(personalized_alleles)) {
        genome_out = paste0(out, 'genomes/')
        make_genome(personalized_alleles[i], genome_out)
    }
    
    # make annotations
    get_sample = function(file) {return(str_split(tail(str_split(file, '/')[[1]], n=1), '_alleles')[[1]][1])}
    samples = lapply(personalized_alleles, get_sample)

    for (sample in samples) {
        annot_out = paste0(out, 'annotations/', sample, '_annotation.gtf')
        make_annotation(sample, annot_out, unique=FALSE)
    }
}