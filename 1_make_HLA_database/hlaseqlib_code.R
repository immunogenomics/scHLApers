# Code and functions originally written by Vitor Aguiar 
# from https://github.com/genevol-usp/hlaseqlib/tree/master/R
# Minorly adapted for this project

code_alleles <- function(genotypes) {
    allele_codes <- genotypes %>%
        dplyr::select(locus, allele) %>%
        dplyr::arrange(locus, allele) %>%
        dplyr::distinct() %>%
        dplyr::group_by(locus) %>%
        dplyr::mutate(code = seq_len(n())) 

    genotypes %>% dplyr::left_join(allele_codes, by = c("locus", "allele"))
}

alleles_to_groups <- function(df, input_fields = 4L) {

    hla_groups <- dplyr::filter(hla_groups, locus %in% df$locus) %>%
        dplyr::mutate(allele = hla_trimnames(allele, input_fields)) %>%
        dplyr::distinct()

    df %>% dplyr::group_by(subject, locus) %>%
        dplyr::mutate(e = 1:n()) %>%
        tidyr::separate_rows(allele, sep = "=") %>%
        dplyr::mutate(h = 1:n()) %>%
        tidyr::separate_rows(allele, sep = "/") %>%
        dplyr::mutate(allele = gsub("^IMGT_", "", allele) %>%
                  hla_trimnames(fields = input_fields)) %>%
        dplyr::left_join(hla_groups, by = c("locus", "allele")) %>%
        dplyr::mutate(group = ifelse(is.na(group), allele, group)) %>%
        dplyr::group_by(subject, locus, h) %>%
        dplyr::summarize(e = unique(e),
                 allele = paste(unique(group), collapse = "/")) %>%
        dplyr::group_by(subject, locus, e) %>%
        dplyr::summarize(allele = paste(unique(allele), collapse = "=")) %>%
        dplyr::ungroup() %>% 
        dplyr::select(-e) %>%
        dplyr::arrange(subject, locus, allele)
}

hla_attribute_seq <- function(incomplete_seq, complete_seq) {
    incomplete_seq %<>% strsplit("") %>% unlist
    complete_seq %<>% strsplit("") %>% unlist

    unseq <- incomplete_seq == "*"

    incomplete_seq[unseq] <- complete_seq[unseq]

    paste(incomplete_seq, collapse = "")
}

hla_format_sequence <- function(cds) {
    cds %>%
        gsub("\\.|\\|", "", .) %>%
        gsub("\\*", "N", .)
}

make_dist_matrix <- function(hla_df) {
    complete_alleles <- hla_df %>% 
        dplyr::filter(!grepl("\\*", cds)) %>%
        dplyr::pull(allele)

    incomplete_alleles <- hla_df %>% 
        dplyr::filter(grepl("\\*", cds)) %>%
        dplyr::pull(allele)

    cds_sequenced <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
        apply(2, function(x) !any(x == "*")) # find region(s) sequenced in all alleles

    run <- rle(cds_sequenced)

    ends <- cumsum(run$lengths)
    starts <- ends - run$lengths + 1L

    run_df <- tibble::tibble(value = run$values, start = starts, end = ends) %>%
        dplyr::filter(value == TRUE) %>%
        dplyr::slice(which.max(end - start + 1)) # find the longest run where all alleles are sequenced
    
    # Record the region of alignment used for each HLA gene:
    # A: 302-3564 (3262)
    # B: 293-3225 (2932)
    # C: 510-3591 (3081)
    # DRB1: 620-18160 (17540)
    # DQA1: 747-6324 (5577)
    # DQB1: 533-7549 (7016)
    # DPA1: 524-5465 (4941)
    # DPB1: 302-10716 (10414)
              
    hla_df_cds_common <- hla_df %>%
        dplyr::mutate(cds = substring(cds, run_df$start, run_df$end))

    hla_df_cds_common_complete <- hla_df_cds_common %>%
        dplyr::filter(allele %in% complete_alleles)

    hla_df_cds_common_incomplete <- hla_df_cds_common %>%
        dplyr::filter(allele %in% incomplete_alleles)

    cds_common_complete <- hla_df_cds_common_complete$cds
    names(cds_common_complete) <- hla_df_cds_common_complete$allele
    
    cds_common_incomplete <- hla_df_cds_common_incomplete$cds
    names(cds_common_incomplete) <- hla_df_cds_common_incomplete$allele

    stringdist::stringdistmatrix(cds_common_incomplete, cds_common_complete,
                                 method = "hamming", useNames = "names", nthread = 1)
}
    
make_closest_allele_df <- function(distmatrix) {

    distmatrix %>%
        split(seq_len(nrow(distmatrix))) %>%
        setNames(rownames(distmatrix)) %>%
        lapply(function(x) which(x == min(x, na.rm = TRUE))) %>% 
        purrr::map(~tibble::tibble(closest = colnames(distmatrix)[.])) %>%
        dplyr::bind_rows(.id = "inc_allele")
}
           
## Function from hlaseqlib

#' Process HLA alignment data in *.nuc or *.gen files and return a data.frame.
#'
#' @param locus character string. The HLA locus to be processed (e.g., "A", "DRB1").
#' @param imgtdb character string. Path to the IMGTHLA directory.
#' @param imgtfile character string. Whether to process "nuc" or "gen" files.  # should be gen for genomic
#' @param exons integer. Range of exons/introns (default is all exons/introns).
#' @param by_exon logical. Whether to return sequences separated for each exon/intron (default FALSE).
#' @param keep_sep logical. Whether to keep '|' as the separator between exons/introns.
#'
#' @return A data.frame.

hla_read_alignment = function(locus, imgtdb, imgtfile = c("nuc", "gen"), 
                              exons = NULL, by_exon = FALSE, keep_sep = FALSE) {
    imgtfile <- match.arg(imgtfile)
    if (imgtfile == "nuc") {
        locus_file <- ifelse(grepl("DRB\\d", locus), "DRB", locus)
        alig_file <- file.path(imgtdb, "alignments", paste0(locus_file, "_nuc.txt"))
    } else if (imgtfile == "gen") { ## Use this one for genomic sequences
        locus_file <- locus
        alig_file <- file.path(imgtdb, "alignments", paste0(locus_file, "_gen.txt"))
    }

    alignments <- readLines(alig_file) %>%
        gsub("\\s{2,}", " ", .) %>% ## 2 or more s's?
        trimws %>% ## trim white space trailing
        .[grepl(sprintf("^%s\\d?\\*\\d{2,3}[:A-Z0-9]*\\s", locus_file), .)] ##?

    hla_df <- tibble::tibble(allele = gsub("^(\\S+)\\s(.*)$", "\\1", alignments),
                             cds = gsub("\\s", "", gsub("^(\\S+)\\s(.*)$", "\\2", alignments))) %>%
        tibble::rowid_to_column() %>% 
        dplyr::group_by(allele) %>%
        dplyr::summarize(rowid = min(rowid), cds = paste(cds, collapse = "")) %>%
        dplyr::arrange(rowid) %>%
        dplyr::select(-rowid)

    hla_df$cds <- stringr::str_split(hla_df$cds, "", simplify = TRUE) %>%
        apply(2, function(i) {i[i == "-"] <- i[1]; i}) %>%
        apply(1, . %>% paste(collapse = ""))

    hla_df <- hla_df %>%
        dplyr::filter(sub("^([^\\*]+).+$", "\\1", allele) == locus) 

    if (by_exon || (!is.null(exons) && is.numeric(exons))) {
        hla_df <- hla_df %>%
            dplyr::mutate(cds = strsplit(cds, "\\|")) %>%
            tidyr::unnest(cds) %>%
            dplyr::group_by(allele) %>%
            dplyr::mutate(exon = seq_len(dplyr::n())) %>%
            dplyr::select(allele, exon, cds) %>% 
            dplyr::ungroup()
    
        if (imgtfile == "gen") {
            hla_df <- hla_df %>%
                dplyr::rename(idx = exon) %>%
                dplyr::mutate(feature = ifelse(idx %% 2 == 0, "exon", "intron"),
                          feature = dplyr::case_when(feature == "intron" & idx == 1 ~ "UTR",
                                                     feature == "intron" & idx == dplyr::last(idx) ~ "UTR",
                                                     TRUE ~ feature)) %>%
                dplyr::group_by(allele, feature) %>%
                dplyr::mutate(idx_grp = seq_len(dplyr::n())) %>% 
                dplyr::ungroup() %>%
                dplyr::select(allele, feature, idx, idx_grp, cds)
        } else if (imgtfile == "nuc") {
            hla_df <- hla_df %>%
            dplyr::mutate(feature = "exon", idx = exon, idx_grp = exon) %>%
            dplyr::select(allele, feature, idx, idx_grp, cds)
        }

        if (!is.null(exons) && is.numeric(exons)) {
            hla_df <- hla_df %>%
            dplyr::filter(feature == "exon", idx_grp %in% exons)
        }
    }
    
    if (!by_exon && !keep_sep) {
        hla_df <- hla_df %>% 
        dplyr::mutate(cds = gsub("\\|", "", cds))
    }
    return(hla_df)
}
               
#' Trim HLA allele names to the specified number of fields.
#'
#' Given a vector of HLA allele names with the format "A*01:01:01:01", returns
#' all or unique allele names at the digit resolution specified by 
#' \code{fields}.
#'
#' @param alleles_vec character string. Vector of allele names.
#' @param fields integer. Number of fields to trim names to.
#'
#' @return All or unique HLA alleles trimmed names.
#' @examples
#' hla_trimnames(c("A*01:01:01:01", "A*02:01:01:01"), fields = 2)
#' @export

hla_trimnames <- function(alleles_vec, fields = 3) {

  regex <- sprintf("^([A-Z]{1,3}\\d?\\*)((:?\\d{2,3}[NQLS]?){%d}).*$", fields)

  strsplit(alleles_vec, "/") %>%
      purrr::map(~sub(regex, "\\1\\2", .)) %>%
      purrr::map(unique) %>%
      purrr::map_chr(~paste(., collapse = "/"))
  return(alleles_vec)
}

               
## Modified this function to work for *.gen files by adding 'gen' argument.

#' Process HLA alignment data in *.nuc files and return a data.frame.
#'
#' @param locus character string.
#' @param imgt.database character string.
#'
#' @return A data.frame.

hla_compile_index <- function(locus, imgt.database, imgtfile = 'gen') {

    message(paste("Processing locus", locus)) 

    # Process alignments
    hla_df <- hla_read_alignment(locus, imgt.database, imgtfile=imgtfile) ## Added gen argument here

    if (all(grepl("\\*", hla_df$cds))) {
        message(paste0("No complete sequence for locus ", locus, "; returning NA"))
        return(NA)
    } else if (nrow(hla_df) == 1L || all(!grepl("\\*", hla_df$cds))) {
        # If N_alleles = 1 or all alleles are complete, output:
        final_df <- hla_df
    } else {
        # find closest complete allele for each incomplete allele
        distmatrix <- make_dist_matrix(hla_df) 

        closest_allele_df <- make_closest_allele_df(distmatrix) %>%
            dplyr::mutate(id = dplyr::group_indices(., inc_allele))
        closest_allele_df_step2 <- closest_allele_df %>% 
            tidyr::gather(ix, allele, 1:2) %>%
            dplyr::distinct(id, allele) %>%
            dplyr::left_join(hla_df, by = "allele") %>%
            split(.$id) %>%
            purrr::map(make_dist_matrix) %>%
            purrr::map(make_closest_allele_df) %>%
            dplyr::bind_rows()

        closest_within_type <- closest_allele_df_step2 %>%
            dplyr::mutate(`1` = hla_trimnames(inc_allele, 1) == hla_trimnames(closest, 1),
                          `2` = hla_trimnames(inc_allele, 2) == hla_trimnames(closest, 2),
                          `3` = hla_trimnames(inc_allele, 3) == hla_trimnames(closest, 3),
                          `4` = hla_trimnames(inc_allele, 4) == hla_trimnames(closest, 4)) %>%
            tidyr::gather(field, value, `1`:`4`) %>%
            dplyr::group_by(inc_allele) %>%
            dplyr::filter(all(value == FALSE) | value == TRUE) %>%
            dplyr::slice(which.max(field)) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(inc_allele) %>%
            dplyr::select(inc_allele, closest)

        inferred_df <- closest_within_type %>%
            dplyr::left_join(hla_df, by = c("inc_allele" = "allele")) %>%
            dplyr::left_join(hla_df, by = c("closest" = "allele")) %>%
            dplyr::mutate(cds = purrr::map2_chr(cds.x, cds.y, hla_attribute_seq)) %>%
            dplyr::select(allele = inc_allele, cds)

        final_df <- hla_df %>%
            dplyr::filter(!grepl("\\*", cds)) %>%
            dplyr::bind_rows(inferred_df) %>%
            dplyr::arrange(allele)
    }

    final_df %>% dplyr::mutate(cds = hla_format_sequence(cds)) %>%
        dplyr::mutate(allele3f = hla_trimnames(allele, 3)) %>%
        dplyr::distinct(allele3f, cds, .keep_all = TRUE) %>%
        dplyr::group_by(allele3f) %>%
        dplyr::mutate(n = dplyr::n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(allele = ifelse(n > 1L, allele, allele3f)) %>%
        dplyr::select(allele, cds) %>%
        dplyr::arrange(allele)
}