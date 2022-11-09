ReadNIHData <- function(path, keep_md = FALSE) {
    # Read in data from GDC Commons
    # Args:
    #   path: path to data (UUID subdirectories are allowed)
    # Returns:
    #   list of data.tables
    if (keep_md) {
        select <- c("gene_id", "gene_name", "gene_type", "unstranded")
    } else {
        select <- c("gene_name", "unstranded")
    }
    
    files <- list.files(path = path, pattern = ".*.tsv$", full.names = TRUE, recursive = TRUE)
    data <- future.apply::future_lapply(files, data.table::fread, select = select)
    x <- tools::file_path_sans_ext(basename(files))
    names(data) <- stringr::str_remove(x, ".rna_seq.augmented_star_gene_counts|.mirnaseq.isoforms.quantification.txt")
    return(data)
}

ReadStJudeData <- function(path) {
    # Read in data from St. Jude Cloud
    # Args:
    #   path: path to data (UUID subdirectories are not allowed)
    # Returns:
    #   list of data.tables
    files <- list.files(path = path, pattern = ".*.txt$", full.names = TRUE)
    data <- future_lapply(files, fread, header = FALSE, col.names = c("gene_name", "unstranded"))
    x <- tools::file_path_sans_ext(basename(files))
    names(data) <- str_remove(str_extract(x, "^.+_"), "_")
    return(data)
}

ReadStJudeAnnotation <- function(path) {
    # Read in annotation from St. Jude Cloud
    # Args:
    #   path: path to annotation file
    # Returns:
    #   data.table
    cols <- c("subject_name", "sample_name", "sj_long_disease_name", "attr_sex", "attr_age_at_diagnosis")
    ann <- fread(path, select = cols)
    new <- c("patient", "sample_id", "disease", "sex", "age")
    setnames(ann, cols, new)
    ann[, .id := sample_id]
}

ReadNIHAnnotation <- function(path) {
    # Read in annotation from GDC Commons
    # Args:
    #   path: path to annotation file
    # Returns:
    #   data.table
    cols <- c("case_id", "case_submitter_id", "primary_diagnosis", "gender", "age_at_diagnosis")
    ann <- fread(path, select = cols, header = TRUE)
    new <- c("patient", "sample_id", "disease", "sex", "age")
    setnames(ann, cols, new)
    ann[, .id := patient][, age := as.numeric(age) / 365]
}