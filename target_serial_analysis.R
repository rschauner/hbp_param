library(data.table)
library(janitor)
library(DESeq2)
library(stringr)
source("lib/import.R")

library(BiocParallel)
register(MulticoreParam(8))

# set up count data
data <- ReadNIHData("data/TARGET/gex", keep_md = TRUE)

dt <- rbindlist(data, use.names = TRUE, idcol = TRUE)
bad_features <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
dt <- dt[!(gene_id %in% bad_features), ]
setnafill(dt, fill=0, cols="unstranded")
x <- dcast(
     dt,
     gene_id ~ .id,
     fill = 0
)

rownames(x) <- x$gene_id
x$gene_id <- NULL

# set up row data
row_data <- dt[
    !duplicated(dt, by = c("gene_name", "gene_type", "gene_id")),
    c("gene_name", "gene_type", "gene_id")
][
    !(gene_type %in% c("pr", "unproces", "")),
]
rownames(row_data) <- row_data$gene_id

# set up col data
col_data <- fread("data/TARGET/gex_information.tsv")
col_data <- clean_names(col_data)
col_data[ , id := str_remove(file_name, ".rna_seq.augmented_star_gene_counts.tsv")]
col_data[ , timepoint := str_split(sample_type, " - ", simplify = TRUE)[, 1]]
col_data[ , type := str_split(sample_type, " - ", simplify = TRUE)[, 2]]
cols <- c("timepoint", "type")
col_data[ , paste0(cols) := lapply(.SD, as.factor), .SDcols = cols]
rownames(col_data) <- col_data$id


# make DESeqDataSet
design <- model.matrix(~ timepoint, data = col_data)
dds <- DESeqDataSetFromMatrix(
    x,
    colData = col_data[id %in% colnames(x),],
    rowData = row_data,
    design = design
)

dds <- DESeq(dds, parallel = TRUE)
rownames(dds) <- rowData(dds)$gene_name
design(dds) <- ~ timepoint
dvr <- results(
    dds,
    contrast = list(
        c("timepointPrimary.Blood.Derived.Cancer"),
        c("timepointRecurrent.Blood.Derived.Cancer")
    ),
    parallel = TRUE
)
dvr <- as.data.frame(dvr)
dvr$gene <- rownames(dvr)
fwrite(dvr, "target_prim_v_recur_de.tsv")

nvd <- results(
    dds,
    contrast = list(
        c("timepointBlood.Derived.Normal"),
        c("timepointPrimary.Blood.Derived.Cancer")
    ),
    parallel = TRUE
)
nvd <- as.data.frame(nvd)
nvd$gene <- rownames(nvd)
fwrite(nvd, "target_prim_v_norm_de.tsv")
