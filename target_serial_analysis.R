library(data.table)
library(janitor)
library(edgeR)
library(limma)
library(stringr)
library(pheatmap)
source("lib/import.R")

# use limma::eBayes
# try correcting for clincial factors (to account for heterogeneity)
# make heatmap with annotations for OGA, OGT
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


# make DGEList

range(colSums(x))[[2]] / range(colSums(x))[[1]]
dge <- DGEList(
    x,
    samples = col_data[id %in% colnames(x),],
    genes = row_data
)

dge <- dge[
    dge$gene$gene_type == "protein_coding",
    dge$samples$timepoint %in% c("Primary Blood Derived Cancer", "Recurrent Blood Derived Cancer")
]
dge <- calcNormFactors(dge)
hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UAP1", "OGT", "OGA", "GFPT1", "GFPT2")
dge$samples$timepoint <- droplevels(relevel(dge$samples$timepoint, ref = "Primary Blood Derived Cancer"))
design <- model.matrix(~ timepoint, data = dge$samples)
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v[v$genes$gene_name %in% hbp_genes, ], design)
fit <- eBayes(fit)

top <- topTable(fit[fit$genes$gene_name %in% c("OGA", "OGT"), ])

exprs <- v[v$genes$gene_name %in% c("OGA", "OGT"), ]
rownames(exprs) <- exprs$genes$gene_name
exprs <- as.data.table(t(exprs$E), keep.rownames=TRUE)

setkey(col_data, "id")
setkey(exprs, "rn")
n_sam_per_pat <- exprs[col_data][ , .(n = .N), by = "case_id"]
setorder(n_sam_per_pat, -n)

merged <- exprs[col_data]
merged <- merged[,
    .(OGT = mean(OGT), OGA = mean(OGA)),
    by = sample_id
][col_data, on = "sample_id"]
merged <- merged[timepoint %in% c(
        "Primary Blood Derived Cancer",
        "Recurrent Blood Derived Cancer"
    ), timepoint := relevel(timepoint, "Primary Blood Derived Cancer")]
mean <- function(...) base::mean(..., na.rm = TRUE)
oga <- dcast(merged, case_id ~ timepoint, value.var = "OGA", fun.aggregate = mean)
oga <- clean_names(oga)
oga[ , fc := log2(primary_blood_derived_cancer / recurrent_blood_derived_cancer)]
oga <- oga[!is.na(fc),]

library(viridisLite)
library(readxl)

dis <- read_excel("data/TARGET/TARGET_AML_ClinicalData_Discovery_20221108.xlsx")
val <- read_excel("data/TARGET/TARGET_AML_ClinicalData_Validation_20221108.xlsx")

sam_info <- rbind(dis, val)
sam_info <- clean_names(sam_info)
setDT(sam_info)

mat <- oga[sam_info, on = .(case_id = target_usi)]
mat <- mat[!is.na(fc), ]
features <- c("flt3_itd_positive", "flt3_itd_allelic_ratio",
              "npm_mutation", "cebpa_mutation", "wt1_mutation",
              "c_kit_mutation_exon_8", "c_kit_mutation_exon_17",
              "bone_marrow_leukemic_blast_percentage_percent", "fab_category",
              "t_6_9", "t_8_21", "t_3_5_q25_q34","t_6_11_q27_q23",
              "t_9_11_p22_q23", "t_10_11_p11_2_q23",
              "t_11_19_q23_p13_1", "inv_16" , "del5q", "del7q", "del9q",
              "monosomy_5", "monosomy_7", "trisomy_8", "trisomy_21",
              "mll", "minus_y", "minus_x"
)
ann_row <- mat[, ..features]

hm_data <- as.matrix(mat[, fc])
rownames(hm_data) <- rownames(ann_row) <- 1:nrow(mat)
pheatmap(
    hm_data,
    color=magma(100),
    cellwidth = 10,
    show_rownames = FALSE,
    cluster_col=FALSE,
    annotation_row=ann_row,
    filename="hm_ann_OGA.pdf"
)

fwrite(oga, "bulk_oga_exp.tsv")