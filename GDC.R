#! /usr/bin/env Rscript

library(data.table)
library(stringr)
library(future)
library(future.apply)
library(ggplot2)
library(here)
library(broom)
library(ggsignif)

plan(multicore)

ReadNIHData <- function(path) {
    # Read in data from GDC Commons
    # Args:
    #   path: path to data (UUID subdirectories are allowed)
    # Returns:
    #   list of data.tables
    files <- list.files(path = path, pattern = ".*.tsv$", full.names = TRUE, recursive = TRUE)
    data <- future_lapply(files, fread, select = c("gene_name", "unstranded"))
    x <- tools::file_path_sans_ext(basename(files))
    names(data) <- str_remove(x, ".rna_seq.augmented_star_gene_counts")
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
    ann[, .id := patient][ , age := as.numeric(age) / 365]
}

SimplifyDiseaseNames <- function(s) {
    # Simplify disease names
    # Args:
    #   s: string
    # Returns:
    #   string
    disease <- str_to_lower(s)
    fcase(
        disease %flike% "acute lymphoblastic leukemia", "ALL",
        disease %flike% "acute myeloid leukemia", "AML",
        disease %flike% "lymphoblastic leukemia", "ALL",
        disease %flike% "chronic myeloid leukemia", "CML",
        disease %flike% "mixed phenotype acute leukemia", "Mixed",
        disease %flike% "myloproliferative neoplasm", "MDS",
        disease %flike% "chronic myelomonocytic leukemia", "CML",
        disease %flike% "acute monoblastic and monocytic leukemia", "AML",
        disease %flike% "acute myelomonocytic leukemia", "AML",
        disease %flike% "acute promyelocytic leukemia", "AML",
        disease %flike% "acute megakaryoblastic leukaemia", "AML",
        disease %flike% "myeloid leukemia", "AML",
        disease %flike% "myelodysplastic", "MDS",
        disease %flike% "normal", "Normal",
        default = "Other"
    )
}

NormalizeByCPM <- function(x) (x / sum(x)) * 1000000
NormalizeByLogCPM <- function(x) log(((x + 1) / sum(x)) * 1000000)

sj_path <- "/fs/ess/PCCF0022/Datasets/Public_Leukemia/sj_counts"
nih_path <- "/fs/ess/PCCF0022/Datasets/Public_Leukemia/nih_counts"
nih_ann <- ReadNIHAnnotation("/fs/ess/PCCF0022/Datasets/Public_Leukemia/nih_patient_information.tsv")
sj_ann <- ReadStJudeAnnotation("/fs/ess/PCCF0022/Datasets/Public_Leukemia/stjude_patient_information.txt")
ann <- rbind(nih_ann, sj_ann)

nih_data <- ReadNIHData(nih_path)
sj_data <- ReadStJudeData(sj_path)

gtex_data <- fread("gene_reads_2017-06-05_v8_whole_blood.gct")
gtex_data[, id := NULL][, Name := NULL]
setnames(gtex_data, "Description", "gene_name")
gtex_data <- melt(gtex_data, id.vars = "gene_name", variable.name = ".id", value.name = "counts")
gtex_data[ , disease := "normal"][gene_name == "MGEA5", gene_name := "OGA"]

data <- c(nih_data, sj_data)

# get list of genes across all datasets
genes <- lapply(data, function(x) x$gene_name)
genes <- Reduce(intersect, genes)

message("Row binding list")
dt <- rbindlist(data, use.names = TRUE, idcol = TRUE)

setnames(dt, "unstranded", "counts")
setkey(dt, .id)
setkey(ann, .id)

# remove genes not in all datasets
message("Removing genes not in all datasets")
dt <- dt[gene_name %in% genes, .(counts = sum(counts)), by = .(gene_name, .id)]

# nest data by sample to allow for annotation
message("Nesting data by sample")
dt2 <- dt[, .(nested = list(.SD)), by = .id]
setkey(dt2, .id)
dt2 <- dt2[ann, roll = TRUE]
message("Unnesting data")
dt <- dt2[, rbindlist(nested), by = .(.id, patient, sample_id, disease, sex, age)]
dt <- rbindlist(list(dt, gtex_data), fill = TRUE)

message("Normalizing Data")
dt[ , cpm := NormalizeByCPM(counts), by = .id]
message("Log Normalizing")
dt[ , log_cpm := NormalizeByLogCPM(counts), by = .id]

message("Simplifying Disease names")
dt[, simple_disease := SimplifyDiseaseNames(disease)]

message("Stratifying Age Groups")
dt[, age_group := fifelse(age <= 29, "Pediatric", "Adult", na = "Unknown")]
dt[, age_disease := paste(age_group, simple_disease, sep = " ")]

message("Plotting")
hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UAP1", "OGT", "OGA", "GFPT2", "GFPT1")
data_to_plot <- dt[simple_disease %in% c("AML", "Normal") & gene_name %in% hbp_genes & (age_group != "Unknown" | simple_disease == "Normal"), ]
data_to_plot[, group := fifelse(simple_disease == "Normal", "Normal", age_group)][, group := factor(group, levels = c("Normal", "Pediatric", "Adult"))]

WrapPlot <- function(ggplot) {
    ggplot +
    geom_violin(mapping = aes(fill = group)) +
    geom_boxplot(fill = "white", width = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(
        values = c(
            "Pediatric" = "#B0BEA9",
            "Adult" = "#92AA83",
            "Normal" = "#F1F7EE"
        )
    ) +
    facet_wrap(~gene_name, strip.position = "left", scales = "free") +
    geom_signif(
        comparisons = list(
            c("Pediatric", "Normal"), 
            c("Adult", "Normal"), 
            c("Pediatric", "Adult")
        ),
        map_signif_level = TRUE,
        step_increase = 0.2,
        vjust = 0.2
    ) +
    theme(
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside"
    ) +
    labs(x = NULL)
}

pdf(here("plots/tumor_v_normal.pdf"), width = 10, height = 10)
ggplot(data = data_to_plot, mapping = aes(x = group, y = cpm)) %>%
    WrapPlot() +
    labs(y = "CPM")

ggplot(data = data_to_plot, mapping = aes(x = group, y = log_cpm)) %>%
    WrapPlot() +
    labs(y = "log(CPM)")
graphics.off()

RunANOVA <- function(data) {
    aov <- aov(log_cpm ~ group, data = data)
    aov_p <- tidy(aov)[["p.value"]][[1]]
    if (aov_p > 0.05) {
        dt <- data.table(adj_p_value = 1, contrast = "Failed ANOVA")
        return(dt)
    } 
    tuk <- TukeyHSD(aov)
    tuk_p <- tidy(tuk)[["adj.p.value"]]
    res <- ifelse(tuk_p == 0, 2e-16, tuk_p)
    res <- as.data.table(res)
    setnames(res, "res", "adj_p_value")
    res[, contrast := tidy(tuk)[["contrast"]]]
    return(res)
}

fwrite(data_to_plot[, RunANOVA(.SD), by = gene_name], file = here("results/tumor_v_normal.tsv"))
fwrite(data_to_plot[ , .N, by = .(gene_name, group)], file = here("results/tumor_v_normal_data.tsv"))
