#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(harmony)
    library(here)
    library(readxl)
    library(SeuratDisk)
    library(patchwork)
    library(msigdbr)
    library(limma)
    library(EnhancedVolcano)
    library(future)
    library(tidyverse)
    library(writexl)
    library(xfun)
})

#plan("multicore")

seurat <- LoadH5Seurat(
    file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat.old"),
    assays = c("RNA", "SCT")
)

DefaultAssay(seurat) <- "SCT"

if (DefaultAssay(seurat) == "RNA") {
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
}

hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UAP1", "OGT", "OGA", "GFPT1", "GFPT2")

seurat <- cache_rds({
    seurat <- GetResidual(seurat, features = hbp_genes)
    seurat <- RunPCA(seurat, verbose = FALSE)
    seurat <- RunUMAP(seurat, dims = 1:20)
    seurat <- FindNeighbors(seurat, dims = 1:20)
    seurat <- FindClusters(seurat, resolution = 0.6)
    },
    file = "total_dim_red.rds"
)

subset_seurat <- list(
    paired_stem = subset(seurat, subset = stemness == "Stem" & paired == "T"),
    paired_nonstem = subset(seurat, subset = stemness == "Nonstem" & paired == "T"),
    diagnosis_stem = subset(seurat, subset = stemness == "Stem" & timepoint == "Diagnosis"),
    diagnosis_nonstem = subset(seurat, subset = stemness == "Nonstem" & timepoint == "Diagnosis")
)
groups <- list(
    paired_stem = "timepoint",
    paired_nonstem = "timepoint",
    diagnosis_stem = "prognosis",
    diagnosis_nonstem = "prognosis"
)

subset_seurat <- cache_rds(
    map(
        subset_seurat,
        ~ {
            .x <- RunPCA(.x, verbose = FALSE)
            .x <- RunUMAP(.x, dims = 1:20)
            .x <- FindNeighbors(.x, dims = 1:20)
            .x <- FindClusters(.x)
            .x
        }
    ),
    file = "subset_dim_red.rds"
)

#source(here("lib", "violin.R"))
#source(here("lib", "feature_plot.R"))
#source(here("lib", "dot_plot.R"))
#source(here("lib", "blend.R"))
source(here("lib", "DE.R"))