library(Seurat)
library(SeuratDisk)
library(here)
library(magrittr)
library(xfun)
library(glue)
library(future)
library(furrr)
library(BiocParallel)
library(BiocParallel.FutureParam)
library(batchelor)
library(sjlabelled)
library(tidyverse)

plan(multicore, workers = 1)

hbp_genes <- c(
    "GFAP", "NAGK", "GNPNAT1", "PGM3",
    "UBAP1", "OGT", "MGEA5", "WNK1",
    "REL", "GFPT2"
    )

#' Filter Cells in a seurat object
#' 
#' @param seu a seurat object
#' @param sd_filtx standard devation to use during filtering
#' @param mt_pct mitocondrial percent to use during filtering
#' 
#' @importClassesFrom Seurat Seurat
#' @export

FilterReads <- function(seu, sd_filtx = 2.5, mt_pct = 10) {

  sc_metric <- data.table::as.data.table(
      seu@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt")]
      )

  nCount_sdv <- sd(seu$nCount_RNA)
  nCount_meanv <- mean(seu$nCount_RNA)

  nFeature_sdv <- sd(seu$nFeature_RNA)
  nFeature_meanv <- mean(seu$nFeature_RNA)
  nc_cutoffs <- c(nCount_meanv - sd_filtx * nCount_sdv, nCount_meanv + sd_filtx * nCount_sdv)
  nf_cutoffs <- c(nFeature_meanv - sd_filtx * nFeature_sdv, nFeature_meanv + sd_filtx * nFeature_sdv)

  if (mt_pct < 0) {
    mt_meanv <- mean(seu$percent.mt)
    mt_sdv <- sd(seu$percent.mt)
    mt_pct <- min(c(mt_meanv + mt_sdv, 25.))
  }

  nf_cutoffs <- c(nFeature_meanv - sd_filtx * nFeature_sdv, nFeature_meanv + sd_filtx * nFeature_sdv)

  seu <- seu[, seu@meta.data$nCount_RNA >= nc_cutoffs[[1]] &
                       seu@meta.data$nCount_RNA <= nc_cutoffs[[2]] &
                       seu@meta.data$nFeature_RNA >= nf_cutoffs[[1]] &
                       seu@meta.data$nFeature_RNA <= nf_cutoffs[[2]] &
                       seu@meta.data$percent.mt < mt_pct]
  return(seu)
}



# load in GSE117498 data

files <- list.files(here("GSE117498"), full.names = TRUE)

file_list <- map(files, ~ read_tsv(.x, col_types = cols()))
file_list <- map2(file_list, 1:length(file_list), ~ set_colnames(.x, paste(colnames(.x), .y, sep = "_")))
md <- map(file_list, colnames)
names(md) <- files %>%
    basename() %>%
    str_remove(".raw_counts.tsv.gz") %>%
    str_remove("GSM[:digit:]+_") 

meta <- md %>% enframe() %>% unnest() %>% as.data.frame() %>% column_to_rownames(var = "value")

names(file_list) <- names(md)

file_list <- map2(
    file_list,
    1:length(file_list),
    ~ rename(.x, Barcode = glue("Barcode_{.y}"))
)
gse <- reduce(file_list, full_join, by = c("Barcode")) %>%
    filter(Barcode != "Library") %>%
    column_to_rownames(var = "Barcode") %>%
    as.matrix()

gse[is.na(gse)] <- 0

gse_seurat <- CreateSeuratObject(counts = gse)
gse_seurat <- AddMetaData(gse_seurat, meta)
gse_seurat[["percent.mt"]] <- PercentageFeatureSet(gse_seurat, pattern = "^MT-")
gse_seurat <- FilterReads(gse_seurat)

seurat <- LoadH5Seurat(
    file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA")
    )

good_features <- intersect(rownames(seurat), rownames(gse_seurat))
seurat <- seurat[good_features, ]

gse_seurat <- gse_seurat[good_features, ]

# Batch Correction
seurat_ls <- SplitObject(seurat, split.by = "patient_id")

seurat_ls[["GSE"]] <- gse_seurat

seurat_ls <- future_map(seurat_ls, NormalizeData)

to_bc <- future_map(seurat_ls, as.SingleCellExperiment)

options(future.globals.maxSize = 40*1024^3) # 40 GB
# MNN was chosen because the goal is to do DE
batch_corrected <- cache_rds(
    batchelor::mnnCorrect(
        to_bc,
        BPPARAM = FutureParam(),
        correct.all = TRUE
    )
)

seurat <- as.Seurat(batch_corrected)

plan(sequential)

seurat <- cache_rds(
    SCTransform(seurat, vars.to.regress = "patient_id", conserve.memory = TRUE),
    file = "int_sct_transform.rds"
)

plan(multicore)

cache_rds({
    seurat <- RunPCA(seurat)
    seurat <- RunUMAP(seurat, dims = 1:30)
    seurat <- FindNeighbors(seurat, dims = 1:30)
    FindClusters(seurat)
    },
    file = "int_dim_red.rds"
)

UMAPPlot(seurat, group.by = "patient_id")
