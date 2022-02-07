suppressPackageStartupMessages({
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
    library(future.batchtools)
    library(doFuture)
    library(batchelor)
    library(sjlabelled)
    library(tidyverse)
})

register(DoparParam()) ## Tell BiocParallel to use a foreach backend
registerDoFuture()     ## Tell foreach to use a future backend

resources_big   <- list(memory = 80, ncpus = 1, walltime = 167 * 60^2)
resources_small <- list(memory = 60, ncpus = 2, walltime = 72 * 60^2)
options(future.globals.maxSize = 40*1024^3) # 40 GB

SetPlan <- function(resources) {
    plan(
        list(
            tweak("batchtools_slurm", resources = resources),
            "multicore"
        )
    )
}

SetPlan <- function(resources) {
    plan(multicore)
} 

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
  nc_cutoffs <- c(
      nCount_meanv - sd_filtx * nCount_sdv,
      nCount_meanv + sd_filtx * nCount_sdv
    )
  nf_cutoffs <- c(
      nFeature_meanv - sd_filtx * nFeature_sdv,
      nFeature_meanv + sd_filtx * nFeature_sdv
  )

  if (mt_pct < 0) {
    mt_meanv <- mean(seu$percent.mt)
    mt_sdv <- sd(seu$percent.mt)
    mt_pct <- min(c(mt_meanv + mt_sdv, 25.))
  }

  seu <- seu[, seu@meta.data$nCount_RNA   >= nc_cutoffs[[1]] &
               seu@meta.data$nCount_RNA   <= nc_cutoffs[[2]] &
               seu@meta.data$nFeature_RNA >= nf_cutoffs[[1]] &
               seu@meta.data$nFeature_RNA <= nf_cutoffs[[2]] &
               seu@meta.data$percent.mt   < mt_pct]
  return(seu)
}

# load in GSE117498 data

files <- list.files(here("GSE117498"), full.names = TRUE)

file_list <- map(files, ~ read_tsv(.x, col_types = cols()))
file_list <- map2(
    file_list, 
    1:length(file_list), 
    ~ set_colnames(.x, paste(colnames(.x), .y, sep = "_"))
)
md <- map(file_list, colnames)
names(md) <- files %>%
    basename() %>%
    str_remove(".raw_counts.tsv.gz") %>%
    str_remove("GSM[:digit:]+_") 

meta <- md %>%
    enframe() %>%
    unnest(cols = c(value)) %>%
    as.data.frame() %>%
    column_to_rownames(var = "value")

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
    assays = c("RNA"),
    verbose = FALSE
    )

good_features <- intersect(rownames(seurat), rownames(gse_seurat))
seurat <- seurat[good_features, ]

gse_seurat <- gse_seurat[good_features, ]

# Batch Correction
print("SPLITTING OBJECT -------")

SetPlan(resources_small)
seurat_ls <- SplitObject(seurat, split.by = "patient_id")

seurat_ls[["GSE"]] <- gse_seurat

seurat <- merge(seurat_ls[[1]], seurat_ls[2:length(seurat_ls)])
seurat
seurat <- subset(seurat, downsample = 50000)

plan(sequential)
seurat <- cache_rds(
    SCTransform(seurat, vst.flavor = "v2", conserve.memory = TRUE)
)

SetPlan(resources_small)
seurat <- cache_rds({
    seurat <- ScaleData(seurat)
    seurat <- FindVariableFeatures(seurat)
    seurat <- RunPCA(seurat)
    seurat <- RunUMAP(seurat, dims = 1:30)
    seurat <- FindNeighbors(seurat, dims = 1:30)
    FindClusters(seurat)
    },
    file = "int_dim_red.rds"
)
pdf("umap_integration.pdf")
UMAPPlot(seurat, group.by = "patient_id")
graphics.off()

seurat <- PrepSCTFindMarkers(seurat)

seurat[["stem_cell_type"]] <- seurat[[c("name", "stemness")]] %>%
    rownames_to_column() %>%
    unite(stem_cell_type, name, stemness) %>%
    column_to_rownames()

Idents(seurat) <- "stemness"
stem <- WhichCells(seurat, idents = 'Stem')
nonstem <- WhichCells(seurat, idents = 'Nonstem')
Idents(seurat) <- "name"

seurat <- SetIdent(seurat, cells = stem,        value = 'AML_Stem')
seurat <- SetIdent(seurat, cells = nonstem, value = 'AML_Nonstem')

markers <- FindAllMarkers(
    seurat,
    features = intersect(rownames(seurat), hbp_genes),
    test.use = "MAST",
    logfc.threshold = 0,
    assay = "SCT"
)
write_tsv(markers, file = "sct_integration_hbp_markers.tsv")

ggplot(data = markers, mapping = aes(x = cluster, y = avg_log2FC)) +
geom_col() +
facet_wrap(~ gene) +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file = "integration_markers_facet_gene.pdf")

ggplot(data = markers, mapping = aes(x = gene, y = avg_log2FC)) +
geom_col() +
facet_wrap(~ cluster) +
theme_classic() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(file = "integration_markers_facet_cluster.pdf")