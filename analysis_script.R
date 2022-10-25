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
library(xfun)

plan("multicore")

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
  #seurat <- RunHarmony(seurat, c("seq_batch", "sort_batch"))
  seurat <- RunUMAP(seurat, dims = 1:20)
  seurat <- FindNeighbors(seurat, dims = 1:20)
  seurat <- FindClusters(seurat, resolution = 0.6)
  },
  file = "total_dim_red.rds"
)

message("Doing Differential Expression ------")
seurat_sub <- subset(seurat, downsample = 1000)

FindMarkers <- function(...) {
  res <- Seurat::FindMarkers(...)
  res <- rownames_to_column(res, var = "gene")
  return(res)
}

seurat_sub <- subset(seurat_sub, subset = paired == "T")
markers <- FindMarkers(
  seurat_sub, 
  test.use = "MAST",
  features = hbp_genes, 
  group.by = "timepoint", 
  ident.1 = "Diagnosis", 
  ident.2 = "Relapse",
  logfc.threshold = 0
)

write_tsv(markers, file = here("results/timepoint_DE_results_SCT_HBP.tsv"))

seurat_split <- SplitObject(seurat_sub, split.by = "patient_id")
seurat_split <- keep(seurat_split, ~nrow(unique(.x[["timepoint"]])) == 2)
markers <- map_dfr(seurat_split,
  FindMarkers,
  test.use = "MAST",
  features = hbp_genes,
  group.by = "timepoint",
  ident.1 = "Diagnosis",
  ident.2 = "Relapse",
  logfc.threshold = 0,
  .id = "patient_id"
)

write_tsv(markers, file = here("results/timepoint_DE_results_SCT_HBP_per_pat.tsv"))

pdf(file = here("plots/violin_plot.pdf"), width = 10, height = 10)
for (i in c("timepoint", "stemness", "prognosis")) {
  VlnPlot(
    seurat,
    slot = "scale.data",
    features = hbp_genes,
    group.by = i,
    split.by = "paitent_id"
  )
}
graphics.off()

message("Subsetting for Plotting ------")
subset_seurat <- list(
  paired_stem = subset(seurat, subset = stemness == "Stem" & paired == "T"),
  paired_nonstem = subset(seurat, subset = stemness == "Nonstem" & paired == "T"),
  diagnosis_stem = subset(seurat, subset = stemness == "Stem" & timepoint == "Diagnosis"),
  diagnosis_nonstem = subset(seurat, subset = stemness == "Nonstem" & timepoint == "Diagnosis"),
  diagnosis_only = subset(seurat, subset = timepoint == "Diagnosis"),
  paired_only = subset(seurat, subset = paired == "T")
)

message("Doing Printing Plots ------")

if (!dir.exists(here("plots"))) dir.create(here("plots"))

## Feature Plots ----
a <- FeaturePlot(
  subset_seurat[["paired_stem"]],
  features = hbp_genes,
  split.by = "timepoint") +
  plot_annotation(title = "Stem | Paired Only | Feature Plot By Timepoint") +
  plot_layout(ncol = 4)

b <- FeaturePlot(
  subset_seurat[["paired_nonstem"]],
  features = hbp_genes) +
plot_annotation(title = "Nonstem | Paired Only | Feature Plot By Timepoint")

c <- FeaturePlot(
  subset_seurat[["diagnosis_stem"]],
  features = hbp_genes,
  split.by = "prognosis") +
plot_annotation(title = "Stem | Diagnosis Only | Feature Plot By Prognosis") +
plot_layout(ncol = 4)

d <- FeaturePlot(
  subset_seurat[["diagnosis_nonstem"]],
  features = hbp_genes) +
plot_annotation(title = "Nonstem | Diagnosis Only | Feature Plot of HBP Genes")

e <- DotPlot(seurat, features = hbp_genes, cluster.idents = TRUE) +
  plot_annotation(title = "Stem + Nonstem | All Samples | Dot Plot By Cluster")

f <- DotPlot(subset_seurat[["paired_only"]], features = hbp_genes, group.by = "timepoint", scale = FALSE) +
NoLegend() +
FontSize(x.text = 0, x.title = 0, y.title = 0)
g <- DotPlot(subset_seurat[["diagnosis_only"]], features = hbp_genes, group.by = "prognosis", scale = FALSE) +
FontSize(x.text = 0, x.title = 0, y.title = 0)
h <- DotPlot(subset_seurat[["diagnosis_only"]], features = hbp_genes, group.by = "stemness", scale = FALSE) +
NoLegend() +
FontSize(y.title = 0)

p <- f / g / h +
  plot_annotation(title = "Stem + Nonstem | All Samples | Dot Plot By Timepoint/Prognosis/Stemness")


md<- mutate(
  seurat[[]],
  pat_time = paste0(patient_id, "_", timepoint)
)
seurat[["pat_time"]] <- md[["pat_time"]]

n <- DotPlot(seurat, features = c("OGA", "OGT"), group.by = "pat_time") +
  plot_annotation(title = "Stem + Nonstem | All Samples | Dot Plot by Timepoint and Patient")

pdf(file = here("plots/feature_plots_SCT.pdf"), width = 10, height = 10)
print(a)
print(b)
print(c)
print(d)
print(e)
graphics.off()

pdf(file = here("plots/dot_plots_SCT.pdf"), width = 10, height = 10)
print(p)
print(n)
graphics.off()


BlendScatter <- function(object) {
  FeaturePlot(
    object,
    features = c("OGA", "OGT"),
    blend = TRUE,
    blend.threshold = 0.1,
    coord.fixed = TRUE,
    order = TRUE,
    ncol = 2
  )
}

subset_seurat <- map(subset_seurat, function(x) {
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:30)
  return(x)
})

pdf(file = here("plots/blend_OGA_OGT.pdf"), width = 10, height = 2.5)
for (i in seq_along(subset_seurat)) {
  BlendScatter(subset_seurat[[i]]) + labs(title = names(subset_seurat)[[i]])
}
graphics.off()
