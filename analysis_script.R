library(Seurat)
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

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                       assays = c("SCT"))

DefaultAssay(seurat) <- "SCT"

if (DefaultAssay(seurat) == "RNA") {
  seurat <- NormalizeData(seurat)
  seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
}

hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UBAP1", "OGT", "MGEA5", "WNK1", "REL")

seurat <- cache_rds(
  DoDimensionReductions(seurat, batch_vars = c("seq_batch", "sort_batch")),
  filename = "total_dim_red.rds")
seurat <- ScaleData(seurat, features = c(VariableFeatures(seurat), hbp_genes))

plan("multiprocess", workers = 16)

seurat_sub <- subset(seurat, downsample = 1000)

#markers <- FindAllMarkers(seurat_sub, test.use = "MAST")
#write_tsv(markers, file = here("cluster_DE_results_SCT.tsv"))

markers <- FindAllMarkers(seurat_sub, test.use = "MAST", features = hbp_genes)
write_tsv(markers, file = here("cluster_DE_results_SCT_HBP.tsv"))

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
markers <- rownames_to_column(markers, var = "gene")
write_tsv(markers, file = here("timepoint_DE_results_SCT_HBP.tsv"))

if (!dir.exists(here("plots"))) dir.create(here("plots"))

## Feature Plots ----
pdf(file = here("plots/feature_plots_SCT.pdf"), width = 10, height = 10)
a <- FeaturePlot(subset(seurat, subset = stemness == "Stem" & paired == "T"),
                 features = hbp_genes,
                 split.by = "timepoint") +
  plot_annotation(title = "Stem | Paired Only | Feature Plot By Timepoint") +
  plot_layout(ncol = 4)

b <- FeaturePlot(subset(seurat, subset = stemness == "Nonstem" & paired == "T"),
                 features = hbp_genes) +
  plot_annotation(title = "Nonstem | Paired Only | Feature Plot By Timepoint")

c <- FeaturePlot(subset(seurat, subset = stemness == "Stem" & timepoint == "Diagnosis"),
                 features = hbp_genes,
                 split.by = "prognosis") +
  plot_annotation(title = "Stem | Diagnosis Only | Feature Plot By Prognosis") +
  plot_layout(ncol = 4)

d <- FeaturePlot(subset(seurat, subset = stemness == "Nonstem" & timepoint == "Diagnosis"),
                 features = hbp_genes) +
  plot_annotation(title = "Nonstem | Diagnosis Only | Feature Plot of HBP Genes")

e <- DotPlot(seurat, features = hbp_genes, cluster.idents = TRUE) +
  plot_annotation(title = "Stem + Nonstem | All Samples | Dot Plot By Cluster")

f <- DotPlot(seurat, features = hbp_genes, group.by = "timepoint") +
  NoLegend() +
  FontSize(x.text = 0, x.title = 0, y.title = 0)
g <- DotPlot(seurat, features = hbp_genes, group.by = "prognosis") +
  FontSize(x.text = 0, x.title = 0, y.title = 0)
h <- DotPlot(seurat, features = hbp_genes, group.by = "stemness") +
  NoLegend() +
  FontSize(y.title = 0)

p <- f / g / h +
  plot_annotation(title = "Stem + Nonstem | All Samples | Dot Plot By Timepoint/Prognosis/Stemness")

print(a)
print(b)
print(c)
print(d)
print(e)
print(p)

VlnPlot(seurat, c("MGEA5", "OGT"), group.by = "patient_id", split.by = "timepoint") +
  plot_annotation(title = "Stem + Nonstem | All Samples | Violin Plot by Timepoint and Patient")

graphics.off()

## GSVA based on MGEA5 expression ----
Idents(seurat) <- "MGEA5_neg"
Idents(seurat, cells = WhichCells(seurat, expression = MGEA5 > 0.5)) <- "MGEA5_lo"
Idents(seurat, cells = WhichCells(seurat, expression = MGEA5 > 2)) <- "MGEA5_hi"

gene_sets <- list(msigdbr(species = "Homo sapiens", category = "H"),
                  msigdbr(species = "Homo sapiens", category = "C7")) %>%
  bind_rows() %>%
  ListGeneSets(gene_set = "name", gene_name = "symbol")

gsva_res <- RunGSVA(seurat, gene_sets = gene_sets, average = TRUE, replicates = 3)

gsva_res %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_set") %>%
  write_tsv(file = here("MGEA5_GSVA_SCT.tsv"))

gsva_stats <- RunStats(gsva_res, "msigdbr")

pdf(here("plots/volcano_SCT.pdf"))

EnhancedVolcano(topTable(fit, coef = "HIvNEG", n = Inf),
                x = "logFC",
                y = "adj.P.Val",
                lab = rownames(topTable(fit, coef = "HIvNEG", n = Inf)),
                FCcutoff = 0.5,
                pCutoff = 0.01,
                title = "MGEA5 hi v neg")

EnhancedVolcano(topTable(fit, coef = "HIvLO", n = Inf),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.01,
                lab = rownames(topTable(fit, coef = "HIvLO", n = Inf)),
                FCcutoff = 0.5,
                title = "MGEA5 hi v lo")

EnhancedVolcano(topTable(fit, coef = "LOvNEG", n = Inf),
                x = "logFC",
                y = "adj.P.Val",
                lab = rownames(topTable(fit, coef = "LOvNEG", n = Inf)),
                selectLab = rownames(topTable(fit, coef = "LOvNEG", n = 5)),
                FCcutoff = 0.5,
                labSize = 2.5,
                labhjust = 2,
                pCutoff = 0.01,
                drawConnectors = TRUE,
                arrowheads = TRUE,
                xlim = c(-1, 1),
                ylim = c(0, 7.5),
                title = "MGEA5 lo v neg")

graphics.off()

## Differential Expression ----
seurat_sub <- subset(seurat, downsample = 1000)
#markers <- FindAllMarkers(seurat_sub, test.use = "MAST")
#write_tsv(markers, file = here("MGEA5_DE_results_SCT.tsv"))
markers <- FindAllMarkers(seurat_sub, test.use = "MAST", features = hbp_genes)
write_tsv(markers, file = here("MGEA5_DE_results_SCT_HBP.tsv"))
