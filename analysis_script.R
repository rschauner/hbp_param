library(Seurat)
library(here)
library(readxl)
library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(msigdbr)

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                       assays = "SCT")

hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UBAP1", "OGT", "MGEA5", "WNK1")
seurat <- ScaleData(seurat, features = c(VariableFeatures(seurat), hbp_genes))

if (!dir.exists(here("plots/HBP"))) dir.create(here("plots/HBP"))

pdf(file = here("plots/HBP/plots-%01d.pdf"), width = 10, height = 10, onefile = FALSE)
a <- FeaturePlot(subset(seurat, subset = stemness == "Stem" & paired == "T"),
                 features = hbp_genes,
                 split.by = "timepoint") +
  plot_annotation(title = "Stem | Paired Only | Feature Plot By Timepoint")
b <- FeaturePlot(subset(seurat, subset = stemness == "Nonstem" & paired == "T"),
                 features = hbp_genes) +
  plot_annotation(title = "Nonstem | Paired Only | Feature Plot By Timepoint")

c <- FeaturePlot(subset(seurat, subset = stemness == "Stem" & timepoint == "Diagnosis"),
                 features = hbp_genes,
                 split.by = "prognosis") +
  plot_annotation(title = "Stem | Diagnosis Only | Feature Plot By Prognosis")
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

# Get MGEA5 hi cells and do GSVA in HALLMARK and C7

gene_sets <-list(msigdbr(species = "Homo sapiens", category = "H"),
                 msigdbr(species = "Homo sapiens", category = "C7")) %>%
  bind_rows() %>%
  ListGeneSets(gene_set = "name", gene_name = "symbol")

gsva_res <- RunGSVA(seurat, gene_sets = gene_sets, average = FALSE)

# compare MGEA5 and OGT expression to HSCs


