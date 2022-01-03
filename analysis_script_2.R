library(Seurat)
library(here)
library(readxl)
library(SeuratDisk)
library(patchwork)
library(msigdbr)
library(rsinglecell)
library(limma)
library(EnhancedVolcano)
library(future)
library(tidyverse)
library(xfun)
library(furrr)

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                       assays = c("SCT"))

DefaultAssay(seurat) <- "SCT"

hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UBAP1", "OGT", "MGEA5", "WNK1", "REL")

objects <- list(stem = subset(seurat, subset = stemness == "Stem"))
objects[["nonstem"]] = subset(seurat, subset = stemness == "Nonstem")

#objects <- map(objects, ~ScaleData(.x, features = c(VariableFeatures(.x), hbp_genes))

plan(tweak("multicore", workers = 2))
markers <- future_map(
    objects,
    ~ FindMarkers(
        .x, 
        test.use = "MAST", 
        features = hbp_genes, 
        group.by = "timepoint", 
        ident.1 = "Diagnosis", 
        ident.2 = "Relapse",
        logfc.threshold = 0
    )
)
markers <- map(markers, rownames_to_column, var = "gene")
markers <- tibble(data = markers, stemness = names(objects))
markers <- unnest(markers, cols = c(data))

library(viridis)

pdf(here("plots/DE_DvR_facet_stemness.pdf"))
ggplot(
    data = markers,
    mapping = aes(x = stemness, y = avg_log2FC, fill = abs(pct.1 - pct.2))
) +
facet_wrap( ~ gene) +
geom_col() + 
scale_fill_viridis() +
theme_classic()
graphics.off()