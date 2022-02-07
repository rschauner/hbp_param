library(Seurat)
library(SeuratDisk)
library(tidyverse)

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                       assays = c("SCT"))

metadata <- seurat[[]]

metadata %>%
    group_by(stemness) %>%
    summarise(
        cells = n(),
        pats  = n_distinct(patient_id),
        rel   = n_distinct(patient_id, timepoint) - pats
    )

#   stemness cells  pats   rel
#   <chr>    <int> <int> <int>
# 1 Nonstem  97755    25    13
# 2 Stem     99200    20    11