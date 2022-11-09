library(xfun)
library(pheatmap)
library(Seurat)
library(future)
library(viridisLite)
library(tidyverse)

plan("multicore")

with_plan <- function(expr, ...) {
    oplan <- plan()
    on.exit(plan(oplan))
    plan(...)
    eval(expr)
}

files <- list.files(
    path = "data/flow_data",
    pattern = "*.csv$",
    full.names = TRUE,
    recursive = FALSE
)

flowLS <- lapply(
    files,
    function(FILE) {
        fileDF <- read.csv(FILE) %>%
            as.data.frame() %>%
            mutate(DonorID = basename(FILE))
        return(fileDF)
    }
)

flowDF <- do.call(bind_rows, flowLS)
flowDF <- mutate(flowDF, DonorID = str_remove(DonorID, pattern = ".csv"))
donor_id <- flowDF$DonorID
flowDF$DonorID <- flowDF$FSC.A <- flowDF$FSC.H <- flowDF$FSC.W <-
    flowDF$SSC.A <- flowDF$SSC.H <- flowDF$SSC.W <-
        flowDF$Viability.Fixable.Near.IR.A <- flowDF$Time <- NULL

colnames(flowDF) <- c("O-GlcNAc", "CD34", "CD117", "CD38")
rownames(flowDF) <- seq(nrow(flowDF))

assay <- CreateAssayObject(counts = t(as.matrix(flowDF)))
assay@data <- t(as.matrix(flowDF))
assay@scale.data <- t(as.matrix(scale(flowDF, scale = TRUE, center = TRUE)))
seurat <- CreateSeuratObject(assay, assay = "FCM")

seurat[["donor_id"]] <- donor_id
seurat[["status"]] <- if_else(sapply(donor_id, nchar) == 4, "AML", "Healthy")

seurat <- subset(seurat, subset = status == "AML")
seurat <- RunPCA(seurat, features = rownames(seurat))
cache_rds(
    expr = seurat <- FindNeighbors(seurat, dims = 1:3, k = 60),
    file = "nn.rds",
    hash = list(seurat)
)
cache_rds(
    with_plan(
        seurat <- FindClusters(seurat, resolution = 0.15),
        "sequential"
    ),
    file = "cluster.rds",
    hash = list(seurat)
)

length(unique(Idents(seurat)))
avg_exp <- AverageExpression(seurat, slot = "counts")[[1]]

pheatmap(
    avg_exp,
    color = magma(100),
    cellwidth = 10,
    cellheight = 10
)