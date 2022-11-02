message("Running Differential Expression Analysis ------")

FindMarkers <- function(...) {
    res <- Seurat::FindMarkers(...)
    res <- rownames_to_column(res, var = "gene")
    return(res)
}

RunMAST <- function(object, group, ...) {
    if (group == "timepoint") {
        ident1 <- "Diagnosis"
        ident2 <- "Relapse"
    } else if (group == "prognosis") {
        ident1 <- "Favorable"
        ident2 <- "Poor"
    }

    hbp_genes <- c("GFAP", "NAGK", "GNPNAT1", "PGM3", "UAP1", "OGT", "OGA", "GFPT1", "GFPT2")
    res <- FindMarkers(
        object = object,
        ident.1 = ident1,
        ident.2 = ident2,
        group.by = group,
        assay = "RNA",
        test.use = "MAST",
        logfc.threshold = 0,
        features = hbp_genes,
        ...
    )
    return(res)
}

RunPerPatientDE <- function(object, group) {
    print("splitting by patient")
    split <- SplitObject(object, split.by = "patient_id")
    print("keeping")
    split <- keep(split, ~ nrow(unique(.x[[group]])) == 2)
    print("running MAST")
    markers <- map_dfr(split,
        RunMAST,
        group = "timepoint",
        .id = "patient_id"
    )
    return(markers)
}

seurat_sub <- map(subset_seurat, subset, downsample = 1000)

# markers <- map2(seurat_sub, groups, RunMAST)
# write_xlsx(markers, path = here("results/bulk_DE.xlsx"))

message("Running Per Patient DE Analysis ------")

markers <- cache_rds(
    map2(seurat_sub, groups, RunPerPatientDE),
    file = "per_pat_de.rds"
)

write_xlsx(markers, path = here("results/per_patient_DE.xlsx"))
