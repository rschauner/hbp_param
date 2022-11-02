message("Making Blended Feature Plots ------")

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

pdf(file = here("plots/blend_OGA_OGT.pdf"), width = 10, height = 2.5)
for (i in seq_along(subset_seurat)) {
    print(BlendScatter(subset_seurat[[i]]) + labs(title = names(subset_seurat)[[i]]))
}
graphics.off()
