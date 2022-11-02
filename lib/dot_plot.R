message("Making Dot Plots ------")
MakeDotPlot <- function(object, group, name, ...) {
    DotPlot(object, features = hbp_genes, group.by = group, scale = FALSE) +
    labs(caption = paste0("Grouped by ", group, " using ", name))
}


pdf(file = here("plots/dot_plots_SCT.pdf"), width = 10, height = 4)
print(pmap(list(subset_seurat, groups, names(subset_seurat)), MakeDotPlot))
graphics.off()
