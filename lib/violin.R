
message("Making Violin Plots ------")
pdf(file = here("plots/violin_plot.pdf"), width = 10, height = 10)
for (i in c("timepoint", "stemness", "prognosis")) {
    print(VlnPlot(
        seurat,
        slot = "scale.data",
        features = hbp_genes,
        group.by = i,
        split.by = "patient_id"
    ))
}
graphics.off()
