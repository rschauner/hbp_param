
message("Making Feature Plots ------")

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

pdf(file = here("plots/feature_plots_SCT.pdf"), width = 10, height = 10)
print(a)
print(b)
print(c)
print(d)
graphics.off()
