# 02_annotation_signature.R
# Cluster annotation + cytotoxic program scoring (portfolio-friendly)
# Inputs: output/pbmc_processed.rds and output/tables/markers_all_clusters.csv
# Outputs: annotated UMAP, cytotoxic score UMAP, cytotoxic score table

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# Load processed object + markers
pbmc <- readRDS("output/pbmc_processed.rds")
markers <- read.csv("output/tables/markers_all_clusters.csv")

# Make sure identities are numeric clusters
Idents(pbmc) <- "seurat_clusters"

# Apply labels based on your observed marker patterns (PBMC3k typical)
pbmc <- RenameIdents(pbmc,
  `0`="CD4 T (CCR7/LEF1, naive/memory)",
  `1`="CD14+ monocytes (S100A8/A9)",
  `2`="Activated T (CD40LG+)",
  `3`="B cells (TCL1A/CD79A)",
  `4`="Cytotoxic T (CD8A/GZMK/CCL5)",
  `5`="Dendritic cells (BATF3+)",
  `6`="NK / cytotoxic (FGFBP2/GZMB)",
  `7`="Ambiguous / low-signal",
  `8`="Dendritic cells (FCER1A/CLEC10A)",
  `9`="Platelets (ITGA2B/GP9)"
)

png("output/figures/05_umap_celltypes.png", width = 900, height = 700, res = 150)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# Cytotoxic program score (example of MoA/biomarker-style scoring)
cytotoxic_genes <- c("NKG7","GNLY","PRF1","GZMB","GZMH","FGFBP2")
pbmc <- AddModuleScore(pbmc, features = list(cytotoxic_genes), name = "CytotoxicScore")

png("output/figures/06_umap_cytotoxic_score.png", width = 900, height = 700, res = 150)
FeaturePlot(pbmc, features = "CytotoxicScore1")
dev.off()

# Table: average score per annotated cell type
score_tab <- data.frame(
  celltype = Idents(pbmc),
  score = pbmc$CytotoxicScore1
) |>
  group_by(celltype) |>
  summarize(mean_score = mean(score), .groups = "drop") |>
  arrange(desc(mean_score))

write.csv(score_tab, "output/tables/cytotoxic_score_by_celltype.csv", row.names = FALSE)

cat("DONE.\nSaved annotated UMAP and cytotoxic score outputs.\n")
print(score_tab)