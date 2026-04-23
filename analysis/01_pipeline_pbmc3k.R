# 01_pipeline_pbmc3k.R
# End-to-end scRNA-seq mini pipeline (PBMC3k) using Seurat (R)
# Outputs: QC plot, PCA elbow, UMAP clusters, canonical marker FeaturePlot, markers table

set.seed(123)

# ---- libraries ----
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(dplyr)
})

# ---- outputs ----
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# ---- load PBMC3k (public SeuratData dataset) ----
# Requires SeuratData installed from GitHub (satijalab/seurat-data)
data("pbmc3k")

# Update object to match local Seurat version (prevents slot/version issues)
pbmc <- UpdateSeuratObject(pbmc3k)

# ---- QC ----
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

png("output/figures/01_qc_violin.png", width = 1400, height = 600, res = 150)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# Filtering (typical PBMC tutorial thresholds)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# ---- Normalize + features + PCA ----
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

png("output/figures/02_pca_elbow.png", width = 900, height = 700, res = 150)
ElbowPlot(pbmc)
dev.off()

# ---- Clustering + UMAP ----
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunUMAP(pbmc, dims = 1:10)

png("output/figures/03_umap_clusters.png", width = 900, height = 700, res = 150)
DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

# ---- Markers ----
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "output/tables/markers_all_clusters.csv", row.names = FALSE)

# Canonical markers (quick sanity check)
canonical <- c("IL7R","CCR7","NKG7","GNLY","MS4A1","CD79A","LYZ","S100A8","FCGR3A","MS4A7","PPBP")
png("output/figures/04_featureplot_markers.png", width = 1600, height = 900, res = 150)
FeaturePlot(pbmc, features = canonical, ncol = 4)
dev.off()

# Save Seurat object for next script
saveRDS(pbmc, file = "output/pbmc_processed.rds")

cat("DONE.\nFigures:", list.files("output/figures"), "\n")
cat("Tables:", list.files("output/tables"), "\n")