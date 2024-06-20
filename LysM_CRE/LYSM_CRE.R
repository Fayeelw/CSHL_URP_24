options(Seurat.object.assay.version = "v3")
library(dplyr)
library(Seurat)
library(patchwork) 
install.packages("R.utils")
library(R.utils)
library(sctransform)
library(ggplot2)
install.packages("harmony")
library(harmony)
install.packages("MAST")
library(MAST)

# Load data
LYSM_C.data <- Read10X(data.dir="/Users/fayelynchwilliams/Desktop/CSHL_URP_24/LysM_CRE_data/control/") # data is expected in the .gz format
LYSM_VP16.data <- Read10X(data.dir="/Users/fayelynchwilliams/Desktop/CSHL_URP_24/LysM_CRE_data/vp16/") # data is expected in the .gz format

# Create Seurat objects

LYSM_C <- CreateSeuratObject(counts = LYSM_C.data, project = "LYSM_C", min.cells = 3, min.features = 200)
LYSM_VP16 <- CreateSeuratObject(counts = LYSM_VP16.data, project = "LYSM_VP16", min.cells = 3, min.features = 200)

# Merge LYSM_C and LYSM_VP16
merged <-  merge(LYSM_C, y = LYSM_VP16, add.cell.ids = c("LYSM_C", "LYSM_VP16")) # ids to identify as control v. vp16 (useful later)

# Pre-processing of nCount_RNA, nFeature_RNA, percent.mt
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "mt-")
vln_plot1 <- VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

# Subset after eyeballing the violin plot to ensure multiple/0 cells are removed
merged <- subset(merged, subset = nCount_RNA < 20000 & nCount_RNA > 100 & nFeature_RNA < 6000 & nFeature_RNA > 100 & percent.mt < 10)
vln_plot2 <- VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

ggsave("distribution_vlnplot.png", plot = vln_plot2, dpi = 300)

# SCTransform
merged <- SCTransform(merged, vars.to.regress = "percent.mt", verbose = FALSE)

# Dimensionality reduction
merged <- RunPCA(merged, verbose = FALSE)

# Batch correction - harmony\
merged$treatment <- factor(merged$orig.ident)
merged <- merged %>%
  RunHarmony(group.by.vars =  "treatment", assay.use="SCT", plot_convergence = TRUE)

# Run UMAP
elbowplot1 <- ElbowPlot(merged, ndims = 50) # Check the inflection point
ggsave("elbowplot.png", plot = elbowplot1, dpi = 300)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE, reduction='harmony')
DimPlot(merged, label = TRUE)

# Cluster and make another DimPlot() to visualise subsets
merged <- FindNeighbors(merged, dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, verbose = FALSE)
dimplot1 <- DimPlot(merged, label = TRUE)
ggsave("cluster_dimplot.png", plot = dimplot1, dpi = 300, width = 10, height = 9)

# Find DE markers for every cluster compared to all remaining cells - report only positive markers
merged <- PrepSCTFindMarkers(merged)
merged.markers <- FindAllMarkers(merged, only.pos = TRUE)

# Arrange DE cluster markers in descending order based on gene expression levels
cluster_markers <- merged.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(cluster_markers, "cluster_markers.csv")

# Use a publicly available list of MOUSE marker genes to manually annotate clusters
cluster0 <- cluster_markers[cluster_markers$cluster == "0", ] # M2 macrophage (Nos2, Arg1)
print(cluster0, n = 100)
cluster1 <- cluster_markers[cluster_markers$cluster == "1", ] # CD8 T cell (Gzmk, Tcrg-C2, Cd8b1/Cd8a, Themis)
print(cluster1, n = 100)
cluster2 <- cluster_markers[cluster_markers$cluster == "2", ] # CD4 T cell (Cd4, Cd40Ig, Il4)
print(cluster2, n = 100)
cluster3 <- cluster_markers[cluster_markers$cluster == "3", ] # CD8 T cell 
print(cluster3, n = 100)
cluster4 <- cluster_markers[cluster_markers$cluster == "4", ] # CD8 T cell
print(cluster4, n = 100)
cluster5 <- cluster_markers[cluster_markers$cluster == "5", ] # Treg 
print(cluster5, n = 100)
cluster6 <- cluster_markers[cluster_markers$cluster == "6", ] # NK cell 
print(cluster6, n = 100)
cluster7 <- cluster_markers[cluster_markers$cluster == "7", ] # M1 Macrophage
print(cluster7, n = 100)
cluster8 <- cluster_markers[cluster_markers$cluster == "8", ] # DC (Cd74 + H2- = Ag presentation; Mgl2 - Ag recognition; Cd209a/e)
print(cluster8, n = 100)
cluster9 <- cluster_markers[cluster_markers$cluster == "9", ] # Neutrophil 
print(cluster9, n = 100)
cluster10 <- cluster_markers[cluster_markers$cluster == "10", ] # NK cell
print(cluster10, n = 100)
cluster11 <- cluster_markers[cluster_markers$cluster == "11", ] # M2 macrophage
print(cluster11, n = 100)
cluster12 <- cluster_markers[cluster_markers$cluster == "12", ] # Monocyte 
print(cluster12, n = 100)
cluster13 <- cluster_markers[cluster_markers$cluster == "13", ] # CD8 T cell  
print(cluster13, n = 100)
cluster14 <- cluster_markers[cluster_markers$cluster == "14", ] # Immature B cell 
print(cluster14, n = 100)
cluster15 <- cluster_markers[cluster_markers$cluster == "15", ] # Th17 
print(cluster15, n = 100)
cluster16 <- cluster_markers[cluster_markers$cluster == "16", ] # Proliferating cell
print(cluster16, n = 100)
cluster17 <- cluster_markers[cluster_markers$cluster == "17", ] # cDC1
print(cluster17, n = 100)
cluster18 <- cluster_markers[cluster_markers$cluster == "18", ] #pDC
print(cluster18, n = 100)
cluster19 <- cluster_markers[cluster_markers$cluster == "19", ] # Plasma cell 
print(cluster19, n = 100)

# Dotplot to check markers
dotplot1 <- DotPlot(merged, features = c("Cd68", "H2-Ab1", "Cd74", "Cd3e", "Cd8a", "Cd4", "Klrg1", "Prf1", "Eomes", "Arg1"))
ggsave("cluster_dotplot.png", plot = dotplot1, width = 10, height = 9)

# Define cluster names based on marker genes 
cluster_annotation <- c("M2 macrophage", "CD8 T cell", "CD4 T cell", "CD8 T cell", "CD8 T cell", "Treg", "NK cell", "M1 macrophage", "DC", "Neutrophil", "NK cell", "M2 macrophage", "Monocyte", "CD8 T cell", "Immature B cell", "Th17", "Proliferating cell", "cDC1", "pDC", "Plasma cell")
names(cluster_annotation) <- levels(merged)
merged <- RenameIdents(merged, cluster_annotation)
merged$cell_type <- Idents(merged) ############!!!!!!!!

# Visualise annotations with a DimPlot()
dimplot2 <- DimPlot(merged, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("annotated_dimplot.png", plot = dimplot2, width = 10, height = 9)

# The above is comparing clusters to one another - now need to perform DE analysis between WT and VP16

# Set the cell identity/classification equal to the metadata column orig.ident - this identity class allows VP16/C classes to be distinguished
Idents(merged) <- "treatment"

# Find DE markers between LYSM_VP16 and LYSM_C across all clusters

de_genes <- FindMarkers(merged, ident.1 = "LYSM_VP16", ident.2 = "LYSM_C", min.pct = 0.25, logfc.thresold = 0.25) 
de_genes <- de_genes %>% arrange(desc(avg_log2FC))

# Find DE markers between LYSM_VP16 and LYSM_C within each cluster
deg_list <- list()

# Loop through each cluster
for (cell_type in levels(merged$cell_type)) {
  # When subsetting by cluster, ensure that the active identity class is set to 'cell_type' before subsetting. 
  Idents(merged) <- "cell_type"
  # Subset the seurat object to the current cell_type
  cell_subset <- subset(merged, idents = cell_type)
  
  # Restore the active identity to 'treatment' for DE analysis
  Idents(cell_subset) <- "treatment"
  # Perform differential expression analysis between LYSM_VP16 and LYSM_C
  degs <- FindMarkers(
    object = cell_subset,
    ident.1 = "LYSM_VP16",
    ident.2 = "LYSM_C",
    group.by = "treatment",
    test.use = "MAST"
  )
  
  # Store the DEGs in the list with the cluster ID as the name 
  deg_list[[cell_type]] <- degs
}
print(deg_list)

# Export to csv
p_val_threshold <- 0.05
for (cell_type in names(deg_list)) {       
  significant_degs <- deg_list[[cell_type]][deg_list[[cell_type]]$p_val < p_val_threshold, ]
  write.csv(significant_degs, file = paste0("DEGs_cluster_", cell_type, ".csv"))
}

# Manual gene level analysis - use the csv files to make plots to deduce different gene expression changes 
DEGs_cluster_CD8T <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_CD8 T cell.csv")
DEGs_cluster_cDC1 <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_cDC1.csv")
DEGs_cluster_DC <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_DC.csv")
DEGs_cluster_developing_CD8 <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Developing CD8 T cell.csv")
DEGs_cluster_immature_B <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Immature B cell.csv")
DEGs_cluster_M2_macrophage <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_M2 Macrophage.csv")
DEGs_cluster_monocyte <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Monocyte.csv")
DEGs_cluster_neutrophil <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Neutrophil.csv")
DEGs_cluster_nk <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_NK Cell.csv")
DEGs_cluster_pDC <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_pDC.csv")
DEGs_cluster_plasma_cell <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Plasma Cell.csv")
DEGs_cluster_proliferating_cell <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Proliferating Cell.csv")
DEGs_cluster_T_cell <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_T cell.csv")
DEGs_cluster_Th2_cell <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Th2 CD4 T cell.csv")
DEGs_cluster_Th17_cell <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Th17 CD4 T cell.csv")
DEGs_cluster_treg <- read_csv("Treatment_DEGs_per_cell_type/DEGs_cluster_Treg.csv")

# Cell Chat/Nichenet - evaluate how cell-cell communication changes with treatment

