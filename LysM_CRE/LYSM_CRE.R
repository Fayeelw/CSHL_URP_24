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
LYSM_C.data <- Read10X(data.dir="/Users/fayelynchwilliams/Desktop/cshl_project_24/LysM_CRE_data/control/") # data is expected in the .gz format
LYSM_VP16.data <- Read10X(data.dir="/Users/fayelynchwilliams/Desktop/cshl_project_24/LysM_CRE_data/vp16/") # data is expected in the .gz format

# Create Seurat objects

LYSM_C <- CreateSeuratObject(counts = LYSM_C.data, project = "LYSM_C", min.cells = 3, min.features = 200)
LYSM_VP16 <- CreateSeuratObject(counts = LYSM_VP16.data, project = "LYSM_C", min.cells = 3, min.features = 200)

# Merge LYSM_C and LYSM_VP16
merged <-  merge(LYSM_C, y = LYSM_VP16, add.cell.ids = c("LYSM_C", "LYSM_VP16")) # ids to identify as control v. vp16 (useful later)

# Pre-processing of nCount_RNA, nFeature_RNA, percent.mt
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "mt-")
vln_plot1 <- VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

# Subset after eyeballing the violin plot to ensure multiple/0 cells are removed
merged <- subset(merged, subset = nCount_RNA < 20000 & nCount_RNA > 100 & nFeature_RNA < 6000 & nFeature_RNA > 100 & percent.mt < 10)
vln_plot2 <- VlnPlot(merged, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

ggsave("distribution_vlnplot.png", plot = vln_plot2, dpi = 300)








