library(Seurat)
library(Matrix)
library(irlba)  # For PCA
library(RcppAnnoy)  # For fast nearest neighbor search
library(dplyr)

# Assuming the PBMC datasets (3k and 10k) are already normalized
# and represented as sparse matrices
# devtools::install_github('satijalab/seurat-data')

library(SeuratData)
AvailableData()

# Set a longer timeout (e.g., 300 seconds or 5 minutes)
options(timeout = 300)

# Then try installing the data again
InstallData("pbmc3k")

pbmc3k<-UpdateSeuratObject(pbmc3k)
pbmc3k@meta.data %>% head()

# The PBMC3k dataset (peripheral blood mononuclear cells, ~3,000 cells) is a classic single-cell RNA-seq dataset that is commonly used for tutorials in Seurat (an R package for single-cell analysis). It typically comes as a Seurat object, which is an S4 object containing multiple slots that store different types of data.

# List all slots
slotNames(pbmc3k)

# View metadata
head(pbmc3k@meta.data)

#nCount_RNA: Total number of RNA molecules (UMIs) per cell.
#nFeature_RNA: Number of genes detected per cell.
#percent.mt: Percentage of reads mapping to mitochondrial genes.
#seurat_clusters: Cluster identity assigned by Seurat.

#What it contains: Raw and processed gene expression matrices.
pbmc3k@assays

pbmc3k[["RNA"]]@counts # Raw UMI counts (sparse matrix).https://help.geneiousbiologics.com/hc/en-us/articles/4781289585300-Understanding-Single-Cell-technologies-Barcodes-and-UMIs 
pbmc3k[["RNA"]]@data # Normalized expression (log-normalized by default).
pbmc3k[["RNA"]]@scale.data # Scaled (Z-scored) data, used for PCA.


# View expression data
pbmc3k[["RNA"]]@counts[1:5, 1:5]    # raw counts
pbmc3k[["RNA"]]@data[1:5, 1:5]      # normalized

# View clustering results
table(Idents(pbmc3k))

#What it contains: History of processing steps (used internally by Seurat).
pbmc3k@commands

#What it contains: Active identity class of each cell (e.g., clusters or cell types).
pbmc3k@active.ident


# Calculate mitochondrial gene percentage: High mitochondrial gene expression often indicates stressed or dying cells.
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")

#     PercentageFeatureSet(pbmc3k, pattern = "^MT-")
#This function scans all gene names in the dataset (pbmc3k) and selects those whose names start with "MT-" using the regular expression ^MT-.
#      pbmc3k[["percent.mt"]] <- ...
# This line adds the resulting vector to the @meta.data slot of the Seurat object as a new metadata column named "percent.mt".
# You can access it later using pbmc3k@meta.data$percent.mt.

#[["percent.mt"]]: This is how you access or add columns to the metadata slot of a Seurat object.
#Seurat objects store various information about each cell in a data frame accessible via pbmc3k@meta.data or, more conveniently, using the [[ operator.
#Cells with an excessively high percent.mt (e.g., >5% or >10%, depending on the dataset and experimental conditions) are typically discarded to ensure that only high-quality, viable cells are analyzed,

# Check first few values
head(pbmc3k@meta.data$percent.mt)

# Visualize mitochondrial percentage vs total counts
VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2 # Requires the 'patchwork' package for combining plots

#Recommended Filtering Step (Before NormalizeData):
pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#nFeature_RNA > 200: filters out empty droplets
#nFeature_RNA < 2500: filters out likely doublets
#percent.mt < 5: removes dying/stressed cells

# routine processing
pbmc3k<- pbmc3k %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  #Normalizes gene expression values per cell.
  #How: Counts are divided by the total expression in each cell, multiplied by 10,000 (scale factor), and then log-transformed.
  #Why: To make expression values comparable across cells (removes library size effects).
  
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  #What it does: Identifies the top 3,000 highly variable genes across cells.
  #Why: Variable genes are more informative for downstream analysis (e.g., clustering, PCA).
  #vst method: Variance stabilizing transformation, a popular method for selecting variable genes.
 
  ScaleData() %>%
  #What it does:Centers and scales gene expression data.
  #For each gene: subtracts the mean, divides by standard deviation.
  #Why: Required for PCA and clustering to ensure all genes contribute equally (removes effects of highly expressed genes).
  #Also allows regressing out unwanted variables (e.g., percent.mt, nCount_RNA) if specified.
  
  RunPCA(verbose = FALSE) %>%
  #What it does: Performs principal component analysis (PCA) on the scaled data.
  #Why: Reduces dimensionality of the dataset while preserving the most variance.
  #PCA results are used for clustering and UMAP.
  
  FindNeighbors(dims = 1:10, verbose = FALSE) %>%
  #What it does: Constructs a K-nearest neighbors (KNN) graph using the first 10 principal components.
  #Why: The graph is used as the foundation for clustering in the next step.
  
  FindClusters(resolution = 0.5, verbose = FALSE) %>%
  #What it does: Applies a graph-based community detection algorithm (Louvain or Leiden) to identify clusters of cells.
  #Resolution: Controls granularity.
  #Higher = more clusters (e.g., 0.8) Lower = fewer clusters (e.g., 0.4)
  #Output: Cluster labels stored in pbmc3k@meta.data$seurat_clusters
 
  RunUMAP(dims = 1:10, verbose = FALSE)
  #What it does: Computes a 2D UMAP (Uniform Manifold Approximation and Projection) for visualization.
  #Input: Uses first 10 PCs for dimensionality.
  #Output: pbmc3k@reductions$umap


#  Get an idea of how the pbmc3k data look like:
p1<- DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.5")
#RNA_snn_res.0.5, This column was automatically added to pbmc3k@meta.data when you ran:
#It stores the unsupervised cluster assignments from the Louvain/Leiden algorithm.
#It's a numeric or factor column like: 0, 1, 2, … — each value is a cluster label.
#It depends on the resolution you used (in this case, 0.5).

p2<- DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = "seurat_annotations", label.size = 3)
#seurat_annotations
#This is manual or reference-based cell type annotation. You (or a function like SingleR, Azimuth, or manual labeling) may have added it.

p1 + p2 #That will show whether your unsupervised clusters align with biological annotations.

#Step: Find all cluster markers

pbmc3k.markers <- FindAllMarkers(pbmc3k, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#only.pos = TRUE: only keep upregulated genes
#min.pct = 0.25: genes must be expressed in at least 25% of cells in either cluster
#logfc.threshold = 0.25: minimum log fold-change threshold

#view top markers per cluster:
pbmc3k.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

#Plot markers:

# Violin plot
VlnPlot(pbmc3k, features = c("MS4A1", "CD79A"))  # B-cell markers

# Feature plot (on UMAP)
FeaturePlot(pbmc3k, features = c("MS4A1", "CD3D", "GNLY"))

