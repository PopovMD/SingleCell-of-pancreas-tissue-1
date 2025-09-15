library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)


# creating a Seurat Object out of csv table of raw counts for pancreas tissue
raw_counts <- read.csv(file = "./pancreas.raw.counts/Pancreas-counts.csv",
                       sep = ",", header = TRUE, row.names = 1)
colnames(raw_counts) <- gsub("\\.", '_', colnames(raw_counts))
sparse_counts <- as.sparse(as.matrix(raw_counts))
pancr <- CreateSeuratObject(counts = sparse_counts,
                            min.cells = 3, 
                            min.features = 200,
                            names.field = 5,
                            project = "pancreas")
head(pancr@meta.data, 10)
pancr@meta.data$orig.ident
Idents(pancr)

## QC and cell selection

# adding a new column to the meta.data, that contains the persantage of reads 
# that map to the mitochondrial genome (by using [[ ]]).
pancr[["percent.mt"]] <- PercentageFeatureSet(pancr, pattern = "^MT-")
head(pancr@meta.data, 5)

plot1 <- FeatureScatter(pancr, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pancr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2

sum(pancr@meta.data$percent.mt != 0) # = 0
# митохонриальных генов в принципе нет, что настораживает.








############
# после того, как стало видно, что много образцов, у которых большое количество генов, решил проверить другую ткань:
############
# checking another tissue
raw_counts_fat <- read.csv(file = "/Users/svetlana/Desktop/SingleCell_test-task_for_lab/pancreas.raw.counts/Pancreas-counts.csv",
                       sep = ",", header = TRUE, row.names = 1)
sparse_counts_fat <- as.sparse(as.matrix(raw_counts_fat))
fat <- CreateSeuratObject(counts = sparse_counts_fat, project = "pancreas", min.cells = 3, min.features = 200)
head(fat@meta.data, 5)

fat[["percent.mt"]] <- PercentageFeatureSet(fat, pattern = "^MT-")
head(fat@meta.data, 5)

plot1 <- FeatureScatter(fat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(fat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

sum(fat@meta.data$percent.mt != 0) # = 0

# and another, just to ensure myself
raw_counts_mg <- read.csv(file = "/Users/svetlana/Desktop/SingleCell_test-task_for_lab/microglia.raw.counts/Brain_Microglia-counts.csv",
                           sep = ",", header = TRUE, row.names = 1)
sparse_counts_mg <- as.sparse(as.matrix(raw_counts_mg))
mg <- CreateSeuratObject(counts = sparse_counts_mg, project = "pancreas", min.cells = 3, min.features = 200)
head(mg@meta.data, 5)

mg[["percent.mt"]] <- PercentageFeatureSet(mg, pattern = "^MT-")
head(mg@meta.data, 5)

plot1 <- FeatureScatter(mg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

sum(mg@meta.data$percent.mt != 0) # = 0
##########
# everything is the same, except for the microglia's amounts of genes:
# for this tissue where are less number of cells, containing more than 4000 of them
# So I would assume, that these counts were already cleared from
# mitochondrial genes or maybe only nucleus of cells were sequenced.

# so i will continue my work with pancreas
##########






# removing samples with less than 200 genes and more that 7500:
pancr <- subset(pancr, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
# maybe it makes sense to remove samples with too many reads

FeatureScatter(pancr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


## Normalizing

pancr_norm <- NormalizeData(pancr, normalization.method = "LogNormalize", scale.factor = 10000)
# LogNormalize: normalizing by total expression, multiplying by scale.factor and then log-transforms the result
# assumption: each cell originally contains the same number of RNA molecules. 



## Selection of highly variable features

sum(pancr_norm@meta.data$nFeature_RNA) # =  6759903
pancr_norm_varfeat <- FindVariableFeatures(pancr_norm, 
                                           selection.method = "vst", 
                                           nfeatures = 3000)
plot1 <- VariableFeaturePlot(pancr_norm_varfeat)
plot1

top10 <- head(VariableFeatures(pancr_norm_varfeat), 5)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

nrow(pancr_norm_varfeat@meta.data)



## Linear transformation of the data:

all.genes <- rownames(pancr_norm_varfeat)
# transforming all of the genes
pancr_linear <- ScaleData(pancr_norm_varfeat, features = all.genes)



## Performing linear dimensional reduction

# using only variable features, determined previously (by default)
pancr_linear <- RunPCA(pancr_linear)
print(pancr_linear[['pca']], nfeatures = 6)

# For the first principal components, Seurat outputs a list of genes 
# with the most positive and negative loadings, representing 
# modules of genes that exhibit either correlation (or anti-correlation) 
# across single-cells in the dataset.

VizDimLoadings(pancr_linear, dims = 1:2, reduction = "pca")

DimPlot(pancr_linear, reduction = "pca") + NoLegend()

DimHeatmap(pancr_linear, dims = 1:6, cells = 500, balanced = TRUE)
#DimHeatmap(pancr_linear, dims = 1:30, cells = 500, balanced = TRUE)

# Choosing the number of principal components by percentage of variance explained by each one
ElbowPlot(pancr_linear, ndims = 50)



## Clustering the cells

pancr_clust23 <- FindNeighbors(pancr_linear, dims = 1:23)
pancr_clust23 <- FindClusters(pancr_clust23, resolution = 0.4)

head(Idents(pancr_clust23))



### Non-linear dimensional reduction

pancr_linear <- RunUMAP(pancr_linear, dims = 1:23)

umap_M_F_plot <- DimPlot(pancr_linear, 
        reduction = "umap")

pancr_clust23 <- RunUMAP(pancr_clust23, dims = 1:23)

umap_clusters_plot <- DimPlot(pancr_clust23, 
                              reduction = "umap",
                              label = TRUE)

umap_clusters_plot + umap_M_F_plot

head(pancr_linear)
head(pancr_clust23)


library(ggplot2)

# Bar plot
ggplot(pancr_clust23@meta.data, 
       aes(x = Idents(pancr_clust23), fill = pancr_clust23@meta.data$orig.ident)) +
       geom_bar(position = "fill") +
       labs(x = "Cluster", y = "Proportion", fill = "Gender") +
  theme_minimal()











