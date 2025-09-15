library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)


# creating a Seurat Object out of csv table of raw counts for pancreas tissue
raw_counts <- read.csv(file = "/Users/svetlana/Desktop/SingleCell_test-task_for_lab/pancreas.raw.counts/pancreas_renamed_subtissue.csv",
                       sep = ",", header = TRUE, row.names = 1)
# making identity classes - by subtissues
colnames(raw_counts) <- gsub("\\.", '_', colnames(raw_counts))
sparse_counts <- as.sparse(as.matrix(raw_counts))
pancr <- CreateSeuratObject(counts = sparse_counts,
                            min.cells = 3, 
                            min.features = 200,
                            names.field = 8, # making identity classes - by subtissues
                            project = "pancreas")


## QC and cell selection

# removing samples with less than 200 genes and more that 7500:
pancr <- subset(pancr, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

scatter <- FeatureScatter(pancr, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


## Normalizing

pancr_norm <- NormalizeData(pancr, normalization.method = "LogNormalize", scale.factor = 10000)


pancr_norm_varfeat <- FindVariableFeatures(pancr_norm, 
                                           selection.method = "vst", 
                                           nfeatures = 3000)




## Linear transformation of the data:

all.genes <- rownames(pancr_norm_varfeat)
pancr_subtiss <- ScaleData(pancr_norm_varfeat, features = all.genes)



## Performing linear dimensional reduction

# using only variable features, determined previously (by default)
pancr_subtiss <- RunPCA(pancr_subtiss)
#print(pancr_linear[['pca']], nfeatures = 6)

#VizDimLoadings(pancr_linear, dims = 1:2, reduction = "pca")

PCA <- DimPlot(pancr_subtiss, reduction = "pca") + NoLegend()

#DimHeatmap(pancr_linear, dims = 1, cells = 500, balanced = TRUE)

#ElbowPlot(pancr_linear, ndims = 50)


## Clustering the cells

pancr_clust23_sub <- FindNeighbors(pancr_subtiss, dims = 1:23)
pancr_clust23_sub <- FindClusters(pancr_clust23_sub, resolution = 0.4)

#head(Idents(pancr_clust23_sub))



pancr_subtiss <- RunUMAP(pancr_subtiss, dims = 1:23)

umap_M_F_plot <- DimPlot(pancr_subtiss, 
                         reduction = "umap")
#

scatter 
PCA + umap_M_F_plot


ggplot(pancr_clust23_sub@meta.data, 
       aes(x = Idents(pancr_clust23_sub), fill = pancr_clust23_sub@meta.data$orig.ident)) +
  geom_bar(position = "fill") +
  labs(x = "Cluster", y = "Proportion", fill = "Type") +
  theme_minimal()



clusters.markers <- FindMarkers(pancr_subtiss, 
                                ident.1 = 'Endocrine', 
                                ident.2 = 'Exocrine')

clusters.markers %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 65) %>%
  ungroup() -> top65
DoHeatmap(pancr_subtiss, 
          features = row.names(top65)) + NoLegend()
