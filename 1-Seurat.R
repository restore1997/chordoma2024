rm(list = ls())

library(dplyr)
library(Seurat)

load("../Raw_total.rda")
table(total$orig.ident)

total[["percent.mt"]] <- PercentageFeatureSet(total, pattern = "^MT-")
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0,group.by = "sample")

total <- subset(total, subset = nFeature_RNA > 500 & percent.mt < 25)
total
table(total$orig.ident)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,pt.size = 0)


total <- NormalizeData(total, normalization.method = "LogNormalize", scale.factor = 10000)

total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(total), 10)

plot1 <- VariableFeaturePlot(total)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(total)
total <- ScaleData(total, features = all.genes)

total <- RunPCA(total, features = VariableFeatures(object = total))
print(total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(total, dims = 1:2, reduction = "pca")
DimPlot(total, reduction = "pca")

total <- FindNeighbors(total, dims = 1:10)
total <- FindClusters(total, resolution = 0.5)

table(total$seurat_clusters)

total <- RunUMAP(total, dims = 1:10)

DimPlot(total, reduction = "umap",label = T,raster=FALSE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
total.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- total.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(total, features = top10$gene) + NoLegend()

save(total,total.markers,file = "total.rda")