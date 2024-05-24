# Build Seurat object from matrix and metadata csv
# Data obtained from Jiang et al, 2023: https://doi.org/10.1038/s41421-022-00513-z

library(Seurat);library("RColorBrewer")

setwd("scRNAseq_Jiang2023/scRNA-seq compressed csv")

data <- read.csv (file ="MTR15682-matrix.csv", header = TRUE)
meta <- read.csv (file = "MTR15682-meta.csv", header = TRUE, row.names = 1)
umap_coordinates <- read.csv (file = "UMAP-MTR15682.csv", header = TRUE, row.names = 1)

View(data)
View(meta)
View(umap_coordinates)

# change the final .number to match -number of metadata 
cellbarcodes <- colnames(data)
cellbarcodes <- gsub('(.+)[.](\\d+)$', '\\1-\\2', cellbarcodes)
cellbarcodes
# add genes to rownames and remove old column
colnames(data) <- cellbarcodes
rownames(data) <- data[, 1]
data <- data[, -1]

# Create Seurat Object
obj <- CreateSeuratObject(counts = data, meta.data = meta)
View(obj)

# Add umap reduction
umap_coordinates_mat <- as(umap_coordinates, "matrix")
obj[['UMAP']] <- CreateDimReducObject(embeddings = umap_coordinates_mat, key = "UMAP_", global = T, assay = "RNA")
obj <- NormalizeData(obj)

DimPlot(obj, group.by = "orig.ident")
DimPlot(obj, group.by = "seurat_clusters")
DimPlot(obj, group.by = "type_19")
BAF_feat<- c ("Smarca4", "Smarca2", "Pbrm1")
png(file='JiangFeatPlot_scaled_res300.png', width=2250, height=2000, res=300)
FeaturePlot(obj, features = BAF_feat)
dev.off()

png(file='JiangDotPlot_scaled_res300.png', width=2000, height=2000, res=300)
p <- Seurat::DotPlot(obj, 
                features = BAF_feat, 
                cols = c("magenta","blue"),
                col.min = 0,
                col.max = 4,
                dot.min = 0,
                dot.scale = 6,
                idents = NULL,
                group.by = "type_19",
                split.by = NULL,
                cluster.idents = FALSE,
                scale = FALSE,
                scale.by = "radius",
                scale.min = NA,
                 scale.max = NA) 
p + labs(y = "Trophoblast populations", x = "BAF complex components") + theme (axis.title.x = element_text(size = 14,
                                                                                                                face = "bold"),
                                                                                    axis.title.y = element_text(size = 14, 
                                                                                                                face = "bold"),
                                                                                    legend.title = element_text(size = 12),
                                                                                    legend.text = element_text(size = 12))
dev.off()
