library (Seurat)
library (SeuratObject)
library (ggrepel)
library (patchwork)
library(ggplot2)

load ("/scRNA analysis published/AllStages_TrophoblastNuclei_obj.Rdata")
mouse.troph.combined.updated = UpdateSeuratObject(object = mouse.troph.combined)

png(file='Marsh_nolabel_res300.png', width=2000, height=2000, res=300)
DimPlot(mouse.troph.combined.updated, reduction = "umap", label = FALSE, repel = TRUE, label.size = 5) + NoLegend() 
dev.off()

DimPlot(mouse.troph.combined, reduction = "umap", label = TRUE, split.by = "GA") + NoLegend()

mouse.troph.combined.updated@meta.data
cluster.averages <- AverageExpression(mouse.troph.combined.updated)

features <- c("Ctnnbip1", "Creg1", "Ada", "1500009L16Rik", "Sema3f", "Hic1", "Hcrtr1", "Prl7a1", "Prl7d1", "Prl2c5", "Lgals9", "Rell2")
BAF_feat<- c ("Smarca2", "Smarca4", "Pbrm1")
stem_feat<- c ("Cdx2", "Eomes", "Elf5", "Esrrb")
Lab_PC<- c ("Ovol2", "Gcm1", "Dlx3", "Cebpa")
Lab_feat <- c ("Syna", "Synb", "Ctsq", "Prl3b1")
SpT_PC<- c ("Prdm1")
JZ_feat <- c ("Gjb3", "Tpbpa", "Ascl2", "Cdkn1c")
Prolactins_feat <- c ("Prl3d1", "Prl3b1", "Prl2c2", "Prl7b1")
feat <- c("Foxj3")

#
RidgePlot(mouse.troph.combined.updated, features =BAF_feat, ncol = 2)
VlnPlot(mouse.troph.combined.updated, features = BAF_feat, pt.size = 0.1)
FeaturePlot(mouse.troph.combined.updated, features = BAF_feat, keep.scale = "feature")
FeaturePlot(mouse.troph.combined.updated, features = "Eomes")
FeaturePlot(mouse.troph.combined.updated, features = "Foxj3")

#BAF/BPTF - all subunits
png(file="/scRNA analysis published/featureplot_BAF.png", width=484, height=484)
par(mfrow=c(1,1))
BAF_feat<- c ("Smarca2", "Smarca4", "Pbrm1", "Arid1a", "Arid1b", "Arid2", "Ss18", "Ss18l1","Smarcc1", "Smarcc2","Smarcd1", 
              "Smarcd2", "Smarcd3", "Smarcb1", "Brd7", "Brd9", "Bcl7a", "Bcl7b", "Bcl7c", "Bcl11a", "Bcl11b",
              "Phf10", "Smarce1", "Actl6a", "Actl6b", "Dpf1", "Dpf2", "Dpf3" )
FeaturePlot(mouse.troph.combined.updated, features = BAF_feat) 
dev.off()
FeaturePlot(mouse.troph.combined, features = BAF_feat, split.by = "GA")

png(file='MarshDotPlot_scaled_res300.png', width=2000, height=2000, res=300)
p <- Seurat::DotPlot(mouse.troph.combined.updated, 
        features = BAF_feat, 
        cols = c("magenta","blue"),
        col.min = 0,
        col.max = 4,
        dot.min = 0,
        dot.scale = 6,
        idents = NULL,
        group.by = NULL,
        split.by = NULL,
        cluster.idents = FALSE,
        scale = FALSE,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) 
p <- p + labs(y = "Trophoblast populations", x = "BAF complex components") + theme (axis.title.x = element_text(size = 14,
                                                                                                      face = "bold"),
                                                                          axis.title.y = element_text(size = 14, 
                                                                                                      face = "bold"),
                                                                          legend.title = element_text(size = 12),
                                                                          legend.text = element_text(size = 12))
dev.off()
p


# reorder clusters
levels(mouse.troph.combined.updated@active.ident)

mouse.troph.combined.updated@active.ident <- factor( mouse.troph.combined.updated@active.ident,
                                                     levels = c("LaTP",
                                                                "LaTP 2",
                                                                "SynTI Precursor",
                                                                "SynTI",
                                                                "SynTII Precursor",
                                                                "SynTII",
                                                                "S-TGC Precursor",
                                                                "S-TGC",
                                                                "JZP 1",
                                                                "JZP 2",
                                                                "SpT Precursor",
                                                                "SpT",
                                                                "Glycogen Cells"))

png(file='MarshVlnPlot_scaled_res300.png', width=4000, height=2000, res=300)
VlnPlot(mouse.troph.combined.updated, features = BAF_feat, pt.size = 0.1)
dev.off()

#Stemness markers
png(file="featureplot_stem.png", width=484, height=484)
FeaturePlot(mouse.troph.combined, features = stem_feat)
dev.off()
FeaturePlot(mouse.troph.combined.updated, features = BAF_feat, split.by ="GA")

# Finding differentially expressed markers 
all_markers<-FindAllMarkers(mouse.troph.combined)
head(all_markers)

write.csv(all_markers, file = "GSEA/all_markers_Marsh_scdata.csv")
