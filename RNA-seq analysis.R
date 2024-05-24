## Load packages as required ##

#install.packages("devtools")
devtools::install_github("JosephCrispell/basicPlotteR")
library(basicPlotteR)
library(devtools) 
library(pheatmap)
library(org.Mm.eg.db)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library (biomaRt)
library(EnsDb.Mmusculus.v79)

number <- 4 # determine the counts/strandedness
ff <- list.files( path = "/Counts/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 ) # need to skip the first 4 lines as we don't want the unmapped reads etc counting as genes + counts
counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/Counts//", "", ff )
ff <- paste0(sapply(strsplit(ff,"/Counts//"),"[", 2))
ff <-paste0(sapply(strsplit(ff, "EKDL"),"[",1))
ff <- gsub("_$", "", ff)
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

# for uniquely mapped reads for each file
ff <- list.files( path = "Counts/", pattern = "*ReadsPerGene.out.tab$", full.names = TRUE )
reads.files <- lapply( ff, read.table, nrows=1 )
reads <- as.data.frame(sapply (reads.files, function(x) x[ , number ]))
reads.rows <- gsub("_", " ", sampleID$x)
row.names(reads) <- reads.rows
colnames(reads) <-  "Uniquely mapped reads"
reads
write.table(reads, file="Counts/reads.txt", sep = "\t")

write.csv(ff, file="/Counts/counts_names.csv")
write.csv(counts, file="/Counts/counts.csv")

sampleID <- read.csv("/Counts/counts_meta_Condition.csv", header=TRUE) # create metadata file as appropriate
colnames(counts) <- sampleID$x
row.names(sampleID) <- sampleID$x 
all(rownames(sampleID) %in% colnames(counts))
all(rownames(sampleID) == colnames(counts))

dds_counts<-DESeqDataSetFromMatrix(countData = counts, colData = sampleID, design = ~ Group)

# merge duplicate read files
dds_counts$run<-paste0(dds_counts$x)
dds_counts$x <- paste(strsplit(dds_counts$x, "_2"))
dds_counts$x <- factor(sub("_$","",dds_counts$x))
ddsColl <- collapseReplicates(dds_counts, dds_counts$x, dds_counts$run)
colData(ddsColl)
colnames(ddsColl)
# View(ddsColl)
# sanity check merging
matchFirstLevel <- dds_counts$x == levels(dds_counts$x)[1]
stopifnot(all(rowSums(counts(dds_counts[,matchFirstLevel])) == counts(ddsColl[,1])))

ddsColl <- DESeq(ddsColl)

head(counts(ddsColl, normalized =TRUE))
countsnorm <- counts(ddsColl, normalized =TRUE)
write.table(countsnorm, file="/Counts/countsnorm.txt")

# regenerate sampleID with collapsed samples
sampleID.Collt <- data.frame(ddsColl@colData)
sampleID.Coll <- sampleID.Collt[,2:6]
write.table (sampleID.Coll, file='/Counts/sampleID.Coll.txt', sep = '\t')
res <- results(ddsColl)
res 
mcols(res, use.names=TRUE)

## PCA plots ##
vst.pca <- vst(ddsColl)
png(file='/Counts/PCA_all.png', width=1100, height=650)
v <- plotPCA(vst.pca, intgroup=c("Condition", "Timepoint"))
v <- v + scale_color_discrete(name="Condition:Timepoint") + geom_label_repel(aes(label=name), show.legend = FALSE)
v
dev.off()

png(file='/Counts/PCA_all_cellMW.png', width=4000, height=2000, res=300)
pcaData <- mwDESeq2PlotPCA(vst.pca, intgroup=c("Condition", "Timepoint","Cell"), returnData=TRUE, pc.a = 1, pc.b = 2) # function created by Dr Mark Walsh and called in
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$CellTimepoint <- paste(pcaData$Cell,pcaData$Timepoint)
ggplot(pcaData, aes(PC1, PC2, color=Cell, shape=Timepoint, fill=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_shape_manual(values = c(21, 24, 22)) + 
  scale_color_manual(values = c("#7197D9","black")) +
  scale_fill_manual(values = c("#D43F88","#FFDD4C", "#7662CA", "#BCF265")) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  guides(color = guide_legend(override.aes=list(shape=1))) + theme (axis.title.x = element_text(size = 14,
                                                                                                face = "bold"),
                                                                    axis.title.y = element_text(size = 14, 
                                                                                                face = "bold"),
                                                                    legend.title = element_text(size = 12),
                                                                    legend.text = element_text(size = 12))
dev.off()

png(file='/Counts/PCA_activin_MW.png', width=1100, height=650)
pcaData <- mwDESeq2PlotPCA(vst.pca, intgroup=c("Activin", "Timepoint"), returnData=TRUE, pc.a = 1, pc.b = 2)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Activin, shape=Timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()

pcaData <- plotPCA(vst.pca, intgroup=c("Condition", "Activin"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= Condition, shape=Activin)) +
  geom_point(size=3) + geom_label_repel(aes(label=name), show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
                               
## Pull DEGs between conditions ##
factor='Group'
reference = 'ACBI.D0.No.Act' ; # options --> Treatment.Timepoint.Act/No.Act
treatment =  'DMSO.D0.No.Act' ; # options --> Treatment.Timepoint.Act/No.Act
res = results(ddsColl, contrast = c(factor, treatment, reference)) 

## Volcano plots - toggle ref and treatment for different contrasts ##
log10.pval <- -log10(res$padj) 
log2.fc <- res$log2FoldChange

plot(log2.fc,log10.pval, xlab="log2 (fold change)", ylab="-log10 (p-value)", xlim=c(-10,10), ylim=c(0,75))

abline(h = -log10(0.05), col='red',lwd = 1.5)
abline(v = -log2(2), col='blue',lwd = 1.5)
abline(v = log2(2), col='blue',lwd = 1.5)

## XY plots ##

# re-generate res tables for XY plots
factor='Group'
reference = 'DMSO.D4.No.Act' ; # options --> DMSO, ACBI, PFI, 
treatment =  'ACBI.D4.No.Act' ; # options --> DMSO, ACBI, PFI/timepoint/Act or No.Act
res = results(ddsColl, contrast = c(factor, treatment, reference)) 

# XY plots continued with countsnorm
Pvalue = 0.05
topGenes<- res$padj < Pvalue
avgRef <- rowMeans( countsnorm[ , which(sampleID.Coll$Group==reference ) ] ) 
avgTreat <- rowMeans( countsnorm[ , which(sampleID.Coll$Group==treatment ) ] )
avglog2 <- log2( cbind(avgRef, avgTreat) + 0.5 ) 
plot(avglog2, xlab=reference, ylab=treatment)
points( avglog2[ which(topGenes) , ], col='red' )

# biomart to annotate with symbol and gene description

# listEnsembl()
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 109)
# mapping <-  getBM(attributes=c("mgi_symbol","ensembl_gene_id", "description"), filters = "ensembl_gene_id", mart=ensembl, values=ACBI_D$ensembl, uniqueRows=TRUE, bmHeader = T)
# map2 <- select (ensembl, keys = ACBI_D$ensembl, columns=c('ensembl_gene_id','mgi_symbol','description'), keytype = "ensembl_gene_id")
map2 <- select (ensembl, keys = ACBI_D$ensembl, keytype = "ensembl_gene_id", column="mgi_symbol")
require(biomaRt)
mart <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
## generate lookup table 
lookuptable <- getBM(
  mart = mart,
  attributes = c(
    'ensembl_gene_id',
    'mgi_symbol',
    'description'),
  uniqueRows = TRUE)
nrow(lookuptable)

write.table (lookuptable, file = "/Counts/lookuptableBM.txt", sep="\t")
                               
####### generate contrasts, add gene info and sort into significantly up-/downregulated genes #######
###### Day 0 DEGs ######
ACBI_D = results(ddsColl, contrast= c("Group", "ACBI.D0.No.Act", "DMSO.D0.No.Act")) 
ACBI_D$ensembl <- sapply( strsplit(rownames(ACBI_D), "\\."), "[", 1)
# ACBI_D$mgi <- mapIds(EnsDb.Mmusculus.v79, keys = ACBI_D$ensembl, keytype = "GENEID", column="SYMBOL")
# ACBI_D$description <- mapIds(EnsDb.Mmusculus.v79, keys = ACBI_D$ensembl, keytype = "ENSEMBL", column="GENENAME")
ACBI_D_all <- merge (as.data.frame(ACBI_D), lookuptable, by.x = "ensembl", by.y = "ensembl_gene_id")
ACBI_D <- merge (as.data.frame(ACBI_D), lookuptable, by.x = "ensembl", by.y = "ensembl_gene_id")
ACBI_D <- ACBI_D[which(ACBI_D$padj < 0.05),]
Ord.ACBI_D <- ACBI_D[order(ACBI_D$padj),] # 1
UOG.ACBI_D<- Ord.ACBI_D[which(Ord.ACBI_D$log2FoldChange>0), ] 
DOG.ACBI_D <- Ord.ACBI_D[which(Ord.ACBI_D$log2FoldChange<0), ] 
sDOG.ACBI_D<-DOG.ACBI_D[order(DOG.ACBI_D$padj),]
sUOG.ACBI_D<-UOG.ACBI_D[order(UOG.ACBI_D$padj), ] 
nrow(sDOG.ACBI_D)
nrow(sUOG.ACBI_D)
                               
ACBI_D_all <- na.omit(ACBI_D_all)
write.table(ACBI_D_all, file = "Counts/D0/D0.ACBI_D_all.txt", sep="\t")
write.table(Ord.ACBI_D, file = "Counts/D0/D0.Ord.ACBI_D.txt", sep="\t")
write.table(sDOG.ACBI_D, file = "Counts/D0/D0.sDOG.ACBI_D.txt", sep="\t")
write.table(sUOG.ACBI_D, file = "Counts/D0/D0.sUOG.ACBI_D.txt", sep="\t")

PFI_D = results(ddsColl, contrast= c("Group", "PFI.D0.No.Act", "DMSO.D0.No.Act")) 
PFI_D$ensembl <- sapply( strsplit(rownames(PFI_D), "\\."), "[", 1)
PFI_D$mgi <- mapIds(org.Mm.eg.db, keys = PFI_D$ensembl, keytype = "ENSEMBL", column="SYMBOL")
PFI_D$description <- mapIds(org.Mm.eg.db, keys = PFI_D$ensembl, keytype = "ENSEMBL", column="GENENAME")
PFI_D <- PFI_D[which(PFI_D$padj < 0.05),]
Ord.PFI_D <- PFI_D[order(PFI_D$padj),] # 0

CIS_D = results(ddsColl, contrast= c("Group", "CIS.D0.No.Act", "DMSO.D0.No.Act")) 
CIS_D$ensembl <- sapply( strsplit(rownames(CIS_D), "\\."), "[", 1)
CIS_D$mgi <- mapIds(org.Mm.eg.db, keys = CIS_D$ensembl, keytype = "ENSEMBL", column="SYMBOL")
CIS_D$description <- mapIds(org.Mm.eg.db, keys = CIS_D$ensembl, keytype = "ENSEMBL", column="GENENAME")
CIS_D <- CIS_D[which(CIS_D$padj < 0.05),]
Ord.CIS_D <- CIS_D[order(CIS_D$padj),] # 0

ACBI_C = results(ddsColl, contrast= c("Group", "ACBI.D0.No.Act", "CIS.D0.No.Act")) 
ACBI_C$ensembl <- sapply( strsplit(rownames(ACBI_C), "\\."), "[", 1)
ACBI_C_all <- merge (as.data.frame(ACBI_C), lookuptable, by.x = "ensembl", by.y = "ensembl_gene_id")
ACBI_C<- ACBI_C[which(ACBI_C$padj < 0.05),]
Ord.ACBI_C <- ACBI_C[order(ACBI_C$padj),]
nrow(Ord.ACBI_C) #0

####### generate as many contrasts as you like #######
                               
###### make a seperate heatmap for each comparison ######

# vst is a quicker transformation! #
rld<-rlog(ddsColl)
hm<-(assay(rld))

### Heatmaps of DEGs ###

# ACBI_D
d0.ACBI_D = results(ddsColl, contrast= c("Group", "ACBI.D0.No.Act", "DMSO.D0.No.Act")) 
d0.ACBI_D <- d0.ACBI_D[which(d0.ACBI_D$padj < 0.05),]
hm.samples<-hm[ , c(d0.ACBI_D_info[,1])]
hm.samples<-hm.samples[(rownames(d0.ACBI_D)), ]
pheatmap(hm.samples, fontsize_row=1, scale="row", treeheight_row = 0, main= "Day 0 ACBI1 vs DMSO, no ActA",legend=TRUE)

d2.ACBI_D = results(ddsColl, contrast= c("Group", "ACBI.D2.No.Act", "DMSO.D2.No.Act")) 
d2.ACBI_D <- d2.ACBI_D[which(d2.ACBI_D$padj < 0.05),]
d2.ACBI_D_info <- read.table(file = "d2.ACBI_D_info.txt",header=TRUE) # you can limit info to relevant samples only
hm.samples<-hm[ , c(d2.ACBI_D_info[,1])]
hm.samples<-hm.samples[(rownames(d2.ACBI_D)), ]
anno_col <- d2.ACBI_D_info["Condition"]
rownames(anno_col) <- colnames(hm.samples)
labels <- sapply(strsplit(c(d2.ACBI_D_info[,1]), "_D2"), "[", 1)
labels <- gsub("_", " ", labels)
png(file='hm_d2_ACBI_D.png', width=2000, height=2500, res=300)
pheatmap(hm.samples, fontsize_row=3, scale="row", treeheight_row = 0, main= "Day 2 ACBI1 vs DMSO, no ActA",legend=TRUE, 
         annotation_col = anno_col,labels_col = labels, fontsize_col = 15)
dev.off()

# includes PFI3 and CIS
d2.ACBI_D_all_info <- read.table(file = "d2.ACBI_D_all_info.txt",header=TRUE)
hm.samples<-hm[ , c(d2.ACBI_D_all_info[,1])]
hm.samples<-hm.samples[(rownames(d2.ACBI_D)), ]
anno_col <- d2.ACBI_D_all_info["Condition"]
rownames(anno_col) <- colnames(hm.samples)
png(file='hm_d2_ACBI_D_all.png', width=2000, height=2500, res=300)
pheatmap(hm.samples, fontsize_row=3, scale="row", treeheight_row = 0, main= "Day 2 treatment",legend=TRUE, annotation_col = anno_col
         ,labels_col = sapply(strsplit(c(d2.ACBI_D_all_info[,1]), "_D2"), "[", 1)) # note different labelling strategies here
dev.off()

d4.ACBI_D = results(ddsColl, contrast= c("Group", "ACBI.D4.No.Act", "DMSO.D4.No.Act")) 
d4.ACBI_D <- d4.ACBI_D[which(d4.ACBI_D$padj < 0.05),]
d4.ACBI_D_info <- read.table(file = "d4.ACBI_D_info.txt",header=TRUE)
hm.samples<-hm[ , c(d4.ACBI_D_info[,1])]
hm.samples<-hm.samples[(rownames(d4.ACBI_D)), ]
anno_col <- d4.ACBI_D_info["Condition"]
rownames(anno_col) <- colnames(hm.samples)
labels <- sapply(strsplit(c(d4.ACBI_D_info[,1]), "_D4"), "[", 1)
labels <- gsub("_", " ", labels)
png(file='hm_d4_ACBI_D.png', width=2000, height=2500, res = 300)
pheatmap(hm.samples, fontsize_row=1, scale="row", treeheight_row = 0, main= "Day 4 ACBI1 vs DMSO, no ActA", legend=TRUE, annotation_col = anno_col,
         labels_col = labels, fontsize_col = 15)
dev.off()

# includes PFI3 and CIS
d4.ACBI_D = results(ddsColl, contrast= c("Group", "ACBI.D4.No.Act", "DMSO.D4.No.Act")) 
d4.ACBI_D <- d4.ACBI_D[which(d4.ACBI_D$padj < 0.05),]
d4.ACBI_D_all_info <- read.table(file = "d4.ACBI_D_all_info.txt",header=TRUE)
hm.samples<-hm[ , c(d4.ACBI_D_all_info[,1])]
hm.samples<-hm.samples[(rownames(d4.ACBI_D)), ]
anno_col <- d4.ACBI_D_all_info["Condition"]
rownames(anno_col) <- colnames(hm.samples)
png(file='=hm_d4_ACBI_D_all.png', width=535, height=625)
pheatmap(hm.samples, fontsize_row=1, scale="row", treeheight_row = 0, main= "Day 4 treatment", legend=TRUE, annotation_col = anno_col,
          labels_col = sapply(strsplit(c(d4.ACBI_D_all_info[,1]), "_D4"), "[", 1))
dev.off()

## All samples/conditions - a little messy!

allres = results(ddsColl) 
allres <- allres[which(allres$padj < 0.05),]
allres_info <- sampleID.Coll
hm.samples<-hm.samples[(rownames(allres)), ]
png(file='hm_allres.png', width=535, height=625)
pheatmap(hm.samples, fontsize_row=2, scale="row", treeheight_row = 0, main= "All conditions", legend=TRUE, annotation_col = sampleID.Coll["Group"])
dev.off()
