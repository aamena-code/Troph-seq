library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(UpSetR)
library(DiffBind)
library(ggplot2)
library(gridExtra)
library(egg)

db <- readRDS ("DiffBind_analysis.rds")

View(db)
plot(db)
DiffBind::dba.plotPCA(db, attributes = DBA_REPLICATE)

DiffBind::dba.plotHeatmap(db, correlations = FALSE)
DiffBind::dba.plotBox(db, contrast = 2)
DiffBind::dba.plotVenn(db, contrast = 9, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
DiffBind::dba.plotVenn(db, contrast = 30, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)

profiles <- DiffBind::dba.plotProfile(db, merge = c(DBA_CONDITION, DBA_REPLICATE)) 
DiffBind::dba.plotProfile(profiles)

repObj_A <- dba.report(db, contrast=9, bDB=TRUE, bGain=TRUE, bLoss=TRUE,bAll=FALSE)
profiles0 <- DiffBind::dba.plotProfile(db, sites=repObj_A,
                                       samples=list(ACBI1=db$masks$D0_ACBI1,
                                                    DMSO=db$masks$D0_DMSO))


repObj_A.d2 <- dba.report(db, contrast=30, bDB=TRUE, bGain=TRUE, bLoss=TRUE,bAll=FALSE)
profiles2 <- DiffBind::dba.plotProfile(db, sites=repObj_A.d2,
                                       samples=list(ACBI1=db$masks$D2_ACBI1,
                                                    DMSO=db$masks$D2_DMSO))

png(file='DB_Profiles_D0.png', width=1800, height=2200, res=300)
DiffBind::dba.plotProfile(profiles0)
dev.off()

png(file='DB_Profiles_D2.png', width=1800, height=2200, res=300)
DiffBind::dba.plotProfile(profiles2)
dev.off()

png(file='PCA_all.png', width=2000, height=1500, res=300)
DiffBind::dba.plotPCA(db, attributes = DBA_CONDITION)
dev.off()

png(file='heatmap_all.png', width=1000, height=750)
plot(db)
dev.off()

hmap <- colorRampPalette(c("blue", "grey", "magenta"))(n = 13)
readscores <- dba.plotHeatmap(db, contrast=9, correlations=FALSE, scale="row", colScheme = hmap)
readscores <- dba.plotHeatmap(db, contrast=30, correlations=FALSE, scale="row", colScheme = hmap)

############################## PEAK ANNOTATION ##############################


atac.con <- list.files( path = "bed_files", pattern = "*.bed.peak.txt", full.names = TRUE ) # replace with relevant file suffix
#atac.files <- lapply (atac.con, read.table)


txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
seqlevelsStyle(txdb) <- "UCSC"
genome(annotations) <- "mm39" # mm39

# Get the current sequence levels
seq_levels <- seqlevels(txdb)
seqlevels (peaks.down.d0.ACBI_D)

# Define a pattern to match the "_random" suffix
pattern <- "_random$"

# Replace "_random" suffix with an empty string
seq_levels <- gsub(pattern, "", seq_levels)

# Update the sequence levels in the TxDb object
seqlevels(txdb) <- seq_levels

# Function to correct chromosome names
correct_chromosome_names <- function(chromosomes) {
  pattern <- "^([^chr]+)\\.(\\d+)$"  # Pattern for matching chromosomes with format like "GL456368.1"
  corrected_names <- chromosomes
  mismatched_indices <- grep(pattern, chromosomes)
  for (i in mismatched_indices) {
    corrected_names[i] <- gsub(pattern, "chrUn_\\1v\\2", corrected_names[i])
  }
  return(corrected_names)
}

chromosomes_data <- unique(seqlevels(peaks.up.d0.ACBI_D))

# Compare with chromosome names from TxDb
chromosomes_txdb <- seqlevels(txdb)
mismatches <- setdiff(chromosomes_txdb, chromosomes_data)

# Print mismatches
print(mismatches)

# downDARs = less accessible regions, upDARs = more accessible regions
##### Day 0 ACBI vs DMSO #####

# downDARs - contrast 9
pf.down.d0.ACBI_D <- atac.con[60]
peaks.down.d0.ACBI_D <- readPeakFile(pf.down.d0.ACBI_D)
seqlevels(peaks.down.d0.ACBI_D) <- correct_chromosome_names(seqlevels(peaks.down.d0.ACBI_D)) 

#check whether there are still mismatches (see above)
seqlevels(peaks.down.d0.ACBI_D) <- sub("chrUn_GL456233v2", "chrX_GL456233v2", seqlevels(peaks.down.d0.ACBI_D))
seqlevels(peaks.down.d0.ACBI_D) <- sub("chrUn_GL456211v1", "chr1_GL456211v1", seqlevels(peaks.down.d0.ACBI_D))

write.table(peaks.down.d0.ACBI_D, file = "peaks.down.d0.ACBI_D.txt", sep="\t", row.names = FALSE, col.names = FALSE)

peakAnno.down.d0.ACBI_D <- annotatePeak(peaks.down.d0.ACBI_D, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d0.ACBI_D)
downDARs.d0.ACBI_D<-as.GRanges(peakAnno.down.d0.ACBI_D)
downDARs.d0.ACBI_D$geneId
downDARs.d0.ACBI_D$SYMBOL

write.table(downDARs.d0.ACBI_D, file = "downDARs.d0.ACBI_D_allfeat.txt", sep="\t")

write.table(downDARs.d0.ACBI_D$SYMBOL, file = "downDARs.d0.ACBI_D.txt", sep="\t")

plotAnnoBar(peakAnno.down.d0.ACBI_D)
plotAnnoPie(peakAnno.down.d0.ACBI_D)
vennpie(peakAnno.down.d0.ACBI_D)
upsetplot(peakAnno.down.d0.ACBI_D) 

png(file='Upset_D0.Down.ACBI_D.png', width=1500, height=2000, res=300)
upsetplot(peakAnno.down.d0.ACBI_D) 
dev.off()

png(file='VennPie_D0.Down.ACBI_D.png', width=4000, height=2000, res=300)
vennpie(peakAnno.down.d0.ACBI_D, cex = 1.25)
dev.off()

plotDistToTSS(peakAnno.down.d0.ACBI_D,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

peakAnnoList <- c(peakAnno.up.d0.ACBI_D,peakAnno.down.d0.ACBI_D)
plotAnnoBar(peakAnnoList)

# upDARs - contrast 9

pf.up.d0.ACBI_D <- atac.con[120]
peaks.up.d0.ACBI_D <- readPeakFile(pf.up.d0.ACBI_D)
seqlevels(peaks.up.d0.ACBI_D) <- correct_chromosome_names(seqlevels(peaks.up.d0.ACBI_D))

peakAnno.up.d0.ACBI_D <- annotatePeak(peaks.up.d0.ACBI_D, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d0.ACBI_D)
upDARs.d0.ACBI_D<-as.GRanges(peakAnno.up.d0.ACBI_D)
upDARs.d0.ACBI_D$geneId
upDARs.d0.ACBI_D$SYMBOL

write.table(upDARs.d0.ACBI_D$geneId, file = "upDARs.d0.ACBI_D.txt", sep="\t")

write.table(upDARs.d0.ACBI_D, file = "upDARs.d0.ACBI_D_allfeat.txt", sep="\t")

plotAnnoBar(peakAnno.up.d0.ACBI_D)
plotAnnoPie(peakAnno.up.d0.ACBI_D)
vennpie(peakAnno.up.d0.ACBI_D)
upsetplot(peakAnno.up.d0.ACBI_D)

png(file='Upset_D0.Up.ACBI_D.png', width=1500, height=2000, res=300)
upsetplot(peakAnno.up.d0.ACBI_D) 
dev.off()

plotDistToTSS(peakAnno.up.d0.ACBI_D,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

##### Day 2 ACBI vs DMSO #####

# down DARs - contrast 30
pf.down.d2.ACBI_D <- atac.con[24]
peaks.down.d2.ACBI_D <- readPeakFile(pf.down.d2.ACBI_D)
seqlevels(peaks.down.d2.ACBI_D) <- correct_chromosome_names(seqlevels(peaks.down.d2.ACBI_D))
seqlevels(peaks.down.d2.ACBI_D) <- sub("chrUn_GL456233v2", "chrX_GL456233v2", seqlevels(peaks.down.d2.ACBI_D))
seqlevels(peaks.down.d2.ACBI_D) <- sub("chrUn_MU069434v1", "chr1_MU069434v1", seqlevels(peaks.down.d2.ACBI_D))

peakAnno.down.d2.ACBI_D <- annotatePeak(peaks.down.d2.ACBI_D, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d2.ACBI_D)
downDARs.d2.ACBI_D<-as.GRanges(peakAnno.down.d2.ACBI_D)
downDARs.d2.ACBI_D$geneId
downDARs.d2.ACBI_D$SYMBOL

write.table(downDARs.d2.ACBI_D, file = "downDARs.d2.ACBI_D_allfeat.txt", sep="\t")
write.table(downDARs.d2.ACBI_D$SYMBOL, file = "downDARs.d2.ACBI_D.txt", sep="\t")

png(file='Upset_D2.Down.ACBI_D.png', width=1500, height=2000, res=300)
upsetplot(peakAnno.down.d2.ACBI_D)
dev.off()

# upDARs - contrast 30

pf.up.d2.ACBI_D <- atac.con[84]
peaks.up.d2.ACBI_D <- readPeakFile(pf.up.d2.ACBI_D)
seqlevels(peaks.up.d2.ACBI_D) <- correct_chromosome_names(seqlevels(peaks.up.d2.ACBI_D))

peakAnno.up.d2.ACBI_D <- annotatePeak(peaks.up.d2.ACBI_D, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d2.ACBI_D)
upDARs.d2.ACBI_D<-as.GRanges(peakAnno.up.d2.ACBI_D)
upDARs.d2.ACBI_D$geneId

write.table(upDARs.d2.ACBI_D$geneId, file = "upDARs.d2.ACBI_D.txt", sep="\t")
write.table(upDARs.d2.ACBI_D, file = "upDARs.d2.ACBI_D_allfeat.txt", sep="\t")

plotAnnoBar(peakAnno.up.d2.ACBI_D)
vennpie(peakAnno.up.d2.ACBI_D)

png(file='Upset_D2.Up.ACBI_D.png', width=1500, height=2000, res=300)
upsetplot(peakAnno.up.d2.ACBI_D)
dev.off()

plotDistToTSS(peakAnno.up.d2.ACBI_D,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

# arrange dist to TSS plots in one panel
TSS.d0.down <- plotDistToTSS(peakAnno.down.d0.ACBI_D,
                             title="Less accessible regions in ACBI1 treated TSCs")
TSS.d0.up <- plotDistToTSS(peakAnno.up.d0.ACBI_D,
                           title="More accessible regions in ACBI1 treated TSCs")
TSS.d2.down <- plotDistToTSS(peakAnno.down.d2.ACBI_D,
                             title="Less accessible regions in ACBI1 treated D2 diff")
TSS.d2.up <- plotDistToTSS(peakAnno.up.d2.ACBI_D,
                           title="More accessible regions in ACBI1 treated D2 diff")

print(TSS.d0.down)
print(TSS.d0.up)
print(TSS.d2.down)
print(TSS.d2.up)
title_gpar <- gpar(fontsize = 16)

png(file='DIST_TO_TSS.png', width=3000, height=2000, res=300)
grid.arrange(TSS.d0.down, TSS.d0.up, TSS.d2.down, TSS.d2.up,
             ncol = 2, top = textGrob("Distribution of DARs relative to TSS", gp = title_gpar))
dev.off()

# alternative arrangements
peakAnnoList <- c(peakAnno.down.d0.ACBI_D, peakAnno.up.d0.ACBI_D, peakAnno.down.d2.ACBI_D, peakAnno.up.d2.ACBI_D)

plotAnnoBar(peakAnnoList[[1]]) # check which is which
plotAnnoBar(peakAnnoList[[2]])
plotAnnoBar(peakAnnoList[[3]])
plotAnnoBar(peakAnnoList[[4]])
plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList)


##### Day 2 ACBI vs DMSO #####

# down DARs - contrast 30
pf.down.d2.ACBI_D <- atac.con[24]
peaks.down.d2.ACBI_D <- readPeakFile(pf.down.d2.ACBI_D)

peakAnno.down.d2.ACBI_D <- annotatePeak(peaks.down.d2.ACBI_D, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d2.ACBI_D)
downDARs.d2.ACBI_D<-as.GRanges(peakAnno.down.d2.ACBI_D)
downDARs.d2.ACBI_D$geneId

write.table(downDARs.d2.ACBI_D$geneId, file = "downDARs.d2.ACBI_D.txt", sep="\t")

# upDARs - contrast 30

pf.up.d2.ACBI_D <- atac.con[84]
peaks.up.d2.ACBI_D <- readPeakFile(pf.up.d2.ACBI_D)

peakAnno.up.d2.ACBI_D <- annotatePeak(peaks.up.d2.ACBI_D, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d2.ACBI_D)
upDARs.d2.ACBI_D<-as.GRanges(peakAnno.up.d2.ACBI_D)
upDARs.d2.ACBI_D$geneId

write.table(upDARs.d2.ACBI_D$geneId, file = "upDARs.d2.ACBI_D.txt", sep="\t")

plotAnnoBar(peakAnno.up.d2.ACBI_D)
vennpie(peakAnno.up.d2.ACBI_D)
upsetplot(peakAnno.up.d2.ACBI_D)


##### Day 2 ACBI vs cis-ACBI #####

# down DARs - contrast 24
pf.down.d2.ACBI_C <- atac.con[17]
peaks.down.d2.ACBI_C <- readPeakFile(pf.down.d2.ACBI_C)

peakAnno.down.d2.ACBI_C <- annotatePeak(peaks.down.d2.ACBI_C, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d2.ACBI_C)
downDARs.d2.ACBI_C<-as.GRanges(peakAnno.down.d2.ACBI_C)
downDARs.d2.ACBI_C$geneId

write.table(downDARs.d2.ACBI_C$geneId, file = "downDARs.d2.ACBI_C.txt", sep="\t")

# upDARs - contrast 24

pf.up.d2.ACBI_C <- atac.con[77]
peaks.up.d2.ACBI_C <- readPeakFile(pf.up.d2.ACBI_C)

peakAnno.up.d2.ACBI_C <- annotatePeak(peaks.up.d2.ACBI_C, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d2.ACBI_C)
upDARs.d2.ACBI_C<-as.GRanges(peakAnno.up.d2.ACBI_C)
upDARs.d2.ACBI_C$geneId

write.table(upDARs.d2.ACBI_C$geneId, file = "upDARs.d2.ACBI_C.txt", sep="\t")

plotAnnoBar(peakAnno.up.d2.ACBI_C)
vennpie(peakAnno.up.d2.ACBI_C)
upsetplot(peakAnno.up.d2.ACBI_C)

##### Day 0 ACBI vs Day 2 ACBI #####

# down DARs - contrast 2

pf.down.d0.A_d2.A <- atac.con[12]
peaks.down.d0.A_d2.A <- readPeakFile(pf.down.d0.A_d2.A)
seqlevels(peaks.down.d0.A_d2.A) <- correct_chromosome_names(seqlevels(peaks.down.d0.A_d2.A))

peakAnno.down.d0.A_d2.A <- annotatePeak(peaks.down.d0.A_d2.A, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d0.A_d2.A)
downDARs.d0.A_d2.A<-as.GRanges(peakAnno.down.d0.A_d2.A)
downDARs.d0.A_d2.A$geneId

# up DARs - contrast 2

pf.up.d0.A_d2.A <- atac.con[72]
peaks.up.d0.A_d2.A <- readPeakFile(pf.up.d0.A_d2.A)
seqlevels(peaks.up.d0.A_d2.A) <- correct_chromosome_names(seqlevels(peaks.up.d0.A_d2.A))

peakAnno.up.d0.A_d2.A <- annotatePeak(peaks.up.d0.A_d2.A, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d0.A_d2.A)
upDARs.d0.A_d2.A<-as.GRanges(peakAnno.up.d0.A_d2.A)
upDARs.d0.A_d2.A$geneId

##### Day 2 ACBI vs DMSO with ActA #####

# down DARs - contrast20

pf.down.d2.ACBI_DAct <- atac.con[13]
peaks.down.d2.ACBI_DAct <- readPeakFile(pf.down.d2.ACBI_DAct)

peakAnno.down.d2.ACBI_DAct <- annotatePeak(peaks.down.d2.ACBI_DAct, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d2.ACBI_DAct)
downDARs.d2.ACBI_DAct<-as.GRanges(peakAnno.down.d2.ACBI_DAct)
downDARs.d2.ACBI_DAct$geneId

write.table(downDARs.d2.ACBI_DAct$geneId, file = "downDARs.d2.ACBI_DAct.txt", sep="\t")

# up DARs - contrast20

pf.up.d2.ACBI_DAct <- atac.con[73]
peaks.up.d2.ACBI_DAct <- readPeakFile(pf.up.d2.ACBI_DAct)

peakAnno.up.d2.ACBI_DAct <- annotatePeak(peaks.up.d2.ACBI_DAct, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d2.ACBI_DAct)
upDARs.d2.ACBI_DAct<-as.GRanges(peakAnno.up.d2.ACBI_DAct)
upDARs.d2.ACBI_DAct$geneId

write.table(upDARs.d2.ACBI_DAct$geneId, file = "upDARs.d2.ACBI_DAct.txt", sep="\t")

##### Day 2 DMSO vs DMSO with ActA #####

# down DARs - contrast66

pf.down.d2.DMSO_DMSOAct <- atac.con[57]
peaks.down.d2.DMSO_DMSOAct <- readPeakFile(pf.down.d2.DMSO_DMSOAct)

peakAnno.down.d2.DMSO_DMSOAct <- annotatePeak(peaks.down.d2.DMSO_DMSOAct, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.down.d2.DMSO_DMSOAct)
downDARs.d2.DMSO_DMSOAct<-as.GRanges(peakAnno.down.d2.DMSO_DMSOAct)
downDARs.d2.DMSO_DMSOAct$geneId

write.table(downDARs.d2.DMSO_DMSOAct$geneId, file = "downDARs.d2.DMSO_DMSOAct.txt", sep="\t")

# up DARs - contrast66

pf.up.d2.DMSO_DMSOAct <- atac.con[117]
peaks.up.d2.DMSO_DMSOAct <- readPeakFile(pf.up.d2.DMSO_DMSOAct)

peakAnno.up.d2.DMSO_DMSOAct <- annotatePeak(peaks.up.d2.DMSO_DMSOAct, tssRegion=c(-3000, 3000), TxDb=edb, annoDb="org.Mm.eg.db")
as.GRanges(peakAnno.up.d2.DMSO_DMSOAct)
upDARs.d2.DMSO_DMSOAct<-as.GRanges(peakAnno.up.d2.DMSO_DMSOAct)
upDARs.d2.DMSO_DMSOAct$geneId

write.table(upDARs.d2.DMSO_DMSOAct$geneId, file = "upDARs.d2.DMSO_DMSOAct.txt", sep="\t")
