#!/usr/bin/env R script

##Mike Mariani UVM 2020 Frietze Lab

##example file name fro Drayman et al data
##After running dropseq pipeline:
##SRR8526684.1.star_gene_exon_tagged.dge.txt

library(grid)
library(gridExtra)
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)

set.seed(seed=1)

mock.dir <- "/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock"
wt.dir <- "/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt"
output.plots.dir <- "/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary"
dropseq_to_seurat <- function(files.in){
  library(Seurat)
  seurat.objects <- list()
  for(i in 1:length(files.in))
  {
    mat.frame <- read.table(file=files.in[i],
                           header = TRUE,
                           sep="\t", 
                           stringsAsFactors = F)
    rownames(mat.frame) <- mat.frame$GENE
    mat.frame <- mat.frame[,-1]
    ##mat <- matrix(mat.file, byrow = TRUE)
    seurat.objects[[i]] <- CreateSeuratObject(counts=mat.frame)
  }
  return(seurat.objects)
}

##read in hsv1 genes by kinetic class:
hsv1.genes.kc <- read.table(file = paste0(output.plots.dir,"/virus.genes.kc.txt"),
           header = FALSE,
           stringsAsFactors = TRUE)$V1
mock.files.in <- list.files(path=mock.dir,
                       pattern=".star_gene_exon_tagged.dge.txt",
                       full.names=TRUE)
mock.objects <- dropseq_to_seurat(files.in=mock.files.in)
merged.mock <- merge(x=mock.objects[[1]],y=mock.objects[c(2:length(mock.objects))])

##merged.object <- SCTransform(merged.object)
##store mitochondrial percentage in object meta data
merged.mock[["percent.mt"]] <- PercentageFeatureSet(merged.mock, pattern = "^MT-")
table(rownames(merged.mock) %in% hsv1.genes.kc)
##No virus counts present in mock
##merged.mock <- PercentageFeatureSet(merged.mock, features = hsv1.genes.kc, col.name = "percent.virus")
mock.qc.violin.plot <- VlnPlot(merged.mock, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
ggsave(filename=paste0(output.plots.dir,"/mock.qc.violin.plot.pdf"),
       plot=mock.qc.violin.plot,
       device="pdf",
       height=12,
       width=8)
mock.qc.corr.plot <- FeatureScatter(merged.mock,
                               feature1 = "nFeature_RNA",
                               feature2 = "nCount_RNA")
ggsave(filename=paste0(output.plots.dir,"/mock.qc.corr.plot.pdf"),
       plot=mock.qc.corr.plot,
       device="pdf",
       height=8,
       width=8)
##run sctransform
##merged.mock <- SCTransform(merged.mock, 
##                             vars.to.regress = "percent.mt", 
##                             verbose = FALSE)
##Above throws error with this
##100 gene dataset so let's try the traditional approach:
merged.mock <- NormalizeData(merged.mock)
all.genes <- rownames(merged.mock)
merged.mock <- ScaleData(merged.mock, features = all.genes)
merged.mock <- RunPCA(merged.mock, features = all.genes)
DimPlot(merged.mock, reduction="pca")

wt.files.in <- list.files(path=wt.dir,
                            pattern=".star_gene_exon_tagged.dge.txt",
                            full.names=TRUE)
wt.objects <- dropseq_to_seurat(files.in=wt.files.in)
merged.wt <- merge(x=wt.objects[[1]],y=wt.objects[c(2:length(wt.objects))])
##merged.object <- SCTransform(merged.object)
##store mitochondrial percentage in object meta data
merged.wt[["percent.mt"]] <- PercentageFeatureSet(merged.wt, pattern = "^MT-")
##merged.wt <- PercentageFeatureSet(merged.wt, features = hsv1.genes.kc, col.name = "percent.virus")
Seurat:::IsMatrixEmpty(GetAssayData(merged.wt, slot = 'counts'))
table(rownames(merged.wt) %in% hsv1.genes.kc)
##Where are the viral counts for PercentageFeatureSet??????
wt.qc.violin.plot <- VlnPlot(merged.wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
ggsave(filename=paste0(output.plots.dir,"/wt.qc.violin.plot.pdf"),
       plot=wt.qc.violin.plot,
       device="pdf",
       height=12,
       width=8)
wt.qc.corr.plot <- FeatureScatter(merged.wt,
                                    feature1 = "nFeature_RNA",
                                    feature2 = "nCount_RNA")
ggsave(filename=paste0(output.plots.dir,"/wt.qc.corr.plot.pdf"),
       plot=wt.qc.corr.plot,
       device="pdf",
       height=8,
       width=8)
##run sctransform
##merged.wt <- SCTransform(merged.wt, 
##                             vars.to.regress = "percent.mt", 
##                             verbose = FALSE)
##Above throws error with this
##100 gene dataset so let's try the traditional approach:
merged.wt <- NormalizeData(merged.wt)
all.genes <- rownames(merged.wt)
merged.wt <- ScaleData(merged.wt, features = all.genes)
merged.wt <- RunPCA(merged.wt, features = all.genes)
DimPlot(merged.wt, reduction="pca")

##Merge two datasets together:
merged.mock$orig.ident <- "mock"
merged.wt$orig.ident <- "wt"
merged.all <- merge(merged.mock, merged.wt)
merged.all <- NormalizeData(merged.all)
merged.all <- FindVariableFeatures(merged.all, selection.method = "vst", nfeatures = 2000)
all.genes  <- rownames(merged.all)
merged.all <- ScaleData(merged.all, features = all.genes)
merged.all <- RunPCA(merged.all, features = all.genes)
pca.dimplot <- DimPlot(merged.all, 
        reduction="pca",
        group.by = "orig.ident")
ggsave(filename = paste0(output.plots.dir,"/pca.dimplot.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=pca.dimplot)
qc.heatmap <- DimHeatmap(merged.all, 
           dims = c(1:10), 
           cells = 500, 
           balanced = TRUE)
ggsave(filename = paste0(output.plots.dir,"/qc.heatmap.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=qc.heatmap)
merged.all <- JackStraw(merged.all, num.replicate = 100)
merged.all <- ScoreJackStraw(merged.all, dims = 1:10)
jackstraw.plot <- JackStrawPlot(merged.all, dims = 1:10)
ggsave(filename = paste0(output.plots.dir,"/jackstraw.plot.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=jackstraw.plot)
elbow.plot <- ElbowPlot(merged.all)
ggsave(filename = paste0(output.plots.dir,"/elbow.plot.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=elbow.plot)
##Looks like about 4 PCs are meaningful let's try UMAP with 4 
##and (2x4=8) dims.
merged.all <- FindNeighbors(merged.all, dims = 1:4)
merged.all <- FindClusters(merged.all, resolution = 0.5)
merged.all <- RunUMAP(merged.all, dims = 1:4)
dimplot.umap.4pc <- DimPlot(merged.all, reduction = "umap")
ggsave(filename = paste0(output.plots.dir,"/dimplot.umap.4pc.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=dimplot.umap.4pc)
merged.all <- FindNeighbors(merged.all, dims = 1:8)
merged.all <- FindClusters(merged.all, resolution = 0.5)
merged.all <- RunUMAP(merged.all, dims = 1:8)
dimplot.umap.8pc <- DimPlot(merged.all, reduction = "umap")
ggsave(filename = paste0(output.plots.dir,"/dimplot.umap.8pc.pdf"),
       height=8,
       width=8,
       device="pdf",
       plot=dimplot.umap.8pc)

panel.mock.1 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526677.1.bam.kneeplot.pdf")
panel.mock.2 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526678.1.bam.kneeplot.pdf")
panel.mock.3 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526679.1.bam.kneeplot.pdf")
panel.mock.4 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526681.1.bam.kneeplot.pdf")
panel.mock.5 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526682.bam.kneeplot.pdf")
panel.mock.6 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526683.1.bam.kneeplot.pdf")
panel.mock.7 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/mock/SRR8526684.1.bam.kneeplot.pdf")

mock.kneeplots.together <- grid.arrange(rasterGrob(panel.mock.1),
                           rasterGrob(panel.mock.2),
                           rasterGrob(panel.mock.3),
                           rasterGrob(panel.mock.4),
                           rasterGrob(panel.mock.5),
                           rasterGrob(panel.mock.6),
                           rasterGrob(panel.mock.7),
                           ncol=2)
ggsave(filename = paste0(output.plots.dir,"/mock.kneeplots.together.pdf"),
       device="pdf",
       height=16,
       width=8,
       plot=mock.kneeplots.together)

panel.wt.1 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt/SRR8526681.1.kneeplot.pdf")
panel.wt.2 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt/SRR8526682.kneeplot.pdf")
panel.wt.3 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt/SRR8526683.1.kneeplot.pdf")
panel.wt.4 <- image_read_pdf("/slipstream/home/mmariani/projects/hsv1_scrna/drop_seq/drayman/summary/wt/SRR8526684.1.kneeplot.pdf")

wt.kneeplots.together <- grid.arrange(rasterGrob(panel.wt.1),
                         rasterGrob(panel.wt.2),
                         rasterGrob(panel.wt.3),
                         rasterGrob(panel.wt.4),
                         ncol=2)

ggsave(filename = paste0(output.plots.dir,"/wt.kneeplots.together.pdf"),
       device = "pdf",
       height = 8,
       width = 8,
       plot = wt.kneeplots.together)
