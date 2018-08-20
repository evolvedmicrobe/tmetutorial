
## 1 Introduction

# Dissecting the composition of non-malignant cells within a tumor is key to
# understanding the interaction between a tumor and its microenvironment,
# predicting the clinical outcome, and assisting in the selection of therapies.
# In the Application Note, [*Characterization of the Tumor
# Microenvironment*](https://www.10xgenomics.com/resources/application-notes/),
# we showed how the Chromium Single Cell Immune Profiling Solution was applied
# to a colorectal cancer (CRC) tumor and a squamous cell non-small cell lung
# carcinoma (NSCLC).

# The bioinformatics community is also actively developing software to analyze
# Chromium Single Cell data, which meet the needs for customized secondary
# analysis towards various biological questions. Here we show how the [NSCLC
# data](https://support.10xgenomics.com/single-cell-vdj/datasets) can be
# analyzed in R with [Seurat](http://satijalab.org/seurat/) and
# [Monocle](http://cole-trapnell-lab.github.io/monocle-release/), which are two
# popular packages for single cell gene expression analysis. Through this
# analysis, we show how to apply Seurat and Monocle to 10x Genomics data, how to
# link gene expression with immune receptor sequencing data for the same single
# cells, and how to communicate results from different packages/analytic
# modules.

## 2 Clustering cells

### 2.1 Load 10x data as a Seurat Object

# We started with the filtered gene-barcode matrices files from the outputs of
# [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger).
# The files used in this tutorial are available
# [here](http://cf.10xgenomics.com/samples/cell-vdj/2.1.0/vdj_v1_nsclc_5gex/vdj_v1_nsclc_5gex_filtered_gene_bc_matrices.tar.gz).
# They have already been downloaded and extracted into the local directory. Cell
# Ranger also provides a set of raw matrix files which include UMI counts from
# both cell and background barcodes. When dealing with highly heterogeneous
# samples (e.g. tumor samples), it is important to examine whether the default
# cell-barcode calling method is optimal (see [Cell Ranger
# AlgorithmsOverview](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/overview)
# for details). The raw matrix files could be useful when some cell barcodes
# appear to fall below the cutoff.

# We start by loading the matrix data as a Seurat object. The Seurat function
# `CreateSeuratObject` also allows for filtering out cells and genes with too
# many zero-UMIs. In this step, 13,463 (40%) genes and 212 (2.7%) cells were
# filtered out.
library(tmetutorial)
library(Seurat)
library(dplyr)
library(Matrix)
library(formatR)
# Load the NSCLC dataset, we've distributed this with R.
gex.data <- Read10X(data.dir = tmetutorial::getSingleCellDataPath())
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 100 detected genes.
gex <- CreateSeuratObject(raw.data = gex.data, min.cells = 3,
                            min.genes = 100, project = "nsclc_5gex")

### 2.2 Standard workflow in Seurat

# The steps below encompass the standard workflow for cell clustering analysis
# based on scRNA-seq data in Seurat, including selection and filtration of cells
# based on QC metrics, data normalization and scaling, detection of highly
# variable genes, dimension reduction, graph-based clustering, and cluster
# visualization. For more details of Seurat workflow, please refer [Seurat's
# tutorial on 10x PBMC data](http://satijalab.org/seurat/pbmc3k_tutorial.html).
# In this tutorial, we only highlight parameters that are customized for tumor
# microenvironment (TME) analysis.

#### 2.2.1 Data normalization

# Here, we plot the distribution of number of genes and UMI counts per cell, and
# fraction of UMIs mapping to mitochondrial genes.

mito.genes <- grep(pattern = "^MT-",
                   x = rownames(x = gex@data), value = TRUE)
percent.mito <-
  Matrix::colSums(gex@raw.data[mito.genes, ]) /
  Matrix::colSums(gex@raw.data)

gex <- AddMetaData(object = gex, metadata = percent.mito,
                     col.name = "percent.mito")
VlnPlot(object = gex,
        features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = gex, gene1 = "nUMI", gene2 = "percent.mito",
         cex.axis = 2, cex.lab = 2)
GenePlot(object = gex, gene1 = "nUMI", gene2 = "nGene",
         cex.axis = 2, cex.lab = 2)

# Based on the plot, we set cutoff to remove cells with abnormal metrics.
# Considering the high heterogeneity of the TME, we allow a wide range of number
# of genes/UMIs per cell. In this case, 287 out of 7,590 cells were filtered
# out. number of filtered cells

sum(gex@meta.data$nGene > 8000 | gex@meta.data$percent.mito > 0.2)
gex <- FilterCells(object = gex, subset.names = c("nGene", "percent.mito"),
                     low.thresholds = c(100, -Inf),
                     high.thresholds = c(8000, 0.2))
gex <- NormalizeData(object = gex, normalization.method = "LogNormalize",
                       scale.factor = 1e4, display.progress = FALSE)

#### 2.2.2 Detection of variable genes across the single cells

# Instead of using all genes, focusing on variable genes increases
# signal-to-noise ratio of the downstream analysis. We here used the
# Seurat-suggested parameters for UMI data that are normalized to a total of 1e4
# molecules, which identified 3199 variable genes.

gex <- FindVariableGenes(object = gex, mean.function = ExpMean,
                           dispersion.function = LogVMR, x.low.cutoff = 0.0125,
                           x.high.cutoff = 3, y.cutoff = 0.5,
                           do.plot = FALSE, display.progress = FALSE)
length(x = gex@var.genes)

#  We adjusted the parameters to generate different numbers of variable genes
#  (2000-4000), and ran the downstream analysis. This resulted in similar cell
#  clustering patterns (not shown here). The number of genes needed will vary
#  depending on the sample, and it will often be helpful to adjust the
#  parameters based on a prior expected degree of heterogeneity and variation of
#  the cells in the sample.

#### 2.2.3 Dimensional reduction

# Principal component analysis (PCA) was used for linear dimensional reduction.
# Here, we select to compute 100 principal components (PCs) instead of the
# default 20 PCs. Using more PCs allows the identification of small clusters
# within a heterogeneous sample.

gex <- ScaleData(object = gex, vars.to.regress = c("nUMI", "percent.mito"),
                   display.progress = FALSE, do.par = TRUE, num.cores = 2)
gex <- RunPCA(object = gex, pcs.compute = 100, pc.genes = gex@var.genes, do.print = FALSE)

# The heatmaps below show the 3 top, middle and bottom PCs.

PCHeatmap(object = gex, pc.use = c(1:3, 48:50, 98:100), cells.use = 500,
          do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)


# Next, we determine the statistically significant principal components.
gex <- JackStraw(object = gex, num.pc = 100,
                   display.progress = FALSE,
                   do.par = TRUE, num.cores = 2)

# The `JackStrawPlot` function in Seurat provides a visualization tool for
# comparing the distribution of p-values for each PC with a uniform distribution
# (dashed line). 'Significant' PCs will show a strong enrichment of genes with
# low p-values (solid curve above the dashed line). In this case PCs 1-58 are
# all significant and the first insignificant PC is PC59.

JackStrawPlot(object = gex, PCs = 1:60, nCol = 6)

# Next, we viewed the plot of the standard deviations of the principal
# components with `PCElbowPlot` and try to find a cutoff. We selected the first
# 50 PCs to be used in the following analysis since it appeared to represent a
# significant portion of variation.

PCElbowPlot(object = gex, num.pc = 100)

#### 2.2.4 Cluster the cells

# Similar to Cell Ranger, the graph-based clustering algorithm as implemented in
# Seurat's `FindClusters` function was adopted to cluster the cells. The number
# of resulting cell clusters could be adjusted via the parameter `resolution`.
# `resolution = 1` was used here according to the instructions from Seurat. We
# also tried `resolution = 2`, and found four more clusters were divided from
# the `resolution = 1` clusters, but their biological identities were not
# clearly distinguished (data not shown).

gex <- FindClusters(object = gex, reduction.type = "pca", dims.use = 1:50,
                      resolution = 1, print.output = 0, save.SNN = TRUE)

# We visualize cell clusters with tSNE.
gex <- RunTSNE(object = gex, dims.use = 1:50, do.fast = TRUE)
TSNEPlot(object = gex, do.label = TRUE)

# We can save the Seurat object, so the analysis can be resumed here without re-running the steps above.
saveRDS(gex, file = "nsclc_tutorial.pc50.noLabel.rds")
# To resume the analysis here, comment the above line and un-comment the line below:
# gex <- readRDS("nsclc_tutorial.pc50.noLabel.rds")

### 2.3 Label clusters

#### 2.3.1 Finding genes that are enriched in a given cluster (cluster biomarkers)

#  We find markers for every cluster compared to cells in all other clusters and
#  only report the genes with higher expression levels within each cluster.
#  `min.pct = 0.25` was used to ensure that the cells in each cluster were well
#  represented by the identified markers.

gex.markers <- FindAllMarkers(object = gex, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25, print.bar = FALSE)
saveRDS(gex.markers, "nsclc_tutorial.gex.markers.rds")

# View the top 2 markers for each cluster in heatmap.
gex.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

#### 2.3.2 Label each cluster based on expression of marker genes

# We take two steps to give each cluster a unique and meaningful name. First, we
# examine the expression of well characterized marker genes to identify major
# cell types that we expect to be associated with tumor samples. In the context
# of the TME, we expect to see tumor cells, lymphocytes (B and T cells) and
# myeloid-derived cells. Here we are looking at B lymphocytes (MS4A1), T
# lymphocytpes (CD3E, CD8A, FOXP3), NK cells (GNLY), and monocytes (CD14,
# FCER1A, FCGR3A, LYZ).

# enable FeaturePlot with cluster identity labels
labelFeaturePlot<- function(object, features.plot, nCol = 2,
                            cols.use = c("grey", "blue"),
                            reduction.use = "tsne", ...) {
  pLabel <- as.data.frame(
    GetDimReduction(object, reduction.type = reduction.use, slot = "cell.embeddings")
  )
  pLabel$ident <- object@ident
  centers <- pLabel %>% dplyr::group_by(ident) %>%
    dplyr::summarize(x = median(tSNE_1),
                     y = median(tSNE_2))
  featureList <- features.plot[which(
    features.plot %in% object@data@Dimnames[[1]]
  )]
  notFound <- features.plot[which(
    !features.plot %in% object@data@Dimnames[[1]]
  )]
  if (length(notFound) > 0) {
    for (nf in notFound) {
      warning(paste0("Cannot find features: "), nf)
    }
  }
  pdf(NULL)
  p <- FeaturePlot(object = object, features.plot = featureList,
                   cols.use = cols.use, reduction.use = reduction.use,
                   do.return = TRUE, ...)
  dev.off()
  cowplot::plot_grid(plotlist = lapply(p, function(x)
    x + geom_text(data = centers, mapping = aes(label = ident))),
    ncol = nCol)
}

# Expression of known markers
labelFeaturePlot(gex, features.plot = c(
  "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
  "FCGR3A", "LYZ", "CD8A", "FOXP3"
), nCol = 3)

# Second, we try to determine the identity of each cluster within the main
# types, based on genes specifically expressed in the cells within each cluster.
# The ideal marker genes should show a clear difference of fractions of cells
# expressing them within different clusters. Based on prior knowledge, a proper
# cluster name representing cellular identity is determined for each cell
# cluster. This is a laborious process, which may not always lead to a
# well-studied cell type, but with a great potential for interesting
# observations and discoveries. Here, we show the identification of Exhausted T
# cells as an example. Based on the expression of CD3E and CD8A, we know
# Cytotoxic T cells should be among Clusters 3, 11 and 18. Then we check the
# genes that are specific to these clusters.

gex.markers %>%
  filter(cluster %in% c(3, 11, 18) & pct.1 - pct.2 > 0.5) %>%
  group_by(cluster) %>%
  top_n(5, avg_logFC)

# CXCL13 is a marker for Exhausted T cells ([Zheng et al,
# 2017](http://dx.doi.org/10.1016/j.cell.2017.05.035)). It was expressed in 79%
# of cells in Cluster 18 but <2% of the rest of the cells. We then view the
# expression of CXCL13 as well as other well-known markers for Exhausted T
# cells, CTLA4, PDCD1 and TIGIT, and they all support the identify of this
# cluster as Exhausted T cells.

labelFeaturePlot(object = gex,
                 features.plot = c("CXCL13", "CTLA4", "PDCD1", "TIGIT"))

# In this case, we iterated this procedure to the other individual clusters, and
# the cluster identities and their marker genes were listed as below. We
# assigned biologically meaningful names to clusters that were well-defined by
# marker genes (e.g. T helper cells). For the other clusters, we assigned the
# names by adding 1-2 of their marker genes to broader, well-defined cell types,
# and also provided related literatures that showed functional roles of these
# marker genes in cancer.

#  Cluster ID | Markers       | Cell Type                         |Reference
#  -----------|---------------|-----------------------------------|---------
#  0          | MS4A1         | B cells                           |
#  1          | IL7R          | T helper cells                    |
#  2          | SERPINB3      | SerpinB3+ Epithelial/tumor cells  | [Catanzaro et al, 2014](https://www.nature.com/articles/ncomms4729)
#  3          | CD8A          | Cytotoxic T cells                 |
#  4          | CD14, GPNMB   | GPNMB+ Macrophages                | [Szulzewsky et al, 2015](https://doi.org/10.1371/journal.pone.0116644)
#  5          | FCN1, S100A12 | FCN1+ Macrophages                 | [Honore et al, 2008](https://doi.org/10.1016/j.molimm.2008.02.005)
#  6          | S100A7        | S100A7+ Epithelial/tumor cells    | [Zhang et al, 2008](http://thorax.bmj.com/content/63/4/352)
#  7          | FOXP3         | Tregs                             |
#  8          | S100A2        | S100A2+ Epithelial/tumor cells    |
#  9          | GNLY, NKG7    | NK cells                          |
#  10         | MS4A2, TPSB2  | Mast cells                        | [Drzewiecka et al, 2013](https://link.springer.com/article/10.1007/s00251-013-0695-8)
#  11         | TRAV1-2       | MAIT cells                        |
#  12         | MZB1, IG*     | Plasma cells                      |
#  13         | FCER1A        | Dendritic cells                   |
#  14         | CXCL8, G0S2   | Monocytes                         | [Huang et al, 2017](https://doi.org/10.1016/j.phrs.2017.04.020)
#  15         | MT*           | Dying cells                       |
#  16         | TCL1A         | TCL1A+ B cells                    | [Punt et al, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4770729/)
#  17         | CD1A, S100B, FCER1A | CD1A+ Dendritic cells       | [Solerlund et al, 2004](http://onlinelibrary.wiley.com/doi/10.1111/j.1600-0684.2003.00053.x/full)
#  18         | CXCL13, CTLA4 | Exhausted T cells                 | [Zheng et al, 2017](http://dx.doi.org/10.1016/j.cell.2017.05.035)
#  19         | IFIT1, IFIT3  | IFIT1+ B cells                    | [Li et al, 2009](http://www.pnas.org/content/106/19/7945.short)
#  20         | FCRL4         | FCRL4+ memory B cells             | [Jourdan et al, 2017](https://doi.org/10.1371/journal.pone.0179793)
#  21         | UBE2C, STMN1  | UBE2C+STMN1+ Epithelial/tumor cells | [Hao et al, 2012 ](https://link.springer.com/article/10.1007/s13277-011-0291-1), [Nie et al, 2015](https://www.nature.com/articles/labinvest2014124)
#  22         | GZMB          | GZMB+ B cells                     | [Hagn et al, 2010](http://onlinelibrary.wiley.com/doi/10.1002/eji.200940113/full)

current.cluster.ids <- 0:22
new.cluster.ids <- c(
  "B cells", "T helper cells", "SerpinB3+ Epithelial/tumor cells",
  "Cytotoxic T cells", "GPNMB+ Macrophages", "FCN1+ Macrophages",
  "S100A7+ Epithelial/tumor cells", "Tregs", "S100A2+ Epithelial/tumor cells",
  "NK cells", "Mast cells", "MAIT cells",
  "Plasma cells", "Dendritic cells", "Monocytes",
  "Dying cells", "TCL1A+ B cells", "CD1A+ Dendritic cells",
  "Exhausted T cells", "IFIT1+ B cells", "FCRL4+ memory B cells",
  "UBE2C+STMN1+ Epithelial/tumor cells", "GZMB+ B cells"
)
gex@ident <- plyr::mapvalues(x = gex@ident,
                             from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = gex, do.label = TRUE, no.legend = TRUE)

# From the `DotPlot` of the marker genes, we can get a sense of the specificity
# of each marker toward their representing clusters. A cluster is labeled as
# "Dying cells" because it had high expression of mitochondrial genes. Taking a
# closer look, we can infer the original identities of these cells by the
# presence of marker genes of several other cell types, such as B cells (MS4A1),
# T helper cells (IL7R) and epithelial/tumor cells (SERPINB3, S100A2, S100A7).
DotPlot(gex,
        genes.plot = c(
          "MS4A1", "IL7R", "SERPINB3", "CD8A", "CD14", "GPNMB",
          "FCN1", "S100A12", "S100A7", "FOXP3", "S100A2", "GNLY",
          "NKG7", "MS4A2", "TPSB2", "TRAV1-2", "MZB1", "FCER1A",
          "CXCL8", "G0S2", "TCL1A", "CD1A", "S100B", "CXCL13",
          "CTLA4", "IFIT1", "IFIT3", "FCRL4", "UBE2C", "STMN1",
          "GZMB"
        ),
        plot.legend = TRUE, x.lab.rot = TRUE)

# We can also tabulate the statistics of cellular composition from the Seurat
# object `ident`. The assignment of the cell identities were subjective, and
# further experiments are usually needed to assess the functionality of the cell
# clusters.
clusterCount <- table(gex@ident)
statCluster <- data.frame(cluster = names(clusterCount), count = as.numeric(clusterCount))
statCluster$Fraction <- round(statCluster$count / sum(clusterCount), 2)
statCluster <- statCluster[order(-statCluster$count), ]
rownames(statCluster) <- 1:nrow(statCluster)
statCluster

# Using a proper number of principal components is important for achieve
# informative clutering. Therefore, we tested how cells were clustered based on
# 21 (around the elbow of the `PCElbowPlot`) or 75 PCs, and compared to result
# from 50 PCs as shown above. The same tSNE projection as derived with 50 PCs
# was used, and clusters were differentiated by colors. The cells were clustered
# similarly between 50 and 75 PCs, whereas an obvious difference was seen
# between 21 and 50 PCs.
gex <- StashIdent(gex, save.name = "labelPC50")
gex <- FindClusters(object = gex, reduction.type = "pca", dims.use = 1:21,
                    resolution = 1, print.output = 0, save.SNN = TRUE)
gex <- StashIdent(gex, save.name = "PC21")
gex <- FindClusters(object = gex, reduction.type = "pca", dims.use = 1:75,
                    resolution = 1, print.output = 0, save.SNN = TRUE)
gex <- StashIdent(gex, save.name = "PC75")
gex <- FindClusters(object = gex, reduction.type = "pca", dims.use = 1:50,
                    resolution = 1, print.output = 0, save.SNN = TRUE)
gex <- StashIdent(gex, save.name = "PC50")
plotPc21 <- TSNEPlot(
  gex, do.return = TRUE, no.legend = TRUE, group.by = "PC21", do.label = TRUE
) + ggtitle("21 PCs")
plotPc50 <- TSNEPlot(
  gex, do.return = TRUE, no.legend = TRUE, group.by = "PC50", do.label = TRUE
) + ggtitle("50 PCs")
plotPc75 <- TSNEPlot(
  gex, do.return = TRUE, no.legend = TRUE, group.by = "PC75", do.label = TRUE
) + ggtitle("75 PCs")
gex <- SetAllIdent(gex, id = "labelPC50")
plot_grid(plotPc50, plotPc21, plotPc75)

# Notably, in the 21-PC result, B cells were almost evenly divided into two
# clusters. We tried to test if there was any marker gene that could inform
# their identity. However, no gene was found expressed preferentially in either
# cluster, which suggested the division of them was not valid. On the other
# hand, compare with the 21-PC result, more clusters with distinguishable
# markers were identified in the 50-PC result.

gex <- SetAllIdent(gex, id = "PC21")
bmarkers <- FindMarkers(gex, ident.1 = 0, ident.2 = 1,
                        print.bar = F)
gex <- SetAllIdent(gex, id = "labelPC50")
bmarkers

## 3 Cell Trajectory Analysis

### 3.1 Import the Seurat object to Monocle

# Transcriptome profiles of single cells reveal differentiation states of each
# cell. *Pseudotime* is commonly used to quantify the progress of cells during
# the transition between differentiation states. It could be considered as an
# aggregation of the expression of a set of genes that change from one
# differentiation state to another. [Monocle
# 2](http://cole-trapnell-lab.github.io/monocle-release) is one of the popular
# packages developed to address this problem. As an example, we use Monocle 2 to
# analyze the trajectory of T cell exhaustion.

# To define genes that order the cells along the trajectory, we need to first
# identify the cells that represent the different states. These cells could be
# harvested from multiple experiments in different conditions (e.g., time
# points, cultured medium) that align to the presumed trajectory. For cells
# collected in a single experiment, as in this case, some clusters could
# represent cells in different states under proper assumptions. Because we
# already performed cell clustering analysis above, rather than starting over
# from the original 10x gene-barcode matrix, we import the Seurat object to
# Monocle using the function `importCDS`.

library(monocle)
HSMM <- importCDS(gex, import_all = TRUE)
colnames(HSMM@phenoData@data)[5] <- 'Cluster'
HSMM@phenoData@data$Cluster <- gex@ident
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

### 3.2 Inferring cell trajectories

# As shown above, the exhausted T cells were CD8 positive, so their earlier
# states could present as cytotoxic T cells. Thus, we identify genes that were
# differentially expressed between cytotoxic T cells and exhausted T cells as
# markers for exhaustion. The exhausted T cells which showed specific expression
# of several marker genes (outlined above),  could be considered as a more
# homogeneous population than cytotoxic T cells, which were more likely a
# mixture of T cells in terms of exhaustion state. Therefore, by applying
# `only.pos = TRUE` in `FindMarkers`, we only use genes that were specifically
# enriched in exhausted T cells to order the cells. Although we did not compare
# T helper cells with regulatory T cells (Tregs), they formed a branched
# trajectory apart from the branch of exhausted T cells, consistent with the
# observation that some T cell exhaustion markers (eg. CTLA4, PDCD1) were also
# expressed in the Tregs cluster.

exh.markers <- FindMarkers(object = gex, ident.1 = c("Exhausted T cells"),
                           ident.2 = c("Cytotoxic T cells"), only.pos = TRUE,
                           min.pct = 0.25, logfc.threshold = 0.5,
                           print.bar = FALSE)
HSMM_ordering_genes <- unique(row.names(exh.markers))
tcellsIdx <- which(gex@ident %in% c("T helper cells", "Cytotoxic T cells", "Tregs",
                                    "MAIT cells", "Exhausted T cells"))
HSMM_t <- HSMM[, tcellsIdx]
HSMM_t <- setOrderingFilter(HSMM_t, ordering_genes = HSMM_ordering_genes)
HSMM_t <- reduceDimension(HSMM_t, method = 'DDRTree')
HSMM_t <- orderCells(HSMM_t)
clusterState <- table(HSMM_t@phenoData@data$Cluster, HSMM_t@phenoData@data$State)

library(reshape2)
clusterState <- melt(clusterState, varnames = c("Cluster", "State"))
colnames(clusterState)[3] <- 'Count'
clusterState$Cluster <- as.character(clusterState$Cluster)
clusterState <- filter(clusterState, Count > 0)
p1 <- ggplot(clusterState, aes(factor(State), Count, fill = Cluster)) +
  geom_bar(stat = "identity")
p2 <- plot_cell_trajectory(HSMM_t, color_by = "Cluster") +
  theme(legend.position = "none")
library(Rmisc)
multiplot(p1, p2, cols = 2)

## 4. Analysis of the immune repertoire

# The 10x Genomics Single Cell Immune Profiling Solution allows the clonotypes
# of the T-cell receptor (TCR) and B cell immunoglobulins (Ig) from the same
# sample to be determined, which enables us to trace the lineages of the immune
# cells. In this case, the corresponding TCR and Ig clonotype data can be
# obtained from https://support.10xgenomics.com/single-cell-vdj/datasets. We can
# then intersect the repertoire sequencing data with the gene expression data.

### 4.1 TCR clonotypes

# We load the T cell clonotype data from the
# [`*filtered_contig_annotations.csv`](http://cf.10xgenomics.com/samples/cell-vdj/2.1.0/vdj_v1_nsclc_t/vdj_v1_nsclc_t_filtered_contig_annotations.csv)
# file. The clonotypes were annotated by Cell Ranger and unique clonotype IDs
# were assigned. We then counted the frequencies of cell barcodes for each
# clonotype and examined the identities of cells supporting the TCR clonotypes
# shared by multiple cells. The cells in the GeX and V(D)J data was linked by
# the cell barcodes, which consists of a 16-nt string and a GEM group number
# including a in both GeX and V(D)J data.

# In this case we noticed that:

# * Most of the expanded clonotypes were CD8+ T cells, especially exhausted T cells
# * Most of the expanded clonotypes were within single cell types
# * The same clonotypes were shared between regular and exhausted T cells

# load TCR clonotype data
tcontig <- read.csv(
  "http://cf.10xgenomics.com/samples/cell-vdj/2.1.0/vdj_v1_nsclc_t/vdj_v1_nsclc_t_filtered_contig_annotations.csv"
)

# Seurat removes the GEM group numbers when they is only one GEM group.
# We do the same thing here to coordinate with the GeX data.
# Do skip this step when working on an aggregated data
tcontig$barcode <- unlist(lapply(as.character(tcontig$barcode),
                                 function (x) strsplit(x, "-")[[1]][1]))
tcontig <- tcontig %>%
  select(barcode, raw_clonotype_id) %>%
  filter(raw_clonotype_id != "None") %>%
  distinct()
# load cell identity data
tcontig$ident <- gex@ident[tcontig$barcode]
tcontig <- filter(tcontig,
                  ident %in% c("T helper cells", "Cytotoxic T cells",
                               "Exhausted T cells", "MAIT cells", "Tregs"))
tcontig$raw_clonotype_id <- as.character(tcontig$raw_clonotype_id)
ctFreq <- table(tcontig$raw_clonotype_id)
tcontig$ctFreq <- as.numeric(ctFreq[tcontig$raw_clonotype_id])
tcontig %>%
  filter(ctFreq >= 10) %>%
  group_by(raw_clonotype_id) %>%
  top_n(3, barcode) %>%
  arrange(desc(ctFreq))
highFreqCT <- unique(names(ctFreq[which(ctFreq > 2)]))
tcontig$raw_clonotype_id <- as.factor(tcontig$raw_clonotype_id)
ggplot(tcontig[which(tcontig$raw_clonotype_id %in% highFreqCT), ],
       aes(reorder(raw_clonotype_id), fill = ident)) + geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# We can also view the clonotypes along the differentiation trajectory. The colored clonotypes are those that were shared by at least three cells.
HSMM_t@phenoData@data[tcontig$barcode, 'clonotype'] <-
  as.character(tcontig$raw_clonotype_id)
HSMM_t@phenoData@data <- HSMM_t@phenoData@data[dimnames(HSMM_t@reducedDimS)[[2]], ]
HSMM_t@phenoData@data[
  which(HSMM_t@phenoData@data$clonotype %in% highFreqCT),
  'expandedClonotype'
  ] <- as.character(HSMM_t@phenoData@data[
    which(HSMM_t@phenoData@data$clonotype %in% highFreqCT),
    'clonotype'
    ])
HSMM_t@phenoData@data <- HSMM_t@phenoData@data[order(
  HSMM_t@phenoData@data$expandedClonotype,
  decreasing = TRUE, na.last = FALSE
), ]
plot_cell_trajectory(HSMM_t, color_by = "expandedClonotype") +
  theme(legend.text = element_text(size = 5), legend.position = "none")


### 4.2 B cell Ig clonotypes

#  Similarly, we can also intersect Ig clonotypes with B cell clusters. We see a
#  high proportion of plasma and FCRL4+ memory B cells with clonotypes shared
#  with other cells. Conversely, a very low fraction of GZMB+ B cells have
#  detected receptor sequences.

# load BCR clonotype data
bcontig <- read.csv(
  "http://cf.10xgenomics.com/samples/cell-vdj/2.1.0/vdj_v1_nsclc_b/vdj_v1_nsclc_b_filtered_contig_annotations.csv"
)
bcontig$barcode <- unlist(lapply(as.character(bcontig$barcode),
                                 function (x) strsplit(x, "-")[[1]][1]))
bcontig <- bcontig %>%
  select(barcode, raw_clonotype_id) %>%
  filter(raw_clonotype_id != "None") %>%
  distinct()
# load cell identity data
bcontig$ident <- gex@ident[bcontig$barcode]
bcontig <- filter(bcontig,
                  ident %in% c(
                    "B cells", "Dying cells", "FCRL4+ memory B cells",
                    "TCL1A+ B cells", "Plasma cells", "IFIT1+ B cells", "GZMB+ B cells"
                  ))
ctFreq <- table(bcontig$raw_clonotype_id)
bcontig$ctFreq <- as.numeric(ctFreq[bcontig$raw_clonotype_id])
fracBcrExp <- table(bcontig$ident) / table(gex@ident)
fracBcrExp <- data.frame(cluster = names(fracBcrExp),
                         fraction = as.numeric(fracBcrExp),
                         stringsAsFactors = FALSE)
fracBcrExp$SharedBetweenCells <- as.numeric(
  table(bcontig[bcontig$ctFreq > 1, 'ident']) / table(gex@ident)
)
fracBcrExp$Unique <- fracBcrExp$fraction - fracBcrExp$SharedBetweenCells
fracBcrExp %>% filter(fraction > 0) %>% arrange(desc(fraction)) -> fracBcrExp
fracBcrExp <- within(
  fracBcrExp,
  cluster <- factor(cluster, levels = cluster)
)
fracBcrExp <- melt(fracBcrExp[, c(1, 3, 4)], id.vars = 'cluster')
colnames(fracBcrExp)[2:3] <- c('clonotype', 'fraction')
ggplot(fracBcrExp, aes(cluster, fraction, fill = clonotype)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# We can overlay the frequencies of clonotypes on the tSNE plots.
# Add the clonotype frequency onto the Seurat object
clonotypeFreq <- c()
clonotypeFreq[gex@data@Dimnames[[2]]] <- NA
clonotypeFreq[tcontig$barcode] <- tcontig$ctFreq
clonotypeFreq[bcontig$barcode] <- bcontig$ctFreq
gex@meta.data$clonotypeFreq <- clonotypeFreq
gex@meta.data$clonotype <- "None"
table(gex@meta.data$clonotypeFreq)
gex@meta.data$clonotype[which(gex@meta.data$clonotypeFreq == 1)] <- "Unique"
gex@meta.data$clonotype[which(gex@meta.data$clonotypeFreq > 1)] <- "SharedBetweenCells"
gex@meta.data$clonotype[which(gex@meta.data$clonotypeFreq > 4)] <- "HighFreq"
colors <- c(HighFreq = "red",
            SharedBetweenCells = "orange",
            Unique = "lightgreen",
            None = "lightgrey")
DimPlot(object = gex, reduction.use = "tsne", do.label = FALSE,
        cols.use = colors, do.hover = TRUE, group.by = "clonotype",
        data.hover = c("ident", "clonotypeFreq"))

## 5. Session Info
sessionInfo()
save.image("nsclc8k_image.Rda")


