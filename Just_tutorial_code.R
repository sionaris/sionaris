# Print the working directory
getwd()

# Optional:
# Set the working directory to a directory you prefer using an absolute path
# e.g. setwd("C:/Users/username/Desktop")

setwd("C:/Users/as3582/OneDrive - University of Cambridge/Desktop/Cambridge/AUTH/IAB06_2025/Tutorial")

# Install and load libraries #####

# Function to check and install CRAN packages
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Function to check and install Bioconductor packages
install_if_missing_bioc <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
  }
}

# List of CRAN packages
cran_packages <- c("dplyr", "ggplot2", "pheatmap", "PoiClaClu", "RColorBrewer", 
                   "glmpca", "ggbeeswarm", "hexbin", "ggpubr", "cowplot",
                   "rcartocolor", "ggrepel", "openxlsx", "ashr", "sva", "pathview")

# List of Bioconductor packages
bioc_packages <- c("airway", "DESeq2", "vsn", "org.Hs.eg.db", "EnhancedVolcano",
                   "apeglm", "clusterProfiler", "ReactomePA")

# Install CRAN packages
lapply(cran_packages, install_if_missing)

# Install Bioconductor packages
lapply(bioc_packages, install_if_missing_bioc)

rm(cran_packages, bioc_packages, install_if_missing, install_if_missing_bioc)

# Load packages
library(airway)
library(sva)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(glmpca)
library(EnhancedVolcano)
library(cowplot)
library(openxlsx)
library(vsn)
library(tibble)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(rcartocolor)
library(PoiClaClu)
library(clusterProfiler)
library(ReactomePA)
library(pathview)

# Clean up unused memory
gc()

# Importing data and preprocessing #####
data("airway")
se <- airway
View(se)

# Use DESeqDataSet() to create a DESeq2 object
# the model design will include treatment and cell line origin
dds <- DESeqDataSet(se, design = ~ cell + dex) 
dds

# Filter the DESeq2 object
# Keep genes with a sum of at least 10 counts
min.count = 10
keep <- rowSums(counts(dds)) >= min.count

# Keep genes for which at least 3 samples with a count of 10 or higher
min.sample = 3
keep <- rowSums(counts(dds) >= min.count) >= min.sample
dds <- dds[keep, ]
dds

# Convert dex to factor and set levels
dds$dex <- factor(dds$dex, levels = c("untrt","trt"), labels = c("untrt", "trt"))
dds$dex <- relevel(dds$dex, ref = "untrt")

# Exploratory analysis #####

# Estimating size factors
dds = estimateSizeFactors(dds)

# Original counts
head(counts(dds))

# Normalized counts
head(counts(dds, normalized  = TRUE))

# Variance stabilizing transformation: vst
# https://github.com/thelovelab/DESeq2/blob/devel/inst/script/vst.pdf
vsd <- vst(dds, blind = FALSE)

# Why blind = FALSE? 
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization

# regularized Log transformation: rlog
# https://doi.org/10.1186/s13059-014-0550-8
rld = rlog(dds, blind = FALSE)

# Simple log2(counts + 1) transformation
ntd = normTransform(dds) # log2(counts + 1) transformation

# Mean - sd plots
# ntd
meanSdPlot(assay(ntd))

# rlog
meanSdPlot(assay(rld))

# vst
meanSdPlot(assay(vsd))

# Exploratory heatmaps ###

# Prepare plot data frame
df = as.data.frame(colData(dds)[,c("cell", "dex")])
df$cell = factor(df$cell)
df$dex = factor(df$dex)

# We are going to use colorblind-friendly palettes here
# You can review them in https://jakubnowosad.com/rcartocolor/

# Annotation colors
ann_colors = list(
  cell = c(N61311 = carto_pal(n = 12, "Bold")[2], 
           N052611 = carto_pal(n = 12, "Bold")[3],
           N080611 = carto_pal(n = 12, "Bold")[9],
           N061011 = carto_pal(n = 12, "Bold")[12]),
  dex = c(trt = carto_pal(n = 7, "Purp")[7], 
          untrt = carto_pal(n = 7, "OrYel")[7])
)

# Sample distance heatmaps (similarities across samples)

# Compute distances (by default, Euclidean distances)
?dist
sampleVSTdists = dist(t(assay(vsd)), method = "euclidean")
sampleDistMatrix = as.matrix(sampleVSTdists)
colors = colorRampPalette(viridisLite::magma(10))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleVSTdists,
         clustering_distance_cols = sampleVSTdists,
         col = colors,
         main = "Sample Euclidean distance heatmap (VST)",
         annotation_col = df,
         annotation_colors = ann_colors)

# Heatmap of sample Poisson distances (more appropriate for count data)
poisd = PoissonDistance(t(counts(dds)))
samplePoisDistMatrix = as.matrix(poisd$dd)
dimnames(samplePoisDistMatrix) = dimnames(sampleDistMatrix)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         main = "Sample Poisson distance heatmap (counts)",
         annotation_col = df,
         annotation_colors = ann_colors)

# Principal Component Analysis
plotPCA(vsd, intgroup = c( "dex", "cell"))
pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# Generalized PCA
gpca <- glmpca(counts(dds), L=2, fam = "nb", minibatch = "none",
               optimizer = "fisher",
               ctl = list(maxIter = 10000, tol = 1e-6))
gpca.dat <- gpca$factors
gpca.dat$dex <- dds$dex
gpca.dat$cell <- dds$cell
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = dex, shape = cell)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# Multidimensional Scaling plots (MDS)
# Useful when we only have a matrix of distances instead of a matrix of data

# We are going to create a plot for the VST distances and Poisson distances
# we estimated before

# VST
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = dex, 
                shape = cell)) +
  geom_point() +
  ggtitle("MDS: VST") +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2")

# Poisson
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = dex, 
                                      shape = cell)) +
  geom_point() +
  ggtitle("MDS: Poisson distances") +
  theme_classic() +
  labs(x = "MDS1", y = "MDS2")

# Differential Gene Expression Analysis #####

# Run DGEA using DESeq() and specify the contrast with results()
dds <- DESeq(dds)

# Results will extract the comparisons of interest as specified in the contrast
# argument along with the corresponding statistics and p-values
res <- results(dds, 
               contrast=c("dex", "trt", "untrt"),
               alpha = 0.05, # default is 0.1. This refers to the adjusted p-value cutoff
               pAdjustMethod = "BH") # default: Benjamini-Hochberg

# Order the results object by ascending adjusted p-value
res = res[order(res$padj), ]

# Inspect results
res
mcols(res) # column info in res

# Get a high-level summary of the results
summary(res)

# Currently, the results have Ensembl feature names
# We can annotate them with more readable gene symbols

# The Ensembl IDs are stored in the rownames of res
# Ensembl IDs contain 15 characters plus a dot "." followed by another number
# which indicates the variant of that gene that is captured
head(rownames(res))

# Extract the 15 main characters for each Ensembl ID
ens.str <- substr(rownames(res), 1, 15)

# We can then use the org.Hs.eg.db package to assign gene symbols and Entrez IDs
# to Ensembl IDs

# Unfortunately a 1:1 relationship between IDs from different nomenclatures cannot
# be guaranteed, so here we assign the first match for each ID

# Add a column for gene symbols 
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first") # assign the first match

# Add a column for Entrez
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first") # assign the first match

# Export the results as an .xlsx or .csv file
write.xlsx(as.data.frame(res), "DGEA_results_trt_vs_untrt.xlsx", overwrite = TRUE)
write.csv(as.data.frame(res), "DGEA_results_trt_vs_untrt.csv")

# Plots based on DGEA #####

# Plot the normalized counts of the top differentially expressed gene (lowest adjusted p-value)
topGene <- rownames(res)[which.min(res$padj)] # first row in ordered object
topGeneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE, normalized = TRUE, transform = FALSE, pc = 0)
# Connecting cell lines
ggplot(topGeneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  geom_point(size = 3) + 
  geom_line(linewidth = 1) +
  labs(y = "Normalized counts")

ggplot(topGeneCounts, aes(x = dex, y = log2(count+1), color = cell, group = cell)) +
  geom_point(size = 3) + 
  geom_line(linewidth = 1) +
  labs(y = "log2-Normalized counts")

# Volcano plot ###
# We will use tidier results object for the volcano
volcano_input = as.data.frame(res) %>%
  dplyr::select(Gene.Symbol = symbol, Entrez.ID = entrez, everything()) %>% # arrange order and names for the first two columns
  dplyr::rename(Wald.statistic = stat, BH.adj.p.value = padj, p.value = pvalue,
                log2FC = log2FoldChange, lFC.SE = lfcSE)

# create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
keyvals.colour <- ifelse(
  volcano_input$log2FC < -1 & volcano_input$BH.adj.p.value < 0.05, 'royalblue',
  ifelse(volcano_input$log2FC > 1 & volcano_input$BH.adj.p.value < 0.05, 'red4',
         ifelse(abs(volcano_input$log2FC) < 1 & volcano_input$BH.adj.p.value < 0.05, 'pink',
                'grey')))

# keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'

volcano_input$aes = keyvals.colour

EnhancedVolcano(volcano_input,
                lab = volcano_input$Gene.Symbol,
                x = 'log2FC',
                y = 'BH.adj.p.value',
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType = "dashed",
                colCustom = keyvals.colour,
                selectLab = volcano_input[1:20, "Gene.Symbol"],
                title = "Volcano plot: contrast = trt vs. untrt",
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                typeConnectors = "closed",
                max.overlaps = Inf)

ggsave(filename = "Volcano_raw_LFCs_trt_vs_untrt.png",
       path = NULL, 
       width = 152, height = 142, device = 'png', units = "mm",
       dpi = 700)
dev.off()

# Post-DGEA heatmaps ###

# We will create two heatmaps:
# 1. Uses Poisson distances to cluster the samples based on the DEGs and also
#    estimates Poisson distances between the DEGs to cluster those too
# 2. Orders samples based on treatment and genes based on the results of DGEA

topDEGs = res[which(res$padj < 0.05),] # 4099 genes
topDEGs$Ensembl.ID = rownames(topDEGs) # Add Ensembl identifiers
# n_degs = 20 # the number of top DEGs we want to use
# topDEGs = topDEGs[1:n_degs, ]

# Create a custom palette with 20*8 = 160 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(160)[80:160],
                             rev(viridisLite::rocket(160)[80:160])))(160)

# Gene clustering with raw counts and Poisson distance
DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Ensembl.ID,])))
DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
rownames(DEGsamplePoisDistMatrix) = dds$SampleName
colnames(DEGsamplePoisDistMatrix) = NULL

# Same for genes
DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Ensembl.ID,]))
DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol # Convert Ensembl to Symbol for the plot
colnames(DEGgenesPoisDistMatrix) = NULL

# Heatmap
heatmat = counts(dds[topDEGs$Ensembl.ID,])
pheatmap(heatmat, annotation_col = df,
         annotation_colors = ann_colors,
         clustering_distance_rows = DEGpoisd_genes$dd,
         clustering_distance_cols = DEGpoisd$dd,
         col = colors2,
         cutree_cols = 2,
         cutree_rows = 2,
         main = "Poisson heatmap, cluster cut = 2",
         legend_labels = c("lower expression", "higher expression"))

# Heatmap with standardized counts

# Matrix normalized with size factors
sfmat = log2(counts(dds, normalized = TRUE) + 1)
sfmat_znorm = (sfmat - rowMeans(sfmat))/sqrt(rowVars(sfmat)) # Standardization

# Set rownames to symbols
sfmat_znorm = sfmat_znorm[topDEGs$Ensembl.ID, ]
rownames(sfmat_znorm) = topDEGs$Gene.Symbol

pheatmap(sfmat_znorm, annotation_col = df,
         annotation_colors = ann_colors,
         clustering_distance_rows = DEGpoisd_genes$dd,
         clustering_distance_cols = DEGpoisd$dd,
         col = colors2,
         cutree_cols = 2,
         cutree_rows = 2,
         main = "Poisson heatmap, cluster cut = 2",
         legend_labels = c("lower expression", "higher expression"))

rm(DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, sfmat,
   sfmat_znorm); gc()

# Pathway analysis #####

# Set download option to libcurl for clusterProfiler
R.utils::setOption("clusterProfiler.download.method","libcurl")

# Extract the log2 fold changes and gene names of the significant genes
log2fc_ensembl = volcano_input[topDEGs$Ensembl.ID, "log2FC"]
names(log2fc_ensembl) = topDEGs$Ensembl.ID

# Rank the list
log2fc_ensembl = sort(log2fc_ensembl, decreasing = TRUE)

# Perform GO gene set enrichment analysis using clusterProfiler
RNGversion("4.2.2")
set.seed(123)
gseaGO = gseGO(geneList = log2fc_ensembl,
               ont = "ALL",
               OrgDb = "org.Hs.eg.db",
               keyType = "ENSEMBL",
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               minGSSize = 3,
               maxGSSize = 800,
               seed = TRUE)

gseaGO = clusterProfiler::setReadable(gseaGO, 'org.Hs.eg.db')

# Perform GO ORA
oraGO = enrichGO(gene = names(log2fc_ensembl),
                 ont = "ALL",
                 universe = rownames(volcano_input),
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENSEMBL",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.1,
                 pAdjustMethod = "BH",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pool = TRUE)

oraGO = clusterProfiler::setReadable(oraGO, 'org.Hs.eg.db')
oraGO_sig = oraGO@result %>% dplyr::filter(p.adjust < 0.05)

# Prepare Entrez and Gene symbols too
log2fc_entrez = log2fc_ensembl
names(log2fc_entrez) = volcano_input[topDEGs$Ensembl.ID, "Entrez.ID"]
log2fc_entrez = sort(log2fc_entrez, decreasing = TRUE)
entrez_universe = volcano_input$Entrez.ID

# GSEA KEGG
RNGversion("4.2.2")
set.seed(123)
gseaKEGG = gseKEGG(geneList = log2fc_entrez,
                   organism = "hsa",
                   keyType = "ncbi-geneid",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 3,
                   maxGSSize = 800,
                   seed = TRUE)

gseaKEGG = clusterProfiler::setReadable(gseaKEGG, 'org.Hs.eg.db', keyType = "ENTREZID")

# KEGG
oraKEGG = enrichKEGG(gene = names(log2fc_entrez),
                     universe = entrez_universe,
                     organism = "hsa",
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     maxGSSize = 800)

oraKEGG = clusterProfiler::setReadable(oraKEGG, 'org.Hs.eg.db', keyType = "ENTREZID")
oraKEGG_sig = oraKEGG@result %>% dplyr::filter(p.adjust < 0.05)

# Example plots ###

# Barplot ORA-WP
barplot(oraKEGG, showCategory = 10)

# Dotplot ORA-KEGG
dotplot(oraKEGG, showCategory = 10)

# Heatplot GSEA-GO
heatplot(gseaGO, foldChange = log2fc_ensembl, showCategory = 10)

# Category netplot GSEA GO
cnetplot(gseaGO, color.params = list(foldChange = log2fc_ensembl), categorySize = "pvalue") 

# Pathview with KEGG
topKEGG = gseaKEGG[1, ]
pathway_id1 = topKEGG$ID
pathway_id2 = gseaKEGG[2, "ID"]
pathway_id3 = gseaKEGG[3, "ID"]

# Generate the pathway diagram
pv1 = pathview(gene.data = log2fc_entrez, pathway.id = pathway_id1, species = "hsa", 
               new.signature = FALSE, kegg.native = TRUE, out.suffix = "pv1")
pv2 = pathview(gene.data = log2fc_entrez, pathway.id = pathway_id2, species = "hsa", 
               new.signature = FALSE, kegg.native = TRUE, out.suffix = "pv2")
pv3 = pathview(gene.data = log2fc_entrez, pathway.id = pathway_id3, species = "hsa", 
               new.signature = FALSE, kegg.native = TRUE, out.suffix = "pv3")

# Tree plot: ORA KEGG
pw = enrichplot::pairwise_termsim(oraKEGG)
enrichplot::treeplot(pw, cluster.params = list(method = "average"))

# Session info #####
sessionInfo()
