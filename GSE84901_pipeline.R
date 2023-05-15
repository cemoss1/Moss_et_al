# loading in GEO dataset

library(GEOquery)
gse <- getGEO("GSE84901", GSEMatrix = TRUE, getGPL = FALSE) #series matrix file will be read in

#check how many platforms were used in the dataset and assign

length(gse)
gse <- gse[[1]]

#check expression, sample information and annotation have been  downloaded

summary(exprs(gse)) #prints summary of expression data
pData(gse) #prints sample information
fData(gse) #prints gene annotation
exprs(gse) #prints expression data

#remove unwanted samples in the dataset
##filter out samples marked with "X" e.g. the treated samples

gsms <- "000000000000XXXXXXXXXXX"
sml <- strsplit(gsms, split = "")[[1]]
sel <- which(sml !="X")
sml <- sml[sel]
gse <- gse[ ,sel]

#method assumes use of log2 scale with values ranging from 0-16 so check distribution from summary(exprs(gse))
##if not in range
#exprs(gse) <- log2(exprs(gse))

#boxplot to check distribution

boxplot(exprs(gse), outline = FALSE)

#neaten up and add labels

par(mar = c(7,5,3,3), cex.main = 1)
boxplot(exprs(gse), outline=FALSE, las =2, ylab = 'Log2 normalised expression', title('Boxplot showing normalised distribution of GSE84901 untreated samples'))

#inspect clinical variables

library(dplyr)
sampleInfo <- pData(gse)
sampleInfo

#select columns for analysis

sampleInfo <- dplyr::select(sampleInfo, `treatment:ch1`,`young or aged:ch1`)

##optionally rename column names for convenience

sampleInfo <- dplyr::rename(sampleInfo,treatment = `treatment:ch1`, Age = `young or aged:ch1`)
sampleInfo

#create correlation matrix and visualise in a heatmap

library(pheatmap)
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix, annotation_col = sampleInfo, main = "Correlation Matrix for GSE84901 Untreated Samples")

#use ggplot2 and ggrepel to produce PCA plot

library(ggplot2)
library(ggrepel)

##make sure to transpose the expression matrix

pca <- prcomp(t(exprs(gse)))

#join the Principle Components to the sample info

cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=Age,label=paste(Age))) + geom_point() + geom_text_repel()

#alternate method
library(ggfortify)
autoplot(pca, data = sampleInfo, colour = 'Age') +ggtitle(label = "Principle Component Analysis of GSE84901 Untreated Samples") + geom_text_repel(aes(col=Age, label = rownames(sampleInfo)))

#differential expression
##allocate samples in dataset to sample groups of interest

library(limma)
design <- model.matrix(~0+sampleInfo$Age)
design

##rename column names

colnames(design) <- c("Aged","Young")

##cut-off level of median expression level

cutoff <- median(exprs(gse))

#find genes expressed in each sample

is_expressed <- exprs(gse) > cutoff

#identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

#check how many genes are removed/retained

table(keep)

#subset to just expressed genes

gse <- gse[keep,]

#lmFit fits the model to the data

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

#define contrasts interested in for differential analysis

contrasts <- makeContrasts(Aged - Young, levels = design)
fit2 <- contrasts.fit(fit, contrasts)

#apply Bayes step to get differential expression

fit2 <- eBayes(fit2)
topTable(fit2)

#see how many genes are differentially expressed
decideTests(fit2)
table(decideTests(fit2))

#calculate relative array weights to deal with outliers

aw <- arrayWeights(exprs(gse),design)
aw
fit <- lmFit(exprs(gse), design, weights = aw)
contrasts <- makeContrasts(Aged - Young, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

#further processing differential expression

anno <- fData(gse)
anno
anno <- dplyr::select(anno,ID,GB_ACC,GENE_SYMBOL,ENSEMBL_ID,CHROMOSOMAL_LOCATION,CYTOBAND)
fit2$genes <- anno
topTable(fit2)

#create data frame for volcano plot, label plots and genes of interest

full_results <- topTable(fit2, number = Inf)
full_results <- tibble::rownames_to_column(full_results,"NAME")
library(ggplot2)
library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>%
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff) %>%
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, GENE_SYMBOL, "")) %>%
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() +geom_text_repel(max.overlaps = 12, col="black")

#heatmaps of DE genes

library(pheatmap)
library(viridis)

#differentially expressed genes

topN <- 20000
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>%
  filter(Rank < topN) %>%
  pull(ID)

#extract corresponding gene symbols

gene_names <- mutate(full_results, Rank = 1:n()) %>%
  filter(Rank < topN) %>%
  pull(GENE_SYMBOL)

#get the rows corresponding to ids_of_interest and all columns

gene_matrix <- exprs(gse)[ids_of_interest,]

#make heatmap simply

pheatmap(gene_matrix, labels_row = gene_names)

##scale rows to highlight differences and add all features

pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, color = plasma(30), fontsize = 8)


#user defined genes of interest
##cell cycle, initiation of mytosis

my_genes <- c("Akt1", "Mnat1", "Usf1", "Mycbp", "Lmnb1", "Hist1h1a", "Cdc25c", "Akt1s1", "Cdk1", "Cdk7", "Ccnh", "Lmnb2", "Cdk10", "Hist1h1b", "Weel", "P1k1", "Myc", "Ccnb1", "Akt2", "Foxm1", "Pkmyt1", "Kif11", "Ccnb2")
ids_of_interest <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(ID)
gene_names <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(GENE_SYMBOL)
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Cell Cycle: Initiation of Mitosis")

#tidied up version

pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Cell Cycle: Initiation of Mitosis", color = plasma(30), border_color = NA, show_colnames = FALSE, fontsize = 14, cluster_cols = FALSE)

#cell cycle metaphase genes

my_genes <- c("Dync1h1", "Mad111", "Mad211b", "Pmf1", "Aurka", "Incenp", "Cenpb", "Spc24", "Plk1", "Birc5", "Kntc1", "Bub1b", "Nuf2", "Bub1", "Cenph", "Aurkb", "Zwilch", "Dsn1", "Cenpf", "Mad211", "Nsl1", "Cdc20", "Cenpa")
ids_of_interest <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(ID)
gene_names <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(GENE_SYMBOL)
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Cell Cycle: Metaphase Genes")

#tidied version

pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Cell Cycle: Metaphase Genes", color = plasma(30), border_color = NA, show_colnames = FALSE, fontsize = 14, cluster_cols = FALSE)

#housekeeping genes

my_genes <- c("Gapdh", "B2m", "Gusb", "Hmbs", "Hnrnpab", "Hprt", "Mau2", "Ppia", "Rpl13a", "Stx5a", "Pum1", "Tbp", "Oaz1", "Tmem199", "Actb", "Rpl32", "Gnb2l1", "Abl1")

ids_of_interest <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(ID)
gene_names <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(GENE_SYMBOL)
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Housekeeping Genes")

#tidied version

pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Housekeeping Genes", color = plasma(30), border_color = NA, show_colnames = FALSE, fontsize = 14, cluster_cols = FALSE)

#conserved genes

my_genes <- c("Cpsf3", "Srp72", "Mapre1", "Api5", "Stt3b", "Caprin1", "Sec31a", "Canx", "Ccni", "Vcp", "Arl6ip1", "Stt3a", "Sec24c", "Psmd1", "Copb1", "Yeats4", "Etf1", "Nrd1", "Hnrnpu", "Bzw1", "Ssr3", "Hnrnpc", "Ddx1", "Psmc6", "Psmb1", "Abcf1", "Tm9sf3", "Eif3d", "Psmc2", "St13", "Nars")


ids_of_interest <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(ID)
gene_names <- filter(full_results,GENE_SYMBOL %in% my_genes) %>%
  pull(GENE_SYMBOL)
gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Conserved Genes")

#tidied version

pheatmap(gene_matrix, labels_row = gene_names, scale = "row", annotation_col = sampleInfo, main = "Conserved Genes", color = plasma(30), border_color = NA, show_colnames = FALSE, fontsize = 14, cluster_cols = FALSE)

cons <- filter(full_results, GENE_SYMBOL %in% my_genes)
library(readr)
cons%>%
  write_csv(path="conserved_de_results.csv")
  
  #define significant genes for ontology analysis...335 genes in output

sig <- filter(full_results, adj.P.Val < 0.05, abs(logFC) > 1.5)
library(readr)
sig%>%
  write_csv(path="filtered_de_results.csv")

#prepare gene lists for enrichR...213 upregulated and 122 downregulated

pos <- filter(sig, logFC > 1.5)
pos <- dplyr::select(pos, GENE_SYMBOL)
neg <- filter(sig, logFC < 1.5)
neg <- dplyr::select(neg, GENE_SYMBOL)

#analysis with enrichR

library(enrichR)
dbs <- c("KEGG_2019_Mouse", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TargetScan_microRNA")
posenriched <- enrichr(pos[,1],dbs)
negenriched <- enrichr(neg[,1],dbs)

#plot tables with enrichR
##upregulated

plotEnrich(posenriched[[1]], showTerms = 20, title = "Upregulated Gene Pathways (KEGG Mouse 2019)")
plotEnrich(posenriched[[2]], showTerms = 20, title = "Encode and ChEA Consesus TFs from Upregulated Genes")
plotEnrich(posenriched[[3]], showTerms = 20, title = "TargetScan microRNA for Upregulated Target Genes")

##downregulated

plotEnrich(negenriched[[1]], showTerms = 20, title = "Downregulated Gene Pathways (KEGG Mouse 2019)")
plotEnrich(negenriched[[2]], showTerms = 20, title = "Encode and ChEA Consesus TFs from Downregulated Genes")
plotEnrich(negenriched[[3]], showTerms = 20, title = "TargetScan microRNA for Downregulated Target Genes")

#pull out TFs in stringent gene lists

upregTFs <- posenriched[[2]]
upregTFs <- upregTFs[c(4,8,12,17,25),]
downregTFs <- negenriched[[2]]
downregTFs <- downregTFs[c(2,3,11,15,20,25),]

#join up and down reg TF lists

TFlist <- rbind(upregTFs, downregTFs)
TFlist %>%
  write_csv(path = "filtered_de_TFs.csv")
