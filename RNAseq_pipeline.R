## specific to M0 MYC knockdown vs M0 but same pipeline for M0 USF1 knockdown vs M0 and M1 USF1 knockdown vs M1
library(tximportData) #quantifies transcript abundance
library(readr) #reads deliminated files
library(tximport) #transcript length estimates
library (rjson) #json data interchange format
library (limma) #differential expression analysis
library (edgeR) #emperical analysis of digital gene expression data
library (gplots) #for plotting data
library(GenomicFeatures) #finctionality for transcript database analysis

#{r tx2gene} #transcript to gene mapping

TxDb <- makeTxDbFromGFF(file ="~/quantification/quantification.dir/geneset_all1.gtf")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDb, k, "GENEID", "TXNAME")

#{r tximport}

library(tximport)
setwd("~/RNA sequencing/Analysis/RNA seq analysis in R/M0vsM0MYCKD")
file = list.files(path ="~/RNA sequencing/Analysis/RNA seq analysis in R/M0vsM0MYCKD", pattern=".sf$")
Files=unlist(file)
name = gsub(".sf", "", Files)
names(Files) <- name
txi <- tximport(Files, type = "salmon", tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM")

df = txi$counts
df_tmp = txi$abundance
colnames(df) = name
txi1 <- tximport(Files, type = "salmon", txOut = TRUE)
KD = txi1$counts
df_tmp = txi1$abundance
colnames(KD) = name


#{r CPM_counts}
kd_group = read.csv("~/RNA sequencing/Analysis/RNA seq analysis in R/M0vsM0MYCKD/M0vsM0MYCKD_Group_KD.csv")

remove_zero = KD[rowSums(KD)>0, ] 
count = DGEList (count = remove_zero) ##
TMM_count = calcNormFactors(count, method = "TMM")

remove_zero = df[rowSums(df)>0, ] 
count = DGEList (count = remove_zero) ##
TMM_count = calcNormFactors(count, method = "TMM")
CPM_Count_kd = cpm(TMM_count, normalized.lib.sizes=TRUE, log=T, prior.count=1)
write.csv(CPM_Count_kd, "CPM_counts.csv")

#PCA plot of individual samples{r,fig.height=2, fig.width=2.5}

library(factoextra)
pca_data = t(CPM_Count_kd)
wdbc.pr <- prcomp(pca_data, center = TRUE, scale = TRUE)
fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 5, 
             fill.ind = kd_group$Group, 
             col.ind = "black", 
             palette = c("#000000","#00B0F0"), 
             addEllipses = FALSE,
             label = "ind",
             labelsize = 3,
             col.var = "black",
             repel = TRUE,
             legend.title = "Conditions") +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),legend.title  = element_text(size = 12),legend.text = element_text(size = 12))

#differential expression analysis

#{r DESEQ2}

library (DESeq2)
rownames(kd_group) = kd_group$Sample
kd_dds <- DESeqDataSetFromTximport(txi, kd_group, design=~Donor+Group)

d <- plotCounts(kd_dds, gene = "ENSG00000136997", intgroup = "Group", returnData = TRUE) #MYC
library("ggplot2")
legend_title <- "Donor"
ggplot(d, aes(x=Group, y=count)) + 
  geom_line(aes(colour=kd_group$Donor, group=kd_group$Donor)) + 
  geom_point(aes(colour=kd_group$Donor),               
             size=4)+
  xlab("Condition")+ylab("Count")


kd_data <- estimateSizeFactors(kd_dds, type = "poscounts", locfunc=genefilter::shorth)
kd_data <- estimateDispersions(kd_data)
kd_data <- nbinomWaldTest(kd_data)
results = results(kd_data)
Signif_kd <- results(kd_data , alpha = 0.9999999999)
summary(Signif_kd)
res_kd = Signif_kd[which(Signif_kd$padj < 0.9999999999 & abs(Signif_kd$log2FoldChange) > 1),]
res_kd_cpm = CPM_Count_kd[rownames(CPM_Count_kd)%in%rownames(res_kd),]
downregulated_RNAs = res_kd[res_kd$log2FoldChange <= 0 ,]
upregulated_RNAs = res_kd[res_kd$log2FoldChange > 0 ,]
write.csv (downregulated_RNAs,"Downregulated_mRNA.csv")
write.csv (upregulated_RNAs,"Upregulated_mRNA.csv")

Signif_all <- results(kd_data, alpha = 0.99999999)
summary(Signif_all)
res_all <- Signif_all[which(Signif_all$padj<0.9999999999 & abs(Signif_all$log2FoldChange)>0),]

library (biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_wt_all = as.character(unlist(rownames(res_all)))
res_wt_all = getBM(attributes= c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = ensembl_wt_all, 
                   mart = ensembl)
all_ign = merge(as.data.frame(res_all),res_wt_all,by.x = "row.names",by.y = "ensembl_gene_id")
rownames(all_ign) = make.names(all_ign$ensembl_transcript_id,unique = T)

library(EnhancedVolcano)
#labels
EnhancedVolcano(all_ign, x="log2FoldChange", y="pvalue", lab=all_ign$hgnc_symbol, 
                title = "Volcano plot", subtitle = bquote(italic("M0 vs M0_MYC_KD")), 
                pCutoff = 0.05,drawConnectors = FALSE, labSize = 1.5, axisLabSize = 8, 
                colAlpha = 1, titleLabSize = 8, subtitleLabSize = 6, captionLabSize = 6, 
                pointSize = 1, legendLabSize = 6, legendIconSize = 3, 
                legendPosition = 'bottom', legendDropLevels = TRUE, borderWidth = 0, 
                col = c("grey30", "forestgreen","red2", "royalblue"), 
                xlim = c(min(all_ign[["log2FoldChange"]], na.rm = TRUE) - 0.5, max(all_ign[["log2FoldChange"]], na.rm = TRUE) +0.5), 
                ylim = c(0, max(-log10(all_ign[["pvalue"]]), na.rm = TRUE) + 0.5))

#no labels
EnhancedVolcano(all_ign, x="log2FoldChange", y="pvalue", lab=NA, title = "Volcano plot", subtitle = bquote(italic("M0 vs M0_MYC_KD")), pCutoff = 0.05,drawConnectors = FALSE, axisLabSize = 8, 
                colAlpha = 1, titleLabSize = 8, subtitleLabSize = 6, captionLabSize = 6, 
                pointSize = 1.5, legendLabSize = 6, legendIconSize = 3, 
                legendPosition = 'bottom', legendDropLevels = TRUE, borderWidth = 0, 
                col = c("grey30", "forestgreen","red2", "royalblue"), 
                xlim = c(min(all_ign[["log2FoldChange"]], na.rm = TRUE) - 0.5, max(all_ign[["log2FoldChange"]], na.rm = TRUE) +0.5), 
                ylim = c(0, max(-log10(all_ign[["pvalue"]]), na.rm = TRUE) + 0.5))

#{r}
library(factoextra)
pca_data = t(as.matrix(res_kd_cpm))
wdbc.pr <- prcomp(pca_data, center = TRUE, scale = TRUE)
fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = kd_group$Group, 
             col.ind = "black", 
             palette = c("#000000","#00B0F0"), 
             addEllipses = TRUE,
             #label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Conditions") +
  theme(text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),legend.title  = element_text(size = 5),legend.text = element_text(size = 5))


# Add transcript ID and gene symbol on signicicant list of RNAseq from kd
#{r add annnotations}

library (biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensembl_wt = as.character(unlist(rownames(res_kd)))
res_wt = getBM(attributes= c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol'), 
               filters = 'ensembl_gene_id', 
               values = ensembl_wt, 
               mart = ensembl)
signif_ign = merge(as.data.frame(res_kd),res_wt,by.x = "row.names",by.y = "ensembl_gene_id")
rownames(signif_ign) = make.names(signif_ign$ensembl_transcript_id,unique = T)


downregulated_RNAs_kd = merge(as.data.frame(downregulated_RNAs),res_wt,by.x = "row.names",by.y = "ensembl_gene_id")
upregulated_RNAs_kd = merge (as.data.frame(upregulated_RNAs),res_wt,by.x = "row.names",by.y = "ensembl_gene_id")
write.csv (downregulated_RNAs_kd,"Downregulated_mRNA.csv")
write.csv (upregulated_RNAs_kd,"Upregulated_mRNA.csv")

#DE heatmap
library(pheatmap)
library(RColorBrewer)
pheatmap(res_kd_cpm, color = colorRampPalette(rev(brewer.pal(n = 6, name ="RdYlBu")))(100), scale = "row", treeheight_col = 35, show_rownames = FALSE, cluster_rows = TRUE, treeheight_row = 0, angle_col = 90)


heatmap_kd <- Signif_kd[which(Signif_kd$padj < 0.999999 & abs(Signif_kd$log2FoldChange) > 1),]
heatmap_kd_cpm = CPM_Count_kd[rownames(CPM_Count_kd)%in%rownames(heatmap_kd),]
pheatmap(heatmap_kd_cpm, color = colorRampPalette(rev(brewer.pal(n = 6, name ="RdYlBu")))(100), scale = "row", treeheight_col = 35, show_rownames = FALSE, cluster_rows = TRUE, treeheight_row = 0, angle_col = 90)

#gene ontology analysis

upreg <- dplyr::select(upregulated_RNAs_kd, hgnc_symbol)
downreg <- dplyr::select(downregulated_RNAs_kd, hgnc_symbol)
upreg <- dplyr::distinct(upreg)
downreg <- dplyr::distinct(downreg)
library(enrichR)
dbs <- c("GO_Biological_Process_2021","GO_Molecular_Function_2021","KEGG_2021_Human", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TargetScan_microRNA")
upregenriched <- enrichr(upreg[,1],dbs)
downregenriched <- enrichr(downreg[,1],dbs)

upprocesses <- upregenriched[[1]]
downprocesses <- downregenriched[[1]]
upfunction <- upregenriched[[2]]
downfunction <- downregenriched[[2]]
upkegg <- upregenriched[[3]]
downkegg <- downregenriched[[3]]
upTFs <- upregenriched[[4]]
downTFs <- downregenriched[[4]]
upMiRs <- upregenriched[[5]]
downMiRs <- downregenriched[[5]]

write.csv(upprocesses, file = "Upregulated_biological_processes.csv")
write.csv(downprocesses, file = "Downregulated_biological_processes.csv")
write.csv(upfunction, file = "Upregulated_molecular_function.csv")
write.csv(downfunction, file = "Downregulated_molecular_function.csv")
write.csv(upkegg, file = "Upregulated_KEGG_pathways.csv")
write.csv(downkegg, file = "Downregulated_KEGG_pathways.csv")
write.csv(upTFs, file = "Upregulated_TFs.csv")
write.csv(downTFs, file = "Downregulated_TFs.csv")
write.csv(upMiRs, file = "Upregulated_MiRs.csv")
write.csv(downMiRs, file = "Downregulated_MiRs.csv")

plotEnrich(upprocesses, showTerms = 20)
plotEnrich(downprocesses, showTerms = 20)
plotEnrich(upfunction, showTerms = 20)
plotEnrich(downfunction, showTerms = 20)
plotEnrich(upkegg, showTerms = 20)
plotEnrich(downkegg, showTerms = 20)
plotEnrich(upTFs, showTerms = 20)
plotEnrich(downTFs, showTerms = 20)
plotEnrich(upMiRs, showTerms = 20)
plotEnrich(downMiRs, showTerms = 20)

