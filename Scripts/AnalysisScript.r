###############################################
# BINF 6110 Assignment 2
#
# Bulk Transcriptomics Differential Expression Analysis
#
# 2026-03-01
#
################################################


# Load packages

library(DESeq2)
library(tximport)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(Biostrings)
library(ashr)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(cowplot)
library(ggrepel)
library(dplyr)

# create metadata table

sample_info <- data.frame(
  sample_id = c("SRR10551665", "SRR10551664", "SRR10551663", "SRR10551662" ,"SRR10551661" ,"SRR10551660", "SRR10551659", "SRR10551658", "SRR10551657"),
  stage = factor(c("Early", "Early", "Early", "Thin", "Thin", "Thin", "Mature" ,"Mature", "Mature"),
                 levels = c("Early","Thin","Mature")),
  row.names = "sample_id"
)


#tx2gene from gtf

gtf_file <- "C:/Users/mahno/Onedrive/Desktop/BINF6110/GCF_000146045.2_R64_genomic.gtf.gz"

txdb <- makeTxDbFromGFF(gtf_file)

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = k,
  columns = "GENEID",
  keytype = "TXNAME"
)


#import salmon quants
quant_dir <- "C:/Users/mahno/Onedrive/Desktop/BINF6110/quants"

files <- file.path(quant_dir, rownames(sample_info), "quant.sf")
names(files) <- rownames(sample_info)
length(files)


txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#DESeq2 with early as reference
dds <- DESeqDataSetFromTximport(txi, sample_info, ~stage)
dds <- DESeq(dds)
resultsNames(dds)

# no thin vs mature so relevel with thin as reference for thin vs mature
dds_tvm <- dds
dds_tvm$stage <- relevel(dds_tvm$stage, ref = "Thin")
dds_tvm <- DESeq(dds_tvm)

#  pairwise contrasts
res_tve   <- results(dds, contrast = c("stage", "Thin", "Early"), alpha = 0.05)
res_mve <- results(dds, contrast = c("stage", "Mature", "Early"), alpha = 0.05)

res_mvt <- results(dds_tvm, contrast = c("stage", "Mature", "Thin"), alpha = 0.05)



#LCF shrinkage to reduce noise

LFC_tve <- lfcShrink(dds,
                     coef = "stage_Thin_vs_Early",
                     type = "apeglm"
)

LFC_mve <- lfcShrink(dds,
                     coef = "stage_Mature_vs_Early",
                     type = "apeglm"
)

LFC_mvt <- lfcShrink(dds_tvm,
                     coef = "stage_Mature_vs_Thin",
                     type = "apeglm"
)


# Save results to CSV
write.csv(as.data.frame(LFC_tve),   "C:/Users/mahno/Onedrive/Desktop/BINF6110/results/thin_vs_early.csv")
write.csv(as.data.frame(LFC_mve), "C:/Users/mahno/Onedrive/Desktop/BINF6110/results/mature_vs_early.csv")
write.csv(as.data.frame(LFC_mvt),  "C:/Users/mahno/Onedrive/Desktop/BINF6110/results/mature_vs_thin.csv")


# DEG summary
summary(res_tve)
summary(res_mve)
summary(res_mvt)

# filter for LFC > 1
count_DEGs <- function(res, padj_cutoff=0.05, log2FC_cutoff=1){
  
  DEGs <- res[which(res$padj < padj_cutoff & abs(res$log2FoldChange) > log2FC_cutoff), ]
  
  n_up <- nrow(DEGs[DEGs$log2FoldChange > log2FC_cutoff, ])
  
  n_down <- nrow(DEGs[DEGs$log2FoldChange < -log2FC_cutoff, ])
  
  n_total <- n_up + n_down
  
  return(list(total=n_total, up=n_up, down=n_down))
}


count_DEGs(res_tve)
count_DEGs(res_mve)
count_DEGs(res_mvt)




# PCA for overall data structure visualization
vsd <- vst(dds) 
pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x=PC1, y=PC2, color=stage)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA of Yeast Biofilm Stages") +
  theme_bw() +
  coord_fixed()


ggsave("C:/Users/mahno/Onedrive/Desktop/BINF6110/results/pca_plot.pdf")


# Volcano plot function for each comparison

Volcano <- function(res_v, title) {
  df_v <- as.data.frame(res_v)
  df_v$gene <- rownames(df_v)
  df_v$sig <- ifelse(df_v$padj < 0.05 & abs(df_v$log2FoldChange) > 1,
                     ifelse(df_v$log2FoldChange > 0, "Up", "Down"), "Not Sig")
  df_v <- na.omit(df_v)
  
  
  # labels for top 5 up and down by padj
  top_labels <- rbind(df_v %>% filter(sig == "Up") %>% slice_min(padj, n = 5), df_v %>% filter(sig == "Down") %>% slice_min(padj, n = 5))
  
  
  ggplot(df_v, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.3, size = 2) +
    scale_color_manual(values = c("Down" = "navy", "Not Sig" = "gray", "Up" = "maroon")) +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value") +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_text_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    theme_bw()
}


Volcano(LFC_tve, "Thin vs Early Biofilm")
Volcano(LFC_mve, "Mature vs Early Biofilm")
Volcano(LFC_mvt, "Mature vs Thin Biofilm")



# Heatmap of top genes by minimum adjusted p-value across all three comparisons

get_top_genes <- function(res1, res2, res3, n = 30) {
  df1 <- as.data.frame(res1) %>% na.omit()
  df2 <- as.data.frame(res2) %>% na.omit()
  df3 <- as.data.frame(res3) %>% na.omit()
  
  all_genes <- union(union(rownames(df1), rownames(df2)), rownames(df3))
  
  padj_matrix <- data.frame(
    r1 = df1[all_genes, "padj"],
    r2 = df2[all_genes, "padj"],
    r3 = df3[all_genes, "padj"]
  )
  rownames(padj_matrix) <- all_genes
  
  min_padj <- apply(padj_matrix, 1, function(x) min(x, na.rm = TRUE))
  top_genes <- names(sort(min_padj))[1:n]
  
  return(top_genes)
}

# Builds and plots the heatmap
plot_heatmap <- function(dds, top_genes, title = "Top Differentially Expressed Genes") {
  vsd <- vst(dds)
  mat <- assay(vsd)[top_genes, ]
  
  ann <- as.data.frame(colData(dds)[, "stage", drop = FALSE])
  ann_colors <- list(stage = c(Early = "limegreen", Thin = "orange", Mature = "darkred"))
  
  pheatmap(mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = ann,
           annotation_colors = ann_colors,
           annotation_names_col = FALSE,
           show_colnames = FALSE,
           main = title)}


top_genes <- get_top_genes(LFC_tve, LFC_mve, LFC_mvt, n = 30)
plot_heatmap(dds, top_genes, "30 most significant genes by adj. p-value")



# GO ORA 
# GO enrichment separate for up and downregulated genes
GO <- function(res_df, ontology = "BP", comparison_name = "") {
  
  df <- as.data.frame(res_df) %>% mutate(ORF = rownames(.))
  
  up   <- df %>% filter(padj < 0.05 & log2FoldChange >  1) %>% pull(ORF) %>% na.omit() %>% unique()
  down <- df %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(ORF) %>% na.omit() %>% unique()
  all  <- df %>% pull(ORF) %>% na.omit() %>% unique()
  
  res <- compareCluster(
    geneClusters  = list(Upregulated = up, Downregulated = down),
    fun           = "enrichGO",
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = ontology,
    universe      = all,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2)
  
  res <- clusterProfiler::simplify(res, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  p <- dotplot(res, showCategory = 10) +
    ggtitle(paste("GO", ontology, "ORA:", comparison_name)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 8))
  
  list(result = res, plot = p)}




# KEGG ORA
#KEGG pathway enrichment separate for up and downregulated genes
KEGG_ora <- function(res_df, comparison_name = "") {
  
  df <- as.data.frame(res_df) %>% mutate(ORF = rownames(.))
  
  up   <- df %>% filter(padj < 0.05 & log2FoldChange >  1) %>% pull(ORF) %>% na.omit() %>% unique()
  down <- df %>% filter(padj < 0.05 & log2FoldChange < -1) %>% pull(ORF) %>% na.omit() %>% unique()
  
  res <- compareCluster(
    geneClusters  = list(Upregulated = up, Downregulated = down),
    fun           = "enrichKEGG",
    organism      = "sce",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2)
  
  p <- dotplot(res, showCategory = 10) +
    ggtitle(paste("KEGG ORA:", comparison_name)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 8))
  
  list(result = res, plot = p)}

# Run and plot functions

go_tve  <- GO(LFC_tve, comparison_name = "Thin vs Early")
go_mve  <- GO(LFC_mve, comparison_name = "Mature vs Early")
go_mvt  <- GO(LFC_mvt, comparison_name = "Mature vs Thin")

plot_grid(go_tve$plot, go_mve$plot, go_mvt$plot, ncol = 3)


kegg_tve <- KEGG_ora(LFC_tve, comparison_name = "Thin vs Early")
kegg_mve <- KEGG_ora(LFC_mve, comparison_name = "Mature vs Early")
kegg_mvt <- KEGG_ora(LFC_mvt, comparison_name = "Mature vs Thin")

plot_grid(kegg_tve$plot, kegg_mve$plot, kegg_mvt$plot, ncol = 3)


