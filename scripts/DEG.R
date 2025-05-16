##___________________ Installation
install.packages("easypackages")

BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
pkgs <- c("EnhancedVolcano", "org.Hs.eg.db", "clusterProfiler", "ReactomePA", "DESeq2", )
if (!requireNamespace(pkgs, quietly = TRUE))
  BiocManager::install(pkgs)
#___________________
library(ggplot2)
library(dplyr)
library(DESeq2)
library(ReactomePA)
library(enrichplot)
library(RSQLite)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

pkgs1 <- c("ggplot2", "dplyr")
install.packages(pkgs1)
#___________________

##---- Analysis start
setwd("/media/bandana/DATA/Bandana/trans/DEGs/Deseq")
sample_ids <- c("SRR26707520", "SRR26707524", "SRR26707526", "SRR26707527", "SRR32259398", "SRR32259400")
count_mat <- NULL
for (sample_id in sample_ids) {
  feature_file <- paste0(sample_id, "_featureCounts.txt")
  feature_data <- read.table(feature_file, header = TRUE, row.names = 1)
  feature_data <- cbind(Genes = rownames(feature_data), feature_data)
  rownames(feature_data) <- NULL
  filtered_data <- feature_data[feature_data[[ncol(feature_data)]] > 10, ]
  final_data <- filtered_data[, c(1, ncol(filtered_data))]
  colnames(final_data)[2] <- sample_id
  
  if (is.null(count_mat)) {
    count_mat <- final_data
  } else {
    count_mat <- merge(count_mat, final_data, by = "Genes", all = TRUE)
  }
}
count_mat[is.na(count_mat)] <- 0
head(count_mat)
rownames(count_mat) <- sub("\\..*", "", rownames(count_mat))
# count_mat <- count_mat[rowSums(count_mat >= 10) >= 2, ]

## Sample info
sample_info <- data.frame(
  sample = colnames(count_mat)[-1],  # remove "Genes" column
  condition = c("cancer", "cancer", "cancer", "cancer", "control", "control")
)
rownames(sample_info) <- sample_info$sample

## preprocessing count data
counts <- count_mat
rownames(counts) <- counts$Genes
counts <- counts[, -1]  # remove "Genes" column
counts <- counts[, rownames(sample_info)]   ## Ensure the columns match the sample metadata

## Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
# resOrd <- res[order(res$pvalue), ]
# head(resOrd)

# Filtration based on threshold
res_filt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]

res_filt <- res_filtered[!is.na(res_filtered$padj), ]

res_rank <- res_filtered[order(res_filtered$padj), ]

write.csv(as.data.frame(res), "deseq2_all_results.csv")
write.csv(as.data.frame(res_rank), "deseq2_significant_DEGs.csv")

###----- Converting ensemble IDs into Gene Symbols
res_rank$symbol <- mapIds(org.Hs.eg.db,
                          keys = rownames(res_ranked),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Save again with gene symbols
write.csv(as.data.frame(res_ranked), "significant_DEGs_with_symbols.csv")
##---------------------------------------- Vizualization ------------------------------------
# MA plot
plotMA(res, ylim = c(-5, 5))
####------------------------------------------------- Volcano plot------------------------------------------------------------------
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# Add a column for significance
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1] <- "Significant"
res_df$Significance[res_df$padj < 0.05 & abs(res_df$log2FoldChange) < 1] <- "Low LFC"
res_df$Significance[res_df$padj >= 0.05 & abs(res_df$log2FoldChange) >= 1] <- "High LFC, not sig"

# Create volcano plot
viol <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey",
                                "Low LFC" = "orange",
                                "High LFC, not sig" = "blue",
                                "Significant" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot: DESeq2 Results",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Legend")
ggsave("Violin_plot_DEG.png", plot = viol, width = 8, height = 6, dpi = 600)
####---------------------------------------- Functional Enrichment analysis ---------------------------------------------------------------
## Gene Ontology
gene_symb <- res_rank$symbol
gene_symb <- gene_symbols[!is.na(gene_symbols)]
gene_symb <- unique(gene_symbols)
gene_symb <- res_ranked$symbol[!is.na(res_rank$symbol)]
go <- enrichGO(gene         = gene_symb,
               OrgDb        = org.Hs.eg.db,
               keyType      = "SYMBOL",
               ont          = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05)
d_go <- data.frame(go)




## KEGG Pathway Enrichment
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)
kegg <- enrichKEGG(gene         = entrez_ids,
                   organism     = 'hsa',   # can be changed based on organism
                   pvalueCutoff = 0.05)

## Reactome Pathway Enrichment
react_p <- enrichPathway(gene         = entrez_ids,
                         organism     = "human",
                         pvalueCutoff = 0.05,
                         readable     = TRUE)

#### Save all results
write.csv(data.frame(go), "Gene_ontology.csv")
write.csv(data.frame(kegg), "KEGG_pathway_Enrivhment.csv")
write.csv(data.frame(react_p), "Reactome_pathway_Enrivhment.csv")

##--------- visualization

# barplot
b_BP <- barplot(go, showCategory = 20, title = "GO Biological Process")   #GO(BP)
ggsave("dotplot_GO Biological Process_600dpi.png", plot = b_BP, width = 8, height = 6, dpi = 600)

b_kegg <- barplot(kegg, showCategory = 20, title = "KEGG Pathway Enrichment")
ggsave("dotplot_KEGG Pathway Enrichment_600dpi.png", plot = b_kegg, width = 8, height = 6, dpi = 600)


b_react_p <- barplot(react_p, showCategory = 20, title = "Reactome Pathway Enrichment")
ggsave("dotplot_Reactome Pathway Enrichment_600dpi.png", plot = b_react_p, width = 8, height = 6, dpi = 600)

# dotplot
d_BP <- dotplot(go, showCategory = 20, title = "GO Biological Process")
ggsave("dotplot_GO Biological Process__600dpi.png", plot = d_BP, width = 8, height = 6, dpi = 600)

d_kegg <- dotplot(kegg, showCategory = 20, title = "KEGG Pathway Enrichment")
ggsave("dotplot_KEGG Pathway Enrichment_600dpi.png", plot = d_kegg, width = 8, height = 6, dpi = 600)

d_react_p <- dotplot(react_p, showCategory = 20, title = "Reactome Pathway Enrichment")
ggsave("dotplot_Reactome Pathway Enrichment_600dpi.png", plot = d_react_p, width = 8, height = 6, dpi = 600)

#cnetplot
c_BP <- cnetplot(
  go,
  layout = igraph::layout_with_kk,
  showCategory = 5,
  color_category = "#E5C494",
  size_category = 1,
  color_item = "#B3B3B3",
  size_item = 1,
  color_edge = "grey",
  size_edge = 0.5,
  node_label = "all",
  foldChange = NULL,
  hilight = "none",
  hilight_alpha = 0.3,)
ggsave("cnetplot_GO Biological Process_600dpi.png", plot = c_BP, width = 8, height = 6, dpi = 600)

c_kegg <- cnetplot(
  kegg,
  layout = igraph::layout_with_kk,
  showCategory = 5,
  color_category = "#E5C494",
  size_category = 1,
  color_item = "#B3B3B3",
  size_item = 1,
  color_edge = "grey",
  size_edge = 0.5,
  node_label = "all",
  foldChange = NULL,
  hilight = "none",
  hilight_alpha = 0.3,)
ggsave("cnet_KEGG Pathway Enrichment_600dpi.png", plot = c_kegg, width = 8, height = 6, dpi = 600)

c_react_p <- cnetplot(
  react_p,
  layout = igraph::layout_with_kk,
  showCategory = 5,
  color_category = "#E5C494",
  size_category = 1,
  color_item = "#B3B3B3",
  size_item = 1,
  color_edge = "grey",
  size_edge = 0.5,
  node_label = "all",
  foldChange = NULL,
  hilight = "none",
  hilight_alpha = 0.3,)
ggsave("cnet_Reactome Pathway Enrichment_600dpi.png", plot = c_react_p, width = 8, height = 6, dpi = 600)