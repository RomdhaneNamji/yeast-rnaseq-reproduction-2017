
library(DESeq2)
library(ggplot2)
library(pheatmap)


# =========================================================
# MAKE ALL PLOTS FOR THE YEAST RNA-SEQ PROJECT
# =========================================================

# -------------------------
# Paths
# -------------------------
project_dir <- "D:/GIS/yeast_metabolomics_2017/R"

counts_file <- file.path(project_dir, "data", "gene_counts.txt")
meta_file   <- file.path(project_dir, "data", "sample_metadata.csv")
results_dir <- file.path(project_dir, "results")
plots_dir   <- file.path(project_dir, "plots")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Load metadata
# -------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[, c("sample", "condition")]
meta$condition <- factor(meta$condition, levels = c("aerial", "root"))

# -------------------------
# Load count matrix
# -------------------------
counts_raw <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

countdata <- counts_raw[, c("Geneid", colnames(counts_raw)[7:ncol(counts_raw)])]

clean_names <- basename(colnames(countdata)[-1])
clean_names <- sub("\\.sorted\\.bam$", "", clean_names)
colnames(countdata) <- c("Geneid", clean_names)

rownames(countdata) <- countdata$Geneid
countdata <- countdata[, -1]
countdata <- as.matrix(countdata)
storage.mode(countdata) <- "integer"

meta <- meta[match(colnames(countdata), meta$sample), ]
stopifnot(all(colnames(countdata) == meta$sample))

# -------------------------
# Filter and run DESeq2
# -------------------------
keep <- rowSums(countdata) >= 10
countdata <- countdata[keep, ]

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData   = meta,
  design    = ~ condition
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "root", "aerial"))
res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)

# Clean gene IDs for plotting labels
res_df$GeneClean <- sub("^gene:", "", res_df$Geneid)

# rlog for PCA / heatmaps
rld <- rlog(dds, blind = TRUE)

# =========================================================
# 1) PCA plot
# =========================================================
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of samples") +
  theme_bw()

ggsave(file.path(plots_dir, "01_PCA_samples.png"), plot = p_pca, width = 7, height = 5, dpi = 300)

# =========================================================
# 2) Sample distance heatmap
# =========================================================
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(rld)
colnames(sample_dist_matrix) <- colnames(rld)

png(file.path(plots_dir, "02_sample_distance_heatmap.png"), width = 900, height = 800)
tryCatch({
  pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    main = "Sample-to-sample distances"
  )
}, finally = {
  dev.off()
})

# =========================================================
# 3) MA plot
# =========================================================
png(file.path(plots_dir, "03_MA_plot.png"), width = 900, height = 800)
plotMA(res, main = "MA plot: root vs aerial", ylim = c(-5, 5))
dev.off()

# =========================================================
# 4) Volcano plot
# =========================================================
volcano_df <- res_df
volcano_df$neglog10padj <- -log10(volcano_df$padj)
volcano_df$neglog10padj[!is.finite(volcano_df$neglog10padj)] <- NA

volcano_df$category <- "Not significant"
volcano_df$category[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 1] <- "Up in root"
volcano_df$category[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange < -1] <- "Up in aerial"

top_labels <- subset(volcano_df, !is.na(padj) & padj < 0.01 & abs(log2FoldChange) > 2)
top_labels <- head(top_labels[order(top_labels$padj), ], 15)

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, neglog10padj, color = category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text(
    data = top_labels,
    aes(label = GeneClean),
    vjust = -0.5,
    size = 3,
    show.legend = FALSE
  ) +
  xlab("log2 fold change (root vs aerial)") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot") +
  theme_bw()

ggsave(file.path(plots_dir, "04_volcano_plot.png"), plot = p_volcano, width = 8, height = 6, dpi = 300)

# =========================================================
# 5) Heatmap of top DE genes
# =========================================================
sig_df <- subset(res_df, !is.na(padj) & padj < 0.05)
top_up_root   <- head(sig_df[order(-sig_df$log2FoldChange), "Geneid"], 20)
top_up_aerial <- head(sig_df[order(sig_df$log2FoldChange), "Geneid"], 20)

top_genes <- unique(c(top_up_root, top_up_aerial))
top_mat <- assay(rld)[top_genes, , drop = FALSE]

# nicer rownames
rownames(top_mat) <- sub("^gene:", "", rownames(top_mat))

annotation_col <- data.frame(condition = meta$condition)
rownames(annotation_col) <- meta$sample

png(file.path(plots_dir, "05_top_DE_heatmap.png"), width = 1000, height = 1200)
tryCatch({
  pheatmap(
    top_mat,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    main = "Top DE genes"
  )
}, finally = {
  dev.off()
})

# =========================================================
# 6) Boxplot of normalized counts
# =========================================================
norm_counts <- counts(dds, normalized = TRUE)
log_norm <- log2(norm_counts + 1)

png(file.path(plots_dir, "06_normalized_counts_boxplot.png"), width = 1000, height = 700)
boxplot(
  as.data.frame(log_norm),
  las = 2,
  main = "Log2 normalized counts per sample",
  ylab = "log2(normalized counts + 1)"
)
dev.off()

# =========================================================
# 7) Dispersion plot
# =========================================================
png(file.path(plots_dir, "07_dispersion_plot.png"), width = 900, height = 800)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# =========================================================
# 8) Marker heatmap based on enrichment-related themes
# translation / stress / sporulation / transport keywords
# =========================================================
theme_keywords <- c("translation", "stress", "sporulation", "transport")

theme_hits <- res_df[grep(paste(theme_keywords, collapse = "|"), res_df$GeneClean, ignore.case = TRUE), ]
theme_genes <- unique(head(theme_hits$Geneid, 30))

if (length(theme_genes) > 1) {
  theme_mat <- assay(rld)[theme_genes, , drop = FALSE]
  rownames(theme_mat) <- sub("^gene:", "", rownames(theme_mat))
  
  png(file.path(plots_dir, "08_theme_gene_heatmap.png"), width = 1000, height = 900)
  tryCatch({
    pheatmap(
      theme_mat,
      scale = "row",
      annotation_col = annotation_col,
      show_rownames = TRUE,
      main = "Selected theme-related genes"
    )
  }, finally = {
    dev.off()
  })
}

# =========================================================
# 9) GO barplots from local GO enrichment results, if present
# =========================================================
aerial_go_file <- file.path(results_dir, "local_GO_BP_aerial.csv")
root_go_file   <- file.path(results_dir, "local_GO_BP_root.csv")

if (file.exists(aerial_go_file)) {
  aerial_go <- read.csv(aerial_go_file, stringsAsFactors = FALSE)
  aerial_go <- aerial_go[order(aerial_go$p.adjust, aerial_go$p.value), ]
  top_aerial_go <- head(aerial_go, 10)
  
  if (nrow(top_aerial_go) > 0) {
    top_aerial_go$label <- factor(top_aerial_go$name, levels = rev(top_aerial_go$name))
    
    p_go_aerial <- ggplot(top_aerial_go, aes(x = label, y = -log10(p.adjust))) +
      geom_col() +
      coord_flip() +
      xlab("") +
      ylab("-log10 adjusted p-value") +
      ggtitle("Top GO BP terms: aerial-up genes") +
      theme_bw()
    
    ggsave(file.path(plots_dir, "09_GO_barplot_aerial.png"), plot = p_go_aerial, width = 9, height = 6, dpi = 300)
  }
}

if (file.exists(root_go_file)) {
  root_go <- read.csv(root_go_file, stringsAsFactors = FALSE)
  root_go <- root_go[order(root_go$p.adjust, root_go$p.value), ]
  top_root_go <- head(root_go, 10)
  
  if (nrow(top_root_go) > 0) {
    top_root_go$label <- factor(top_root_go$name, levels = rev(top_root_go$name))
    
    p_go_root <- ggplot(top_root_go, aes(x = label, y = -log10(p.adjust))) +
      geom_col() +
      coord_flip() +
      xlab("") +
      ylab("-log10 adjusted p-value") +
      ggtitle("Top GO BP terms: root-up genes") +
      theme_bw()
    
    ggsave(file.path(plots_dir, "10_GO_barplot_root.png"), plot = p_go_root, width = 9, height = 6, dpi = 300)
  }
}

# =========================================================
# 10) Save a concise plotting summary
# =========================================================
sink(file.path(results_dir, "plot_summary.txt"))
cat("Plot summary\n")
cat("============\n")
cat("Samples:\n")
print(meta)
cat("\nGenes after filtering:", nrow(countdata), "\n")
cat("Significant genes (padj < 0.05):", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n")
cat("Top root-up gene examples:\n")
print(head(sig_df[order(-sig_df$log2FoldChange), c("GeneClean", "log2FoldChange", "padj")], 10))
cat("\nTop aerial-up gene examples:\n")
print(head(sig_df[order(sig_df$log2FoldChange), c("GeneClean", "log2FoldChange", "padj")], 10))
sink()

cat("All plots generated in:\n", plots_dir, "\n")