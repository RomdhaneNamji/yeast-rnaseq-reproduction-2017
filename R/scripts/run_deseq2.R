
library(DESeq2)
library(ggplot2)
library(pheatmap)

# ========= Paths ==========
counts_file <- "data/gene_counts.txt"
meta_file <- "data/sample_metadata.csv"

dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# =========== Read metadata ========
meta <- read.csv("data/sample_metadata.csv", stringsAsFactors = FALSE) 

# Keep only the columns needed for deseq2
meta <- meta[,c("sample","condition")]

# Set factor levels : aerial= reference
meta$condtion <- factor(meta$condition, levels = c("aerial","root"))

# ============= Read featureCounts matrix ======
counts_raw <- read.table(
  "data/gene_counts.txt",
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE)

count_data <- counts_raw[,c("Geneid", colnames(counts_raw)[7:ncol(counts_raw)])]

# Clean sample column names:
# /mnt/d/.../AE1.sorted.bam  -> AE1
clean_names <- basename(colnames(count_data)[-1])
clean_names <- sub("\\.sorted\\.bam$", "", clean_names)
colnames(count_data) <- c("Geneid", clean_names)

rownames(count_data) <- count_data$Geneid
count_data <- count_data[,-1]

count_data <- as.matrix(count_data)
storage.mode(count_data) <- "integer"

#======= MAtch metadata to count matrix ===
meta[match(colnames(count_data), meta$sample),]

stopifnot(all(colnames(count_data)== meta$sample))
stopifnot(!any(is.na(meta$sample)))

# ==== Filter low-count genes ====
# Keep genes with at least 10 counts total


keep <- rowSums(count_data) >= 10
count_data <- count_data[keep,]

cat("Number of genes after filtering :", nrow(count_data),"\n")


# ==== Building DESeq2 object =====

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData   = meta,
  design    = ~ condition
)

dds <- DESeq(dds)


# ==== Results ========
res <-  results(dds, contrast = c("condition", "root", "aerial"))
res <- res[order(res$padj),]

# Add gene IDs as a column
res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)
res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]


# save the results
write.csv(res_df, "results/root_vs_aerial.csv", row.names = FALSE)

# save significant genes 
sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange)>=1)
write.csv(sig, "results/root_vs_aerial_sig.csv", row.names = FALSE)



# =========== Normalized counts ==========
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_df <-  as.data.frame(norm_counts)
norm_counts_df$Geneid <- rownames(norm_counts_df)
norm_counts_df <- norm_counts_df[, c("Geneid", setdiff(colnames(norm_counts_df), "Geneid"))]

write.csv(norm_counts_df, "results/normalized_counts.csv", row.names = FALSE)

# ====== PCA ===========
rld <- rlog(dds, blind = TRUE)

pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 2) +
  geom_text(vjust = -1)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("DESeq2 PCA")

ggsave("plots/pca_deseq2.png", plot=pca_plot , width = 7, height = 5, dpi = 300 )


# ============= heatmap ===================
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(rld)
colnames(sample_dist_matrix) <- colnames(rld)

png("plots/sample_distance_heatmap.png", width = 800, height = 700)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         main = "Sample-to-sample distances")
dev.off()


#====== MA plot =======================
png("plots/ma_plot.png", width = 800, height = 700)
plotMA(res, main = "DESeq2 MA plot: root vs aerial", ylim = c(-5, 5))
dev.off()

# -------------------------
# Summary
# -------------------------
sink("results/deseq2_summary.txt")
cat("DESeq2 summary\n")
cat("====================\n")
cat("Samples:\n")
print(meta)
cat("\nGenes after filtering:", nrow(countdata), "\n\n")
cat("Results summary:\n")
print(summary(res))
cat("\nSignificant genes (padj < 0.05 and |log2FC| >= 1): ", nrow(sig), "\n")
sink()

cat("DESeq2 analysis complete.\n")
cat("Results written to results/\n")
cat("Plots written to plots/\n")




# =========== GO analysis =================

# =========================================================
# PREPARE CLEAN GENE LISTS FOR GO ENRICHMENT
# Assumes these already exist in the R session:
#   - res_df   (DESeq2 results data frame with Geneid, log2FoldChange, padj)
#   - dds      (DESeqDataSet object used for DESeq2)
# =========================================================

# -----------------------------
# 1) Create directional DE lists
# -----------------------------
aerial_up <- subset(
  res_df,
  !is.na(padj) & padj < 0.05 & log2FoldChange < -1
)

root_up <- subset(
  res_df,
  !is.na(padj) & padj < 0.05 & log2FoldChange > 1
)

cat("Aerial-up genes:", nrow(aerial_up), "\n")
cat("Root-up genes:", nrow(root_up), "\n")

# -----------------------------
# 2) Clean DESeq2 gene IDs
# Example: gene:YAL069W -> YAL069W
# -----------------------------
aerial_ids <- sub("^gene:", "", aerial_up$Geneid)
root_ids   <- sub("^gene:", "", root_up$Geneid)

background_ids <- rownames(dds)
background_ids <- sub("^gene:", "", background_ids)

# Remove duplicates just in case
aerial_ids <- unique(aerial_ids)
root_ids <- unique(root_ids)
background_ids <- unique(background_ids)

cat("Background genes:", length(background_ids), "\n")
cat("Example aerial IDs:\n")
print(head(aerial_ids, 10))
cat("Example root IDs:\n")
print(head(root_ids, 10))

# -----------------------------
# 3) Save raw directional lists
# -----------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

write.table(
  aerial_ids,
  "results/aerial_up_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(
  root_ids,
  "results/root_up_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(
  background_ids,
  "results/background_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

# -----------------------------
# 4) Read GTF
# IMPORTANT: adjust this path if needed
# -----------------------------
gtf_file <- "/mnt/d/GIS/yeast_metabolomics_2017/ref/downloads/sacCer3.gtf"

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

cat("GTF feature counts:\n")
print(table(gtf$feature))

# -----------------------------
# 5) Use CDS-containing genes as protein-coding genes
# -----------------------------
gtf_cds <- subset(gtf, feature == "CDS")

cat("Number of CDS rows:", nrow(gtf_cds), "\n")
cat("Example CDS attributes:\n")
print(head(gtf_cds$attribute, 10))

# Attribute format in your file looks like:
# transcript_id transcript:YAL069W_mRNA; gene_id gene:YAL069W;
#
# So extract gene IDs in two simple steps:
pc_gene_ids <- sub(".*gene_id gene:", "", gtf_cds$attribute)
pc_gene_ids <- sub(";.*", "", pc_gene_ids)

# Clean and deduplicate
pc_gene_ids <- trimws(pc_gene_ids)
pc_gene_ids <- pc_gene_ids[nchar(pc_gene_ids) > 0]
pc_gene_ids <- unique(pc_gene_ids)

cat("Protein-coding gene IDs from CDS:", length(pc_gene_ids), "\n")
cat("Example protein-coding IDs:\n")
print(head(pc_gene_ids, 20))

# -----------------------------
# 6) Intersect DESeq2-tested genes with protein-coding genes
# -----------------------------
background_pc <- intersect(background_ids, pc_gene_ids)
aerial_ids_pc <- intersect(aerial_ids, pc_gene_ids)
root_ids_pc   <- intersect(root_ids, pc_gene_ids)

cat("Protein-coding background:", length(background_pc), "\n")
cat("Protein-coding aerial-up:", length(aerial_ids_pc), "\n")
cat("Protein-coding root-up:", length(root_ids_pc), "\n")

cat("Example aerial protein-coding IDs:\n")
print(head(aerial_ids_pc, 10))
cat("Example root protein-coding IDs:\n")
print(head(root_ids_pc, 10))

cat("Are aerial and root protein-coding lists identical? ",
    identical(aerial_ids_pc, root_ids_pc), "\n")
cat("Overlap size between aerial and root protein-coding lists: ",
    length(intersect(aerial_ids_pc, root_ids_pc)), "\n")

# -----------------------------
# 7) Save clean protein-coding lists
# -----------------------------
write.table(
  background_pc,
  "results/background_pc_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(
  aerial_ids_pc,
  "results/aerial_up_pc_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

write.table(
  root_ids_pc,
  "results/root_up_pc_ids.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)

cat("\nDone. Clean gene lists saved in results/.\n")
