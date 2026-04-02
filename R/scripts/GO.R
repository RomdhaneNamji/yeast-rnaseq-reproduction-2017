if (!file.exists(gaf_file)) stop("Missing file: ", gaf_file)
if (!file.exists(obo_file)) stop("Missing file: ", obo_file)

gaf <- read.delim(
  gaf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "!",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(gaf) <- c(
  "DB",
  "DB_Object_ID",
  "DB_Object_Symbol",
  "Qualifier",
  "GO_ID",
  "DB_Reference",
  "Evidence_Code",
  "With_or_From",
  "Aspect",
  "DB_Object_Name",
  "DB_Object_Synonym",
  "DB_Object_Type",
  "Taxon",
  "Date",
  "Assigned_By",
  "Annotation_Extension",
  "Gene_Product_Form_ID"
)

# Biological Process only
gaf_bp <- subset(gaf, Aspect == "P")

# This is the column that matched your genes
TERM2GENE <- unique(gaf_bp[, c("GO_ID", "DB_Object_Symbol")])
colnames(TERM2GENE) <- c("term", "gene")

# Read OBO
obo_lines <- readLines(obo_file)
term_starts <- which(obo_lines == "[Term]")
go_terms <- list()

for (i in seq_along(term_starts)) {
  start <- term_starts[i]
  end <- if (i < length(term_starts)) term_starts[i + 1] - 1 else length(obo_lines)
  block <- obo_lines[start:end]
  
  id_line <- grep("^id: GO:", block, value = TRUE)
  name_line <- grep("^name: ", block, value = TRUE)
  namespace_line <- grep("^namespace: ", block, value = TRUE)
  obsolete_line <- grep("^is_obsolete: true", block, value = TRUE)
  
  if (length(id_line) == 1 &&
      length(name_line) == 1 &&
      length(namespace_line) == 1 &&
      length(obsolete_line) == 0) {
    
    go_id <- sub("^id: ", "", id_line)
    go_name <- sub("^name: ", "", name_line)
    go_namespace <- sub("^namespace: ", "", namespace_line)
    
    go_terms[[length(go_terms) + 1]] <- data.frame(
      GO_ID = go_id,
      GO_NAME = go_name,
      NAMESPACE = go_namespace,
      stringsAsFactors = FALSE
    )
  }
}

go_df <- do.call(rbind, go_terms)
go_bp <- subset(go_df, NAMESPACE == "biological_process")

TERM2NAME <- unique(go_bp[, c("GO_ID", "GO_NAME")])
colnames(TERM2NAME) <- c("term", "name")

common_terms <- intersect(TERM2GENE$term, TERM2NAME$term)
TERM2GENE <- subset(TERM2GENE, term %in% common_terms)
TERM2NAME <- subset(TERM2NAME, term %in% common_terms)

write.csv(TERM2GENE, file.path(results_dir, "TERM2GENE_BP.csv"), row.names = FALSE)
write.csv(TERM2NAME, file.path(results_dir, "TERM2NAME_BP.csv"), row.names = FALSE)

cat("\nGO TERM2GENE dimensions:", dim(TERM2GENE), "\n")
cat("GO TERM2NAME dimensions:", dim(TERM2NAME), "\n")
cat("GO coverage in background:", sum(unique(TERM2GENE$gene) %in% background_pc), "\n")
cat("GO coverage in aerial-up:", sum(unique(TERM2GENE$gene) %in% aerial_ids_pc), "\n")
cat("GO coverage in root-up:", sum(unique(TERM2GENE$gene) %in% root_ids_pc), "\n")

# =========================================================
# 9) GO enrichment with enrichR
# =========================================================
library(enrichR)
dbs <- enrichR::listEnrichrDbs()
dbs[grep("yeast|saccharomyces|go", dbs$libraryName, ignore.case = TRUE), c("libraryName", "categoryId")]

db_use <-  c("GO_Biological_Process_2025")

aerial_enrich <- enrichR::enrichr(aerial_ids_pc, db_use)
aerial_res <-  aerial_enrich[[1]]
head(aerial_res,20)

root_enrich <- enrichR::enrichr(root_ids_pc, db_use)
root_res <- root_enrich[[1]]
head(root_res, 20)

write.csv(aerial_res, "results/enrichr_aerial_GO_BP_2025.csv", row.names = FALSE)
write.csv(root_res,   "results/enrichr_root_GO_BP_2025.csv",   row.names = FALSE)

aerial_res <- aerial_res[order(aerial_res$Adjusted.P.value), ]
root_res   <- root_res[order(root_res$Adjusted.P.value), ]

head(aerial_res[, c("Term", "Adjusted.P.value", "Overlap", "Odds.Ratio", "Combined.Score")], 20)
head(root_res[, c("Term", "Adjusted.P.value", "Overlap", "Odds.Ratio", "Combined.Score")], 20)


run_go_fisher <- function(gene_list, universe, TERM2GENE, TERM2NAME, min_genes = 2) {
  gene_list <- unique(gene_list)
  universe <- unique(universe)
  
  # keep only genes in universe
  gene_list <- intersect(gene_list, universe)
  term2gene <- TERM2GENE[TERM2GENE$gene %in% universe, ]
  
  terms <- unique(term2gene$term)
  
  results <- lapply(terms, function(term) {
    term_genes <- unique(term2gene$gene[term2gene$term == term])
    
    overlap_genes <- intersect(gene_list, term_genes)
    
    a <- length(overlap_genes)
    b <- length(gene_list) - a
    c <- length(term_genes) - a
    d <- length(universe) - a - b - c
    
    if (a < min_genes) return(NULL)
    if (any(c(a, b, c, d) < 0)) return(NULL)
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    pval <- fisher.test(mat, alternative = "greater")$p.value
    
    data.frame(
      term = term,
      overlap = a,
      gene_set_size = length(term_genes),
      query_size = length(gene_list),
      universe_size = length(universe),
      p.value = pval,
      genes = paste(sort(overlap_genes), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  
  results <- do.call(rbind, results)
  
  if (is.null(results) || nrow(results) == 0) {
    return(data.frame())
  }
  
  results$p.adjust <- p.adjust(results$p.value, method = "BH")
  
  results <- merge(results, TERM2NAME, by = "term", all.x = TRUE)
  results <- results[order(results$p.adjust, results$p.value), ]
  
  rownames(results) <- NULL
  results
}

aerial_go <- run_go_fisher(
  gene_list = aerial_ids_pc,
  universe = background_pc,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  min_genes = 2
)

root_go <- run_go_fisher(
  gene_list = root_ids_pc,
  universe = background_pc,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME,
  min_genes = 2
)

dim(aerial_go)
dim(root_go)

head(aerial_go[, c("term", "name", "overlap", "gene_set_size", "p.value", "p.adjust", "genes")], 20)
head(root_go[, c("term", "name", "overlap", "gene_set_size", "p.value", "p.adjust", "genes")], 20)

write.csv(aerial_go, "results/local_GO_BP_aerial.csv", row.names = FALSE)
write.csv(root_go,   "results/local_GO_BP_root.csv",   row.names = FALSE)
