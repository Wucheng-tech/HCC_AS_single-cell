############################################
# Differential alternative splicing function
############################################
# Parameters:
# psi_matrix   : PSI matrix (rows = AS events, columns = cells)
# meta         : metadata dataframe (rows = cells)
# group.by     : column name in meta defining comparison groups
# group1       : group name (e.g. "T")
# group2       : group name (e.g. "N")
# min_cells    : minimum number of valid cells per group
# pval_thresh  : p-value cutoff
# delta_thresh : ΔPSI cutoff
# method       : statistical test ("wilcox" or "t.test")
# pseudo       : pseudo count for fold change calculation
# expr_thresh  : threshold used for pct calculation
library("Seurat")
FindASMarkers <- function(psi_matrix, meta, group.by, group1, group2, min_cells = 3, pval_thresh = 0.05, delta_thresh = 0.1, method = "wilcox", pseudo = 0.001, expr_thresh = 0) {
  if (!(group.by %in% colnames(meta))) {
    stop(paste("Group variable", group.by, "not found in meta data"))
  }
  group_vector <- meta[, group.by]
  names(group_vector) <- rownames(meta)
  group1_cells <- intersect(names(group_vector)[group_vector == group1], colnames(psi_matrix))
  group2_cells <- intersect(names(group_vector)[group_vector == group2], colnames(psi_matrix))
  if (length(group1_cells) < min_cells | length(group2_cells) < min_cells) {
    stop("Too few cells in one of the groups")
  }
  results <- apply(psi_matrix, 1, function(x) {
    x1 <- as.numeric(x[group1_cells])
    x2 <- as.numeric(x[group2_cells])
    x1_nonNA <- x1[!is.na(x1)]
    x2_nonNA <- x2[!is.na(x2)]
    if (length(x1_nonNA) < min_cells || length(x2_nonNA) < min_cells) {
      return(c(NA, NA, NA, NA, NA, NA))
    }
    pct1 <- sum(x1_nonNA > expr_thresh, na.rm = TRUE) / length(x1_nonNA)
    pct2 <- sum(x2_nonNA > expr_thresh, na.rm = TRUE) / length(x2_nonNA)
    if (method == "t.test") {
      test <- t.test(x1_nonNA, x2_nonNA)
      pval <- test$p.value
    } else if (method == "wilcox") {
      test <- wilcox.test(x1_nonNA, x2_nonNA)
      pval <- test$p.value
    } else {
      stop("Unknown method. Use 't.test' or 'wilcox'")
    }
    mean1 <- mean(x1_nonNA)
    mean2 <- mean(x2_nonNA)
    delta <- mean1 - mean2
    fc <- (mean1 + pseudo) / (mean2 + pseudo)
    logFC <- log2(fc)
    return(c(pval, delta, fc, logFC, pct1, pct2))
  })
  results_df <- as.data.frame(t(results))
  colnames(results_df) <- c("p_val", "delta_PSI", "FC", "log2FC", "pct.1", "pct.2")
  results_df$p_val <- as.numeric(results_df$p_val)
  results_df$delta_PSI <- as.numeric(results_df$delta_PSI)
  results_df$FC <- as.numeric(results_df$FC)
  results_df$log2FC <- as.numeric(results_df$log2FC)
  results_df$pct.1 <- as.numeric(results_df$pct.1)
  results_df$pct.2 <- as.numeric(results_df$pct.2)
  results_df$status <- "no_change"
  results_df$status[results_df$p_val < pval_thresh & results_df$delta_PSI > delta_thresh] <- "up"
  results_df$status[results_df$p_val < pval_thresh & results_df$delta_PSI < -delta_thresh] <- "down"
  return(results_df)
}

############################################
# Cell-type specific AS analysis
############################################

setwd("/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/celltype_bin_100/Seurat")
psi_matrix <- read.csv("AS_abs.csv",row.names = 1,check.names = FALSE)
meta <- read.table("Meta_Gene_AS.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
celltypes <- c("Hepatocyte","T","Mye","Fib","Edno","B","pB","NK","Pro")
for (ct in celltypes) {
  others <- setdiff(celltypes, ct)
  res_status <- FindASMarkers(
    psi_matrix, meta,
    group.by = "celltype",
    group1 = ct,
    group2 = others
  )  
  res_status1 <- res_status[res_status$status == "up" &(res_status$pct.1 >= 0.1 | res_status$pct.2 >= 0.1),]  
  write.table(res_status1,paste0(ct, ".txt"),sep = "\t",quote = FALSE)  
  cat(ct, "done | significant:", nrow(res_status1), "\n")
}

############################################
# Tumor vs Normal AS analysis
############################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
psi_matrix <- read.csv("AS_abs.csv",row.names = 1,check.names = FALSE)
meta <- read.table("Meta_Gene_AS.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
meta$status <- ifelse(grepl("T1|T2", meta$sample),"T", ifelse(grepl("PN", meta$sample), "N", NA))
clusters <- c("Hepatocyte", "T", "NK", "Mye", "B", "pB", "Fib", "Edno", "Pro")
base_dir <- "/public/home/wucheng/Analysis/AS/scRNA/scFAST/Sample/AA/celltype_bin_100/Seurat"
for (j in clusters) {
  cat("\nProcessing:", j, "\n")
  outdir <- file.path(base_dir, j)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  cells_use <- rownames(meta)[meta$celltype == j]
  if (length(cells_use) < 10) {
    cat("skip:", j, "(too few cells)\n")
    next
  }
  psi_sub <- psi_matrix[, cells_use, drop = FALSE]
  meta_sub <- meta[cells_use, , drop = FALSE]
  keep <- !is.na(meta_sub$status)
  psi_sub <- psi_sub[, keep, drop = FALSE]
  meta_sub <- meta_sub[keep, , drop = FALSE]
  if (length(unique(meta_sub$status)) < 2) {
    cat("skip:", j, "(missing T or N)\n")
    next
  }
  res_status <- FindASMarkers( psi_sub,meta_sub,group.by = "status",group1 = "T", group2 = "N"
  )
  res_up <- res_status[
    res_status[,7] == "up" &
    (res_status[,5] >= 0.1 | res_status[,6] >= 0.1),
  ]
  write.table(res_up,file.path(outdir, paste0(j, "_T_vs_N_up.txt")),sep = "\t",quote = FALSE)
  res_down <- res_status[
    res_status[,7] == "down" &
    (res_status[,5] >= 0.1 | res_status[,6] >= 0.1),
  ]
  write.table(res_down,file.path(outdir, paste0(j, "_T_vs_N_down.txt")),sep = "\t",quote = FALSE)
}  



