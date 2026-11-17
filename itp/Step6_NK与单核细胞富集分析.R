#######################################################
# 拓展KEGG和GSEA富集分析
# Monocytes + NK: KEGG扩展 + GSEA
#######################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(msigdbr)
  library(ggridges)
  library(DOSE)
})

# ==== Configuration ====
BASE_DIR <- getwd()
RESULT_DIR <- file.path(BASE_DIR, "差异富集 result")

CELL_TYPES <- c("Monocytes", "NK")
FOLDER_NAMES <- c("Monocytes" = "单核细胞", "NK" = "NK细胞")

# ==== Helper: KEGG enrichment with expanded plots ====
run_kegg_expanded <- function(gene_list, logFC_vec, ct_name, table_dir, plot_dir) {
  names(logFC_vec) <- gene_list

  sig_genes <- names(logFC_vec)[abs(logFC_vec) > 0.25]
  if (length(sig_genes) < 5) {
    cat("  KEGG ORA: insufficient genes (<5), skipping\n")
    return(NULL)
  }

  gene_ids <- tryCatch({
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)

  if (is.null(gene_ids) || nrow(gene_ids) < 5) {
    cat("  KEGG ORA: ID conversion failed\n")
    return(NULL)
  }

  ekegg <- tryCatch({
    enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    kegg_df <- as.data.frame(ekegg)
    cat(sprintf("  KEGG ORA: %d enriched pathways\n", nrow(kegg_df)))

    write.csv(kegg_df, file.path(table_dir, "KEGG_enrichment.csv"), row.names = FALSE)

    if (nrow(kegg_df) >= 3) {
      p_dot <- dotplot(ekegg, showCategory = 15, title = paste(ct_name, "- KEGG Dotplot")) +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
      ggsave(file.path(plot_dir, "KEGG_dotplot.pdf"), p_dot, width = 12, height = 8)

      p_bar <- barplot(ekegg, showCategory = 15, title = paste(ct_name, "- KEGG Barplot")) +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
      ggsave(file.path(plot_dir, "KEGG_barplot.pdf"), p_bar, width = 12, height = 8)
    }

    if (nrow(kegg_df) >= 3 && nrow(kegg_df) <= 30) {
      tryCatch({
        gene_fc <- logFC_vec[sig_genes]
        p_cnet <- cnetplot(ekegg, showCategory = min(5, nrow(kegg_df)),
                          foldChange = gene_fc, circular = FALSE,
                          colorEdge = TRUE, cex_label_gene = 0.5,
                          cex_label_category = 0.7) +
          ggtitle(paste(ct_name, "- KEGG Network")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        ggsave(file.path(plot_dir, "KEGG_cnetplot.pdf"), p_cnet, width = 14, height = 10)
      }, error = function(e) {
        cat(sprintf("  KEGG cnetplot failed: %s\n", e$message))
      })
    }

    if (nrow(kegg_df) >= 5) {
      tryCatch({
        ekegg_pair <- pairwise_termsim(ekegg)
        p_emap <- emapplot(ekegg_pair, showCategory = min(15, nrow(kegg_df)),
                          cex_label_category = 0.7) +
          ggtitle(paste(ct_name, "- KEGG Enrichment Map")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        ggsave(file.path(plot_dir, "KEGG_emapplot.pdf"), p_emap, width = 12, height = 10)
      }, error = function(e) {
        cat(sprintf("  KEGG emapplot failed: %s\n", e$message))
      })
    }

    return(ekegg)
  } else {
    cat("  KEGG ORA: no enriched pathways\n")
    return(NULL)
  }
}
# ==== Helper: Load human MSigDB collections once for both cell types ====
read_msigdb_collection <- function(label, collection, subcollection = NULL) {
  tryCatch({
    sets <- msigdbr(species = "Homo sapiens", collection = collection,
                    subcollection = subcollection)
    gene_sets <- split(sets$gene_symbol, sets$gs_name)
    cat(sprintf("  MSigDB %s: %d gene sets (%s)\n",
                label, length(gene_sets), paste(unique(sets$db_version), collapse = ", ")))
    gene_sets
  }, error = function(e) {
    message("  无法加载 MSigDB ", label, "，跳过这一集合：", conditionMessage(e))
    NULL
  })
}

# ==== Helper: GSEA with fgsea ====
run_gsea_fgsea <- function(ranked_genes, ct_name, table_dir, plot_dir, gene_sets) {
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  ranked_genes <- ranked_genes[!is.na(ranked_genes) & is.finite(ranked_genes)]
  cat(sprintf("  GSEA: %d genes in ranked list\n", length(ranked_genes)))

  all_results <- list()
  for (gs_name in names(gene_sets)) {
    gs <- gene_sets[[gs_name]]
    if (is.null(gs) || length(gs) < 5) next

    cat(sprintf("  Running fgsea for %s (%d gene sets)...\n", gs_name, length(gs)))
    fgsea_res <- tryCatch({
      fgsea(pathways = gs, stats = ranked_genes,
            minSize = 10, maxSize = 500, nPermSimple = 10000)
    }, error = function(e) {
      cat(sprintf("  fgsea %s failed: %s\n", gs_name, e$message))
      return(NULL)
    })
    if (is.null(fgsea_res) || nrow(fgsea_res) == 0) next

    fgsea_res <- fgsea_res[order(padj), ]
    fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = "/")
    write.csv(as.data.frame(fgsea_res)[, !colnames(as.data.frame(fgsea_res)) %in% "leadingEdge_raw"],
              file.path(table_dir, sprintf("GSEA_%s_fgsea.csv", gs_name)), row.names = FALSE)

    n_sig <- sum(fgsea_res$padj < 0.05, na.rm = TRUE)
    cat(sprintf("  fgsea %s: %d/%d significant pathways (padj<0.05)\n", gs_name, n_sig, nrow(fgsea_res)))
    all_results[[gs_name]] <- fgsea_res

    # Enrichment plot for top pathways
    top_sig <- fgsea_res[padj < 0.05][order(padj)]
    if (nrow(top_sig) > 0) {
      top_n <- min(20, nrow(top_sig))
      top_pathways <- top_sig$pathway[1:top_n]
      for (pw in head(top_pathways, 6)) {
        tryCatch({
          safe_pw <- gsub("[^A-Za-z0-9_]", "_", substr(pw, 1, 50))
          pdf(file.path(plot_dir, sprintf("GSEA_%s_%s.pdf", gs_name, safe_pw)),
              width = 8, height = 6)
          p <- plotEnrichment(gs[[pw]], ranked_genes) +
            labs(title = pw, subtitle = sprintf("NES=%.2f, padj=%.2e",
                  fgsea_res[pathway == pw, NES], fgsea_res[pathway == pw, padj])) +
            theme_minimal()
          print(p)
          dev.off()
        }, error = function(e) NULL)
      }

      # Ridge plot
      if (length(top_pathways) >= 3) {
        tryCatch({
          ridge_data <- data.frame()
          for (pw in top_pathways) {
            nes_val <- fgsea_res[pathway == pw, NES]
            padj_val <- fgsea_res[pathway == pw, padj]
            label <- sprintf("%s (NES=%.2f, padj=%.1e)", pw, nes_val, padj_val)
            gs_genes <- intersect(gs[[pw]], names(ranked_genes))
            if (length(gs_genes) > 0) {
              ranks_sub <- ranked_genes[gs_genes]
              ridge_data <- rbind(ridge_data, data.frame(
                pathway = label, rank = ranks_sub, stringsAsFactors = FALSE
              ))
            }
          }
          if (nrow(ridge_data) > 0) {
            p_ridge <- ggplot(ridge_data, aes(x = rank, y = pathway, fill = after_stat(x))) +
              geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
              scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                                  midpoint = 0, name = "log2FC") +
              labs(x = "log2 Fold Change", y = "",
                   title = paste(ct_name, "- GSEA", gs_name, "Ridge Plot")) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5, size = 12))
            ggsave(file.path(plot_dir, sprintf("GSEA_%s_ridge.pdf", gs_name)),
                   p_ridge, width = 14, height = max(4, 2 + length(top_pathways) * 0.6))
          }
        }, error = function(e) NULL)
      }
    }
  }
  return(all_results)
}

# ==== Helper: GSEA with clusterProfiler (gseKEGG/gseGO) ====
run_gsea_clusterprofiler <- function(ranked_genes_entrez, ct_name, table_dir, plot_dir) {
  gsea_kegg <- tryCatch({
    gseKEGG(geneList = ranked_genes_entrez, organism = "hsa",
            minGSSize = 10, maxGSSize = 500,
            pvalueCutoff = 0.05, pAdjustMethod = "BH")
  }, error = function(e) {
    cat(sprintf("  gseKEGG failed: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    gsea_kegg_df <- as.data.frame(gsea_kegg)
    cat(sprintf("  gseKEGG: %d enriched pathways\n", nrow(gsea_kegg_df)))
    write.csv(gsea_kegg_df, file.path(table_dir, "GSEA_KEGG_clusterProfiler.csv"), row.names = FALSE)

    if (nrow(gsea_kegg_df) >= 3) {
      tryCatch({
        p_ridge <- ridgeplot(gsea_kegg, showCategory = 20) +
          ggtitle(paste(ct_name, "- gseKEGG Ridge Plot")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        ggsave(file.path(plot_dir, "GSEA_KEGG_ridge.pdf"), p_ridge, width = 14, height = 10)
      }, error = function(e) NULL)
    }

    top_paths <- head(gsea_kegg_df$ID[gsea_kegg_df$p.adjust < 0.05], 6)
    for (pw_id in top_paths) {
      tryCatch({
        pdf(file.path(plot_dir, sprintf("GSEA_KEGG_%s.pdf", gsub(":", "_", pw_id))),
            width = 8, height = 6)
        p <- gseaplot2(gsea_kegg, geneSetID = pw_id,
                       title = gsea_kegg_df$Description[gsea_kegg_df$ID == pw_id][1],
                       pvalue_table = TRUE)
        print(p)
        dev.off()
      }, error = function(e) NULL)
    }
  } else {
    cat("  gseKEGG: no enriched pathways\n")
  }

  gsea_go <- tryCatch({
    gseGO(geneList = ranked_genes_entrez, ont = "BP", OrgDb = org.Hs.eg.db,
          minGSSize = 10, maxGSSize = 500,
          pvalueCutoff = 0.05, pAdjustMethod = "BH")
  }, error = function(e) {
    cat(sprintf("  gseGO failed: %s\n", e$message))
    return(NULL)
  })

  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    gsea_go_df <- as.data.frame(gsea_go)
    cat(sprintf("  gseGO: %d enriched terms\n", nrow(gsea_go_df)))
    write.csv(gsea_go_df, file.path(table_dir, "GSEA_GO_BP_clusterProfiler.csv"), row.names = FALSE)

    if (nrow(gsea_go_df) >= 3) {
      tryCatch({
        p_ridge_go <- ridgeplot(gsea_go, showCategory = 20) +
          ggtitle(paste(ct_name, "- gseGO BP Ridge Plot")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14))
        ggsave(file.path(plot_dir, "GSEA_GO_BP_ridge.pdf"), p_ridge_go, width = 14, height = 10)
      }, error = function(e) NULL)
    }
  } else {
    cat("  gseGO: no enriched terms\n")
  }

  return(list(KEGG = gsea_kegg, GO = gsea_go))
}
# ==== Main: Process each cell type ====
msigdb_sets <- NULL
for (ct in CELL_TYPES) {
  folder <- FOLDER_NAMES[ct]
  table_dir <- file.path(RESULT_DIR, folder, "05_differential")
  plot_dir <- file.path(RESULT_DIR, folder, "06_visualization")

  cat(sprintf("\n========================================\n"))
  cat(sprintf("Cell Type: %s\n", ct))
  cat(sprintf("========================================\n"))

  pb_file <- file.path(table_dir, "DE_pseudobulk_DESeq2.csv")
  if (!file.exists(pb_file)) {
    cat(sprintf("  ERROR: Pseudobulk DE file not found: %s\n", pb_file))
    next
  }
  pb_de <- read.csv(pb_file, stringsAsFactors = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("  Loaded pseudobulk DE: %d genes\n", nrow(pb_de)))

  if (!all(c("gene", "log2FoldChange", "padj") %in% names(pb_de)) ||
      !is.numeric(pb_de$log2FoldChange) || !is.numeric(pb_de$padj)) {
    stop("配对差异表应包含 gene、数值列 log2FoldChange 和 padj：", pb_file)
  }

  # KEGG ORA 和 GSEA 都使用 Step5 的配对差异结果。
  # 保留原有基因符号，只转换 Ensembl ID。重复时优先使用原有符号的结果。
  is_ensembl <- grepl("^ENSG[0-9]+$", pb_de$gene)
  native_rows <- which(!is_ensembl & !is.na(pb_de$gene) & nzchar(pb_de$gene))
  symbol_map <- data.frame(ENSEMBL = character(), SYMBOL = character())
  if (any(is_ensembl)) {
    symbol_map <- tryCatch({
      bitr(pb_de$gene[is_ensembl], fromType = "ENSEMBL", toType = "SYMBOL",
           OrgDb = org.Hs.eg.db)
    }, error = function(e) {
      cat(sprintf("  ENSG->SYMBOL mapping failed: %s\n", e$message))
      data.frame(ENSEMBL = character(), SYMBOL = character())
    })
  }

  mapped_rows <- c(native_rows, match(symbol_map$ENSEMBL, pb_de$gene))
  symbols <- c(pb_de$gene[native_rows], symbol_map$SYMBOL)
  if (length(mapped_rows) > 0) {
    fold_changes <- setNames(pb_de$log2FoldChange[mapped_rows], symbols)

    cat("\n--- KEGG Expanded (ORA) ---\n")
    adjusted_p <- pb_de$padj[mapped_rows]
    significant <- is.finite(adjusted_p) & adjusted_p < 0.05 &
      is.finite(fold_changes) & abs(fold_changes) > 0.25
    ora_fc <- fold_changes[significant]
    ora_fc <- ora_fc[!duplicated(names(ora_fc))]
    cat(sprintf("  ORA: %d mapped symbols with padj < 0.05 and |log2FC| > 0.25\n",
                length(ora_fc)))
    run_kegg_expanded(names(ora_fc), ora_fc, ct, table_dir, plot_dir)

    cat("\n--- GSEA with fgsea ---\n")
    pb_ranked <- fold_changes[is.finite(fold_changes)]
    pb_ranked <- pb_ranked[!duplicated(names(pb_ranked))]
    pb_ranked <- sort(pb_ranked, decreasing = TRUE)
    cat(sprintf("  Ranked %d genes (from pseudobulk log2FC)\n", length(pb_ranked)))
    if (is.null(msigdb_sets)) {
      msigdb_sets <- list(
        Hallmark = read_msigdb_collection("Hallmark", "H"),
        KEGG_MEDICUS = read_msigdb_collection("KEGG_MEDICUS", "C2", "CP:KEGG_MEDICUS"),
        GO_BP = read_msigdb_collection("GO_BP", "C5", "GO:BP")
      )
    }
    run_gsea_fgsea(pb_ranked, ct, table_dir, plot_dir, msigdb_sets)

    # ==== 3. GSEA with clusterProfiler ====
    cat("\n--- GSEA with clusterProfiler ---\n")
    entrez_map <- tryCatch({
      bitr(names(pb_ranked), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) NULL)

    if (!is.null(entrez_map) && nrow(entrez_map) > 0) {
      entrez_ranked <- pb_ranked[entrez_map$SYMBOL]
      names(entrez_ranked) <- entrez_map$ENTREZID
      entrez_ranked <- entrez_ranked[!duplicated(names(entrez_ranked))]
      entrez_ranked <- sort(entrez_ranked, decreasing = TRUE)
      cat(sprintf("  Converted %d genes to ENTREZ IDs\n", length(entrez_ranked)))
      run_gsea_clusterprofiler(entrez_ranked, ct, table_dir, plot_dir)
    } else {
      cat("  Skipping clusterProfiler GSEA: ENTREZ ID conversion failed\n")
    }
  } else {
    cat("  Skipping enrichment: gene ID mapping failed\n")
  }

  cat(sprintf("\n%s complete!\n", ct))
}

cat("\n\n========================================\n")
cat("Expanded KEGG & GSEA Analysis COMPLETE!\n")
cat("Results saved to:", RESULT_DIR, "\n")
cat("========================================\n")