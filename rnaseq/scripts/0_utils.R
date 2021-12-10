#' Utility functions for help with Vanderbilt RNA-seq analysis

suppressPackageStartupMessages({
  require(magrittr)
  require(ggplot2)
  require(org.Hs.eg.db)
  require(glue)
})

is_almost_ok <- function(vec, tol = .5) {
  (sum(is.na(vec)) <= tol * length(vec)) & (length(unique(vec)) > 1) & (sd(vec) > 0) & (sum(vec != 0) > 1)
}

filter_eval <- function(df, filter_string) {
  dplyr::filter(df, !!friendlyeval::treat_string_as_expr(filter_string))
}

wrap_plot_multi <- function(objs, FUN, ...) {
  lapply(objs, function(x) FUN(x, ...))
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `ep_str_wrap`
#' @noRd
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

custom.dotplot <- function(object, x= "Cluster", colorBy="p.adjust",
                           showCategory=5, by="geneRatio", size="geneRatio",
                           split=NULL, includeAll=TRUE,
                           font.size=12, title="", label_format = 30,
                           group = FALSE, shape = FALSE) {
  color <- NULL
  df <- fortify(object, showCategory=showCategory, by=by,
                includeAll=includeAll, split=split)
  if (by != "geneRatio")
    df$GeneRatio <- parse_ratio(df$GeneRatio)
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  if (is.null(size)) size <- by
  by2 <- switch(size, rowPercentage = "Percentage",
                count         = "Count",
                geneRatio     = "GeneRatio")
  p <- ggplot(df, aes_string(x = x, y = "Description", size = by2))
  if (group) {
    p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size=.3) +
      ggnewscale::new_scale_colour()
  }

  if (shape) {
    ggstar <- "ggstar"
    require(ggstar, character.only=TRUE)
    # p <- p + ggsymbol::geom_symbol(aes_string(symbol = "Cluster", fill = colorBy)) +
    p <- p + ggstar::geom_star(aes_string(starshape="Cluster", fill=colorBy)) +
      scale_fill_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))
  }  else {
    p <- p +  geom_point(aes_string(color = colorBy))
  }
  suppressMessages(print(
    p + scale_color_continuous(guide=guide_colorbar(reverse=TRUE)) +
      ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
      scale_size_continuous(range=c(3, 8)) +
      scale_y_discrete(labels = label_func)
  ))
}

#' Get ENSEMBL -> ENTREZ -> GENE ID / GENE NAME mappings using annotables
get_mapping_annotables <- function(keys, keytype = "ensgene", columns = c("symbol")) {
  unique(annotables::grch37[
    annotables::grch37[[keytype]] %in% keys,
    c(keytype, columns)
  ])
}

#' Detecting elbow point from density plot.
#'
find_elbow <- function(x, y){
  n <- length(x)
  lineVec = c(x[n]-x[1], y[n]-y[1])
  lineVecNorm = lineVec/(sqrt(sum(lineVec^2)))
  vecFromFirst = cbind(x-x[1], y-y[1])
  scalaProd =rowSums(vecFromFirst * cbind(rep(lineVecNorm[1], n), rep(lineVecNorm[2], n)))
  vecFromFirstParallel = outer(scalaProd, lineVecNorm)
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2))
  idx = which.max(distToLine)
  return(x[idx])
}


many_de_heatmap <- function(de_objects, fdr = 0.05) {
  print(glue("Using FDR threshold of: {fdr}"))
  tmp <- do.call(rbind, lapply(de_objects, function(x) {
    tmp <- x$res[which(x$res$padj < fdr), ]
    data.frame(design = x$design, gene = rownames(tmp), pvalue = tmp$pvalue, log2FC = tmp$log2FoldChange)
  }))
  rownames(tmp) <- c()
  return(tmp)
}

# get_multi_mapping <- function(keys, keytype = "ENSEMBL", columns = c("SYMBOL")) {
#   suppressMessages({
#     AnnotationDbi::select(org.Hs.eg.db, keys = keys, keytype = keytype, columns = columns)
#   })
# }


get_mapping <- function(keys, keytype = "ENSEMBL", column = "SYMBOL") {
  suppressMessages({
    df <- stack(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = keys, keytype = keytype, column = column, multiVals = "first")
    )
  })
  colnames(df) <- c(column, keytype)
  return(df)
}


#' Normalization methods
#'
#' EdgeR's CPM normalization
normalize_edger <- function(counts) {
  dge <- edgeR::DGEList(counts, lib.size = NULL, norm.factors = NULL)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = T)
}

#' DESeq2's Variance Stabilizing Transformation
normalize_deseq <- function(counts, ...) {
  # Outputs a "log-like" transformed value
  DESeq2::varianceStabilizingTransformation(counts, ...)
}

#' TPM normalization
normalize_tpm <- function(counts, lengths = NULL) {
  if (!"Length" %in% colnames(counts)) {
    if (is.null(lengths)) {
      print("Gene lengths not specified.")
      return()
    }
  } else {
    lengths <- counts$Length
  }
  scater::calculateTPM(counts, lengths)
}

make_grid_plot <- function(pca, meta, k = 4,
                           cols.of.interest = c("Age", "BMI", "Disease.Status", "Sex", "Batch"),
                           alpha = .5, color = NULL, ...) {
  col_names <- colnames(meta)
  tmp_plots <- lapply(cols.of.interest, function(x) {
    if (!x %in% col_names) {
      message(glue::glue("{x} not found, skipping."))
      next
    }

    if (!is.null(color)) {
      if (color %in% col_names) {
        color <- celltype_meta[[color]]
      }
    } else {
      color <- NULL
    }

    lapply(paste0("PC", 1:k), function(y) {
      geom <- "auto"

      if (is.character(celltype_meta[[x]]) | is.factor(celltype_meta[[x]])) {
        geom <- c("jitter", "boxplot")
        qplot(celltype_pca$x[, y], celltype_meta[[x]], geom = geom, alpha = .5, color = color, ...)
      } else if (is.numeric(celltype_meta[[x]])) {
        geom <- c("point", "smooth")
        qplot(celltype_pca$x[, y], celltype_meta[[x]], method = "lm", geom = geom, alpha = .5, color = color, ...)
      }
    })
  })

  GGally::ggmatrix(
    purrr::flatten(tmp_plots),
    ncol = 4,
    nrow = length(cols.of.interest),
    xAxisLabels = paste0("PC", 1:k),
    yAxisLabels = cols.of.interest,
    progress = FALSE
  )
}

plot_genes_across_models <- function(obj, gene_symbols, fdr = 0.05, ...) {
  symbols <- get_mapping(rownames(obj[[1]]$res))$SYMBOL
  idx <- which(symbols %in% gene_symbols)

  res <- do.call(rbind, lapply(obj, function(x) {
    tmp_r <- x$res[idx, ]
    tmp_r$Gene <- symbols[idx]
    tmp_r$Design <- x$design
    tmp_r
  }))
  res %>% as_tibble() %>%
    ggplot(aes(Design, Gene)) + geom_tile(aes(fill = log2FoldChange)) +
    geom_point(aes(size = ifelse((padj < fdr) & (!is.na(padj)), 'Sig.', 'N.S')), shape = 21) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
    theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0)) +
    scale_size_manual(values = c(0, 4)) +
    scale_fill_distiller(palette = "RdBu") +
    labs(size = glue("FDR {fdr}"), x = '', y = '')
  # %>% tidyHeatmap::heatmap(Gene, Design, log2FoldChange, border = TRUE, ...)
}

plot_compare_foldchange_models <- function(obj) {
  df <- do.call(cbind, lapply(obj, function(x) x$res$log2FoldChange))
  plot_list <- lapply(1:ncol(df), function(i) {
    lapply(1:ncol(df), function(j) {
      if (i == j || i < j) {
        ggplot(NULL)
      } else {
        ggplot(NULL, aes(df[, i], df[, j])) +
          geom_point(size = .1) + geom_abline(slope = 1, intercept = 0)
      }
    })
  })
  GGally::ggmatrix(
    purrr::flatten(plot_list),
    ncol = ncol(df),
    nrow = ncol(df),
    xAxisLabels = colnames(df),
    yAxisLabels = colnames(df),
    progress = FALSE
  )
}

plot_count_heatmap <- function(mat, transform = NULL, name = 'matrix', n = 50, symbol = TRUE) {
  mat_ <- mat[order(rowMeans(mat), decreasing = TRUE)[1:n], ]
  labels <- rownames(mat_)
  if (symbol) {
    labels <- get_mapping(keys = labels)$SYMBOL
  }
  ComplexHeatmap::Heatmap(mat_, name = name,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 7),
    row_labels = labels
  )
}


get_variances <- function(pca_fit) {
  (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
}

#' Redundant. Use `dplyr::pull`
select_col <- function(df, col) {
  dplyr::select(df, !!col)[, 1]
}

get_pca <- function(counts_mat, rank = 10, center = T, ...) {
  # Transpose if Genes x Samples given
  if (nrow(counts_mat) > ncol(counts_mat)) {
    print("info: transposing matrix")
    counts_mat <- t(counts_mat)
  }
  # Samples x Genes
  result <- BiocSingular::runPCA(counts_mat, rank = rank)
  rownames(result$x) <- rownames(counts_mat)
  return(result)
}

plot_pca <- function(pca, x = "PC1", y = "PC2", metadata = NULL, color_by = NULL, ...) {
  # Return if rownames is not set on PCA object
  if (is.null(rownames(pca$x))) {
    print("row.names is not set on x")
    return()
  }
  if (!"Sample.ID" %in% names(metadata)) {
    print("could not find Sample.ID in metadata")
    return()
  }

  vars <- round(get_variances(pca) * 100, 1)

  df <- data.frame(pca$x) %>% tibble::rownames_to_column(var = "Sample.ID") %>%
    left_join(metadata %>% dplyr::select(-starts_with("PC")), by = c("Sample.ID"))

  bp <- ggplot(df, aes_string(x, y), ...)

  if(!(is.null(color_by) || is.null(metadata))) {
    bp <- bp + geom_point(aes(col = !!friendlyeval::treat_string_as_expr(color_by)))
  } else {
    bp <- bp + geom_point()
  }

  bp + labs(x = glue("{x} ({vars[as.integer(gsub('PC', '', x))]} %)"),
            y = glue("{y} ({vars[as.integer(gsub('PC', '', y))]} %)"),
            col = '')
}


plot_pca_grid <- function(pca, metadata = NULL,
    additional.cols = NULL, k = 4, color_by = NULL, method = "kendall", ...) {
  df <- data.frame(pca$x[, 1:k])
  vars <- round(get_variances(pca) * 100, 1)[1:k]
  colnames(df) <- sapply(1:k, function(i) glue::glue("PC{i} ({vars[i]}%)"))

  if (!(is.null(metadata) | is.null(color_by))) {
    df[[color_by]] <- metadata[match(rownames(pca$x), metadata$Sample.ID), ][[color_by]]
  }

  if (!is.null(additional.cols) & !is.null(metadata)) {
    for (i in additional.cols) {
      df[[i]] <- metadata[match(rownames(pca$x), metadata$Sample.ID), ][[i]]
    }
  }

  if (is.null(additional.cols)) {
    GGally::ggpairs(df,
      lower = list(continuous = GGally::wrap("cor", size = 6, method = method), combo = "blank"),
      upper = list(continuous = "points"),
      mapping = ggplot2::aes_string(col = color_by, ...), progress = F)
  } else {
    GGally::ggpairs(df,
      lower = list(continuous = GGally::wrap("cor", size = 6, method = method), combo = "blank"),
      upper = list(continuous = "smooth"),
      mapping = ggplot2::aes_string(col = color_by, ...), progress = F)
  }
}


get_correlation <- function(df, method = "kendall") {
  if (!"Sample.ID" %in% names(df)) {
    print("Sample.ID not found in data.")
    return()
  }
  dplyr::select(df, -any_of(c("Type", "Sample.ID", "Donor.ID"))) %>%
    DataExplorer::dummify() %>%
    dplyr::select_if(is_almost_ok, ~.x) %>%
    psych::corr.test(., use = "na.or.complete", method = method)
}

plot_correlation <- function(df, label, absolute = FALSE) {
  if (absolute) {
    df <- abs(df)
  }
  ComplexHeatmap::Heatmap(df, name = label, border = T, show_column_dend = F, cluster_rows = T, cluster_columns = T,
    row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8))
}

inv_norm <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

qqplot_custom <- function (observedPValues, main = NULL) {
  qplot(-log10(1:length(observedPValues) / length(observedPValues)),
    -log10(sort(observedPValues)),
    xlab = "Expected P-value (-log10)",
    ylab = "Observed P-value (-log10)",
    main = main) + geom_abline(data = NULL, intercept = 0, slope =  1, col = 'red') +
    theme(plot.title = element_text(size = 10))
}

combine_by_sampleid <- function(df, ...) {
  dplyr::left_join(df, ..., by = c("Sample.ID"))
}

combine_pca_with_metadata <- function(pca, ...) {
  dplyr::left_join(
    data.frame(pca$x[, 1:4]) %>% tibble::rownames_to_column("Sample.ID"),
    ...,
    by = c("Sample.ID")
  )
}

plot_counts <- function(de_obj, gene = NULL) {
  if (is.null(gene))
    gene = which.min(de_obj$res$padj)
  DESeq2::plotCounts(de_obj$dds, gene = gene, intgroup = "Disease.Status")
}

#' Create DDS object.
#' This object rarely changes, so do not compute it everytime.
create_dds_obj <- function(counts_mat, metadata) {
  sample_ids <- colnames(counts_mat)
  message(glue("info: counts matrix: {paste(dim(counts_mat), collapse = ', ')}"))
  message(glue("info: covariate matrix: {paste(dim(metadata), collapse = ', ')}"))

  col_data <- metadata %>%
    mutate(Sample.ID = factor(Sample.ID, levels = sample_ids, ordered = T)) %>%
    arrange(Sample.ID) %>%
    mutate_if(is.numeric, ~scale(.x)) %>%
    mutate_if(is.factor, ~forcats::fct_drop(.x)) %>%
    column_to_rownames(var = "Sample.ID")

  if (any(rownames(col_data) != colnames(counts_mat))) {
    print("names don't match. check or error!")
    return()
  }

  message("info: creating dds object..")
  DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~1)
}


run_many_designs_deseq <- function(dds, design, additional_covs, contrast, ...) {
  result <- list()
  for(i in 0:length(additional_covs)) {
    design_str <- paste(c(design, additional_covs[i]), collapse = "+")
    suppressMessages({
      resDe <- run_differential(
        dds,
        design = design_str,
        contrast = contrast_vec,
        ...
      )
    })
    result[[design_str]] <- resDe
  }
  return(result)
}

#' Run the differential expression logic for the given design formula
run_differential <- function(dds, design, contrast, interaction = NULL, fdr = 0.05, shrink = FALSE) {
  if (!is.null(interaction))
    design <- paste0(design, "+", interaction)

  writeLines(glue("info: running design {design}"))
  design(dds) <- as.formula(design)  # replace / update design

  dds <- DESeq2::DESeq(dds)
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))

  # extract results for the specified FDR threshold
  result <- DESeq2::results(dds, contrast = contrast, alpha = fdr)

  # if shrink is specifies, then report shrunken log2FC using apeglm method
  if (shrink) {
    writeLines("info: estimating shrunk log2FC")
    coef <- gsub(' ', '.', glue("{contrast[1]}_{contrast[2]}_vs_{contrast[3]}"))
    result <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm")
  }

  return(list(dds = dds, vst = vst, res = result, design = design, shrink = shrink))
}


run_ruvseq <- function(dds, design, contrast, k = 4, p.val.thresh = 0.5, method = "ruvg", ...) {
  # get counts from DESeq2 -- don't normalize yet
  counts <- DESeq2::counts(dds, normalized = FALSE)

  set <- EDASeq::newSeqExpressionSet(counts)
  set <- EDASeq::betweenLaneNormalization(set, which="upper")

  if (method == "ruvg") {
    # Run first-pass DESeq scan to get empirical control genes
    suppressMessages({
      de_res <- run_differential(dds, design = design, contrast = contrast, ...)
    })

    # Use non-significant genes as as empirical control genes
    not.sig <- rownames(de_res$res)[which(de_res$res$pvalue > p.val.thresh)]
    empirical <- rownames(set)[rownames(set) %in% not.sig]
  } else if (method == "ruvs") {
    empirical <- rownames(set)
  }

  print(glue("info: using {length(empirical)} control genes to estimate variation.."))

  result <- list()
  for (num_latent in 1:k) {
    # Run specified method of unwanted variation control
    if (method == "ruvg") {
      print(glue("info: estimating k={num_latent} factors using RUVg.."))
      new_set <- RUVSeq::RUVg(set, empirical, k = num_latent) # Run RUVSeq
    } else if (method == "ruvs") {
      print(glue("info: estimating k={num_latent} factors using RUVs.."))
      groups <- RUVSeq::makeGroups(SummarizedExperiment::colData(dds)[, contrast[1]])
      new_set <- RUVSeq::RUVs(set, empirical, k = num_latent, groups) # Run RUVSeq
    }

    ddsruv <- dds
    for (w in 1:num_latent) {
      col_idx <- glue("W_{w}")
      ddsruv[[col_idx]] <- new_set[[col_idx]]
    }

    new_design <- glue("{design} + {paste('W', 1:num_latent, collapse = '+', sep = '_')}")

    suppressMessages({
      de_res <- run_differential(ddsruv, design = new_design, contrast = contrast, ...)
    })

    result[[new_design]] <- list(de = de_res, W = pData(new_set), normCounts = new_set@assayData$normalizedCounts)
  }

  return(result)
}

plot_ruv_diagnostics <- function(normCounts, metadata, group = "Disease.Status", main = NULL, ...) {
  ord_idx <- match(colnames(normCounts), metadata$Sample.ID)
  colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(metadata[ord_idx, ][[group]])]
  par(mfrow=c(1,2))
  EDASeq::plotRLE(normCounts, col = colors, outline = FALSE, las = 3, cex.axis = .5, ylab = "Relative Log Expression", main = main, cex.main = .5)
  EDASeq::plotPCA(normCounts, col = colors, cex = .75, cex.lab = 1, cex.axis = 1, ...)
}

run_sva <- function(raw_counts, metadata, contrast, k = "auto",
    base_model = "~ Disease.Status + Age + Sex + BMI",
    additional.covs = c("TIN_mean", "Batch", "Islet.Isolation.Center", "Race"), ...) {

  full_model <- paste(base_model, paste0(additional.covs, collapse = "+"), sep = '+')
  metadata <- metadata %>%
    dplyr::select(all_of(
      c("Sample.ID", unlist(stringr::str_extract_all(full_model, "([A-Za-z.-_]+)", simplify = T)))
    ))

  print("Starting svaseq..")
  result <- list()

  # full model contains both adjustment variables and variables of interest (e.g., disease.status)
  # null model contains only adjustment variables
  size.factors <- DESeq2::estimateSizeFactorsForMatrix(raw_counts)
  norm_counts <- t(apply(raw_counts, 1, function(x) x/size.factors))

  sva_res <- sva::svaseq(norm_counts, mod = model.matrix(as.formula(full_model), metadata),
                         mod0 = model.matrix(as.formula(base_model), metadata))

  if (k == "auto") {
    k <- sva_res$n.sv
  }

  print("Starting DESeq..")
  design <- glue("~ {gsub('~', '', base_model)} + {paste('V', 1:k, collapse = '+', sep = '')}")

  print(design)

  sva_df <- data.frame(sva_res$sv[, 1:k])
  colnames(sva_df) <- paste0("V", 1:k)
  sva_df$Sample.ID <- colnames(counts)

  meta_df <- dplyr::left_join(metadata, sva_df, by = c("Sample.ID"))

  suppressMessages({
    de_res <- run_differential(
      raw_counts,
      meta_df,
      design = design,
      contrast = contrast,
      ...
    )
  })

  result[[design]] <- list(de = de_res, sva = sva_res)
  return(result)
}

merge_multi_go_results <- function(go_list) {
  list_names <- names(go_list)
  do.call(rbind, lapply(1:length(go_list), function(i) {
    cbind(go_list[[i]], name = list_names[i])
  }))
}

plot_multi_go_results <- function(df, fdr = FDR) {
  pathways_of_interest <- c("hsa04950", "hsa04940", "hsa04930", "hsa04911", "hsa04910", "hsa04972")
  test <- df %>% mutate(odds.ratio = odds.ratio - 1) %>% replace(., is.na(.), 0) %>%
  dplyr::select(Concept.name, Concept.ID, status, name, odds.ratio, FDR) %>%
    group_by(name, status) %>% filter((FDR < fdr) | (Concept.ID %in% pathways_of_interest))

  ggplot(test, aes(str_trunc(Concept.name, 30), name)) +
    geom_raster(aes(fill = odds.ratio)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
    theme(axis.text.y = element_text(size=8), axis.text.x = element_text(angle = 45, size = 8, hjust = 1)) +
    coord_flip() +
    geom_point(aes(size = -log10(FDR)), color = 'black', fill = 'white', shape = 21, alpha = .5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(x = '', y = 'Pathways')
}

plot_volcano <- function(de_object, fdr = 0.05, gene.symbol = TRUE) {
  xs <- data.frame(de_object$res)

  topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
  topxs$gene <- get_mapping(topxs$geneid)$SYMBOL

  plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
    scale_color_manual(values = c("gray", "firebrick")) +
    labs(col = "", title = de_object$design) +
    theme(plot.title = element_text(size = 10))

  if (gene.symbol) {
    plot <- plot +
    ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = gene), size = 3)
  } else {
    plot <- plot +
    ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)
  }

  return(plot)
}

plot_gene_counts <- function(de_obj, metadata, num.genes = 16) {
  res <- de_obj$res[order(de_obj$res$padj)[1:num.genes], ]
  counts <- log2(as.data.frame(DESeq2::counts(de_obj$dds, normalized = TRUE)[rownames(res), ]) + 1)
  counts$Geneid <- rownames(counts)

  if (nrow(res) < 1) {
    print("No genes to plot")
    return()
  }

  counts %>%
    mutate(Gene = get_mapping(Geneid)$SYMBOL) %>%
    tidyr::pivot_longer(cols = -c(Gene, Geneid), names_to = "Sample.ID", values_to = "Count") %>%
    dplyr::left_join(., dplyr::select(metadata, "Disease.Status", "Sample.ID"), by = "Sample.ID") %>%
    ggplot(aes(Disease.Status, Count, col = Disease.Status)) + geom_boxplot() +
    geom_jitter() + facet_wrap(~Gene, scales = "free") +
    theme(axis.text.x = element_blank()) +
    labs(y = "log2 abundance")
}

plot_deseq_heatmap <- function(de_object, metadata, split_by = NA, fdr = 0.05, ...) {
  sig.genes <- which(de_object$res$padj < fdr)
  norm_counts <- as.data.frame(de_object$vst[sig.genes, ])

  norm_counts$Gene <- get_mapping(rownames(norm_counts))$SYMBOL

  if (nrow(norm_counts) < 1) {
    print("No genes to plot. Adjust FDR if neccessary.")
    return()
  }

  meta_cols <- na.omit(c(split_by, "Sample.ID", "Sex"))

  norm_counts %>%
    dplyr::filter(!is.na(Gene)) %>%
    tidyr::pivot_longer(cols = -c(Gene), names_to = "Sample.ID", values_to = "Expression") %>%
    dplyr::left_join(., dplyr::select(metadata, any_of(meta_cols)), by = "Sample.ID") %>%
    { if (!is.na(split_by)) dplyr::group_by(., !!!friendlyeval::treat_strings_as_exprs(split_by)) else . } %>%
    tidyHeatmap::heatmap(
      Gene,
      Sample.ID,
      Expression,
      .scale = "row",
      palette_value = c("blue", "white", "red"),
      ...
    ) %>%
    add_tile(Sex)
}

quantify_gene_cov_associations <- function(tpm, metadata,
  covs = c("Age", "Sex", "BMI", "PC1", "PC2", "RIN", "Islet.Isolation.Center", "Disease.Status",
           "TIN_mean", "Batch", "Islet.Days.Since.First")) {
  res <- apply(tpm, 1, function(x) {
    x_transform <- inv_norm(x)
    right_side_x <- paste0(glue::glue("metadata${covs}"), collapse = "+")
    lm(formula = as.formula(glue::glue("x_transform ~ 0 + {right_side_x}")))
  })

  res <- do.call(cbind, lapply(res, function(x) {
    data.frame(summary(x)$coef[, 4])
  }))

  data.frame("pct_genes_associated" = apply(res, 1, function(x) {
    sum(x < 0.05, na.rm = T) * 100 / length(na.omit(x))
  }))
}


kegg_overrep_analysis <- function(de_res, fdr = 0.05, pvalue.cutoff = 0.05) {
  res <- de_res[which(de_res$padj < fdr), ]
  gene_list <- get_mapping(rownames(res), column = "ENTREZID")$ENTREZID
  clusterProfiler::enrichKEGG(gene = na.omit(gene_list), organism = 'hsa', pvalueCutoff = pvalue.cutoff)
}

kegg_gsea <- function(de_res, fdr = 0.05) {
  res <- de_res[which(de_res$padj < fdr), ]
  geneList <- res$log2FoldChange
  names(geneList) <- hash::values(geneid_to_entrez, rownames(res))
  geneList <- sort(geneList, decreasing = T)

  clusterProfiler::gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
}

run_go_analysis <- function(deseq_obj, database = "KEGG", directional = FALSE) {
  if (!"res" %in% names(deseq_obj)) {
    print("res item not found")
    return()
  }

  if (directional) {
    direction <- deseq_obj$res$log2FoldChange
  } else {
    direction <- NULL
  }

  rna_enrich(deseq_obj$res$pvalue,
             get_mapping(rownames(deseq_obj$res), column = "ENTREZID")$ENTREZID,
             deseq_obj$res$baseMean,
             "hsa",
             direction = direction,
             database = database)
}

format_rnaenrich_results <- function(result, n = 20) {
  result %>% dplyr::group_by(status) %>% dplyr::top_n(n, -p.value) %>% dplyr::ungroup
}

plot_go_volcano <- function(go, fdr = 0.05) {
  ggplot(go, aes(odds.ratio, -log10(p.value), size = n.genes)) +
    geom_point(aes(col = ifelse(FDR < fdr, "red", "gray"))) +
    ggrepel::geom_text_repel(aes(label = ifelse(FDR < fdr, Concept.name, NA)), size = 3) +
    labs(x = "Odds ratio", y = "-log10(P-value)", size = "Genes", col = glue("P-value (FDR < {fdr})")) +
    scale_color_manual(values  = c("gray", "red"), labels = c("N.S", "Signif."))
}

many_de_summary <- function(de_objects, fdr = 0.05, log2fc.abs = 0) {
  print(glue("Using FDR={fdr} and log2FC={log2fc.abs} to filter genes"))
  do.call(rbind, lapply(de_objects, function(x) {
    tmp <- x$res[which(x$res$padj < fdr & abs(x$res$log2FoldChange) > log2fc.abs), ]
    if (purrr::is_empty(tmp)) {
      return(
        data.frame(total_signif = 0, up = 0, down = 0)
      )
    }
    data.frame(total_signif = nrow(tmp), up = sum(tmp$log2FoldChange > 0), down = sum(tmp$log2FoldChange < 0))
  }))
}


plot_ruv_loadings <- function(obj) {
  tmp_df <- obj$W
  rownames(tmp_df) <- colnames(obj$de$vst)
  ComplexHeatmap::Heatmap(tmp_df, name = " ", cluster_columns = F,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.3f", tmp_df[i, j]), x, y, gp = grid::gpar(fontsize = 10))
  })
}
