library(hash)
library(BiocSingular)

setwd("~/analyses/2020-01_vanderbilt_rna/")

counts <- readRDS("work/2020-02-14_qc_summary/counts.rds")
metadata <- read.table("sample_info/covariates.txt", header = T, sep = "\t")
gtf <- readRDS('data/gtf/gencode.v19.annotation.gtf.rds')

get_variances <- function(pca_fit) {
  pca_fit$sdev^2/sum(pca_fit$sdev^2)
}

get_col <- function(df, col) {
  df %>% select(!!col) %>% .[[1]]
}

get_pca <- function(counts_mat, rank = 20, center = T, ...) {
  # Transpose if Genes x Samples given
  if (nrow(counts_mat) > ncol(counts_mat)) {
    print("info: transposing matrix")
    counts_mat <- t(counts_mat)
  }
  # Samples x Genes
  result <- runPCA(counts_mat, rank = rank, center = center, ...)
  rownames(result$x) <- rownames(counts_mat)
  result
}

plot.pca <- function(pca, x = "PC1", y = "PC2", metadata = NULL, color_by = NULL) {
  # Return if rownames is not set on PCA object
  if(is.null(rownames(pca$x))) return()
  if(!"Sample.ID" %in% names(metadata)) return()

  variances <- round(get_variances(pca)*100, 1)

  df <- pca$x %>% data.frame() %>% rownames_to_column(var = "Sample.ID") %>%
    left_join(metadata %>% select(-starts_with("PC")), by = c("Sample.ID"))

  bp <- ggplot(df, aes_string(x, y))

  if(!(is.null(color_by) || is.null(metadata))) {
    bp <- bp + geom_point(aes(col = !!friendlyeval::treat_string_as_expr(color_by)))
  } else {
    bp <- bp + geom_point()
  }

  bp + labs(x = glue("{x} ({variances[as.integer(gsub('PC', '', x))]} %)"),
            y = glue("{y} ({variances[as.integer(gsub('PC', '', y))]} %)"),
            col = '')
}

#################### DATA

protein_coding_genes <- as.data.frame(gtf[gtf$type == 'gene' & gtf$gene_type == 'protein_coding'])[, c(1:5, 10, 14)] # 20345 protein-coding genes

pc_df <- counts[which(protein_coding_genes$gene_id %in% counts$Geneid), ]
pc_tpm <- scater::calculateTPM(pc_df %>% select(-Geneid, -Length) %>% as.matrix(), pc_df$Length)

# all non-mitochondrial genes --> 57783
genes <- as.data.frame(gtf[gtf$type == 'gene' & seqnames(gtf) != "chrM"])[, c(1:5, 10, 12, 14)]
genes_df <- counts[which(genes$gene_id %in% counts$Geneid), ]
genes_tpm <- scater::calculateTPM(genes_df %>% select(-Geneid, -Length) %>% as.matrix(), genes_df$Length)

geneid_to_name <- hash::hash(gtf[gtf$type == 'gene']$gene_id, gtf[gtf$type == 'gene']$gene_name)

### 3. Normalize using edgeR
normalize_edger <- function(fc_counts) {
  dge <- edgeR::DGEList(fc_counts[, -c(1,2)], lib.size = NULL, norm.factors = NULL)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = T)
}

### 3. Normalize using DESea2
normalize_deseq <- function(fc_counts, sample_info) {
  samples <- names(fc_counts[, -c(1,2)])

  colData <- metadata %>% select(Sample.ID, Type, Disease.Status) %>%
    mutate(Disease.Status = ifelse(grepl("Normal|normal", Disease.Status), "Normal", "T2D"),
           Sample.ID = factor(Sample.ID, levels = samples, ordered = T)) %>%
    arrange(Sample.ID) %>%
    column_to_rownames(var = "Sample.ID")

  stopifnot(rownames(colData) == samples)

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = fc_counts[, -c(1,2)],
                                colData = colData,
                                design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  DESeq2::getVarianceStabilizedData(dds)
}

plot.pca.grid <- function(pca_obj, k = 4, metadata = NULL, color_by = NULL) {
  plots <- lapply(seq(1, k, 2), function(i) return(plot.pca(pca_obj, x = glue("PC{i}"), y = glue("PC{i+1}"), metadata = metadata, color_by = color_by)))
  patchwork::wrap_plots(plots) + patchwork::plot_layout(guides = 'collect')
}

##################### Explore Principal Components
### 1. Bulk PCA
pca1 <- get_pca(log2(pc_tpm + 1))
saveRDS(file = "work/2020-02-14_qc_summary/pca_all-protein-coding-tpm.rds", pca1)

pca2 <- get_pca(normalize_edger(pc_df))
saveRDS(file = "work/2020-02-14_qc_summary/pca_all-protein-coding-tmm.rds", pca2)

saveRDS(file = "work/2020-02-14_qc_summary/pca_islet-protein-coding-tmm.rds",
        get_pca(normalize_edger(pc_df)[, metadata[metadata$Type == "Whole.Islet.RNA", ]$Sample.ID]))

saveRDS(file = "work/2020-02-14_qc_summary/pca_beta-protein-coding-tmm.rds",
        get_pca(normalize_edger(pc_df)[, metadata[metadata$Type == "Islet.RNA"]$Sample.ID]))

pdf("figures/04-30_figs/pre-qc/pca_all-tpm.pdf", width = 14, height = 4)
plot.pca.grid(pca1, k = 6, metadata, color_by = "Type")
dev.off()

pdf("figures/04-30_figs/pre-qc/pca_all-tmm.pdf", width = 14, height = 4)
plot.pca.grid(pca2, k = 6, metadata = metadata, color_by = "Type")
dev.off()

exclude_samples <- c("MG-0066_i-ND-15m", "MG-0124_b-ND-3", "MG-0054_a-ND-11", "MG-0055_b-ND-11",
                     "MG-0045_a-T2D-51", "MG-0104_a-ND-19", "MG-0008_b-ND-10m")

### 2. Cell-type specific PCA

is.almost.ok <- function(vec, tol = 1) {
  (sum(is.na(vec)) <= tol) & (length(unique(vec)) > 1)
}

make_celltype_plots <- function(set_of_interest, label) {
  samples_of_interest <- set_of_interest$Sample.ID
  set_of_interest <- set_of_interest %>% select(-starts_with("PC"))

  message("running pca..")
  pca_obj <- get_pca(normalize_edger(pc_df)[, samples_of_interest])
  df_with_pc <- cbind(set_of_interest, pca_obj$x[, 1:4])

  message("making heatmap..")
  ht <- df_with_pc %>%
    select(Age, Sex, Disease.Status, BMI, READ_PAIR_OK, RIN, TIN_mean, sry_exp, xist_exp, Islet.Isolation.Center,
           chrY_mean_exp, chrM_mean_exp, Sequencing.Batch, PC1, PC2, PC3, PC4, Islet.Days.Since.First) %>%
    mutate(Disease.Status = ifelse(grepl("Normal", Disease.Status), 0L, 1L)) %>%
    dummify() %>%
    select_if(is.almost.ok, ~.x) %>%
    cor(., use = "na.or.complete", method = "kendall") %>%
    Heatmap(name = label, border = T, show_column_dend = F, cluster_rows = T, cluster_columns = T,
            row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10))

  list(
    df = df_with_pc,
    pca_plot = plot.pca.grid(pca_obj, k = 6, df_with_pc, color_by = "Sequencing.Batch"),
    heat_plot = ht
  )
}

pdf("figures/pca_islet-tmm.pdf", width = 14, height = 5)
islet_res <- make_celltype_plots(metadata %>% filter(Type == "Whole.Islet.RNA", QC_OK), "Whole.Islet")
dev.off()

saveRDS(file = "work/2020-03-25_post-qc-covariates/islet_info.rds", islet_res)

pdf("figures/pca_beta-tmm.pdf", width = 14, height = 5)
beta_res <- make_celltype_pca_plots(metadata %>% filter(Type == "Beta.cell.RNA", QC_OK), "Beta")
dev.off()

saveRDS(file = "work/2020-03-25_post-qc-covariates/beta_info.rds", beta_res)

pdf("figures/pca_alpha-tmm.pdf", width = 14, height = 5)
alpha_res <- make_celltype_pca_plots(metadata %>% filter(Type == "Alpha.cell.RNA", QC_OK), "Alpha")
dev.off()

saveRDS(file = "work/2020-03-25_post-qc-covariates/alpha_info.rds", alpha_res)

## Make PC - vs - all scatter plots
vars_to_correlate <- c("Age", "BMI")
make_grid <- function(df, pc = "PC1", vars = vars_to_correlate, point.col = NULL) {
  plot_list <- lapply(vars, function(x) {
    if(is.character(tmp[[x]])) {
      ggplot(df, aes_string(x, pc)) + geom_half_boxplot(notch = F, outlier.colour = 'white') +
        geom_quasirandom(aes_string(col = point.col), size = .5) + #geom_smooth(method = "lm", col = 'black', aes(group=1)) +
        guides(col = F) + coord_flip()
    } else {
      ggplot(df, aes_string(pc, x)) + geom_point(aes_string(col = point.col), size = .5) + geom_smooth(method = "lm", col = 'black') +
        ggpubr::stat_cor(label.y.npc = "top", label.x.npc = "left") + guides(col = F)
    }
  })
  plot_grid(plotlist = plot_list)
}

select_PC <- "PC2"
make_grid(tmp, pc = select_PC)

save_plot(glue("./figures/qc/whole_{select_PC}-vs-all"), base_width = 14, base_height = 10)