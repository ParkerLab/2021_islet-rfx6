#' This script compares the WGCNA results for Beta.cell for three values of k
#' to ensure that our results are consistently reproduced across
#' different choices of parameter k (power).
#'
go_module_res12 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/12-limma-TRUE/GO-summary.rds")
go_module_res13 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/13-limma-TRUE/GO-summary.rds")
go_module_res14 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/14-limma-TRUE/GO-summary.rds")

# ME1
module_res12 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/12-limma-TRUE/module-assignment.rds")
# ME4
module_res13 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/13-limma-TRUE/module-assignment.rds")
# ME3
module_res14 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/14-limma-TRUE/module-assignment.rds")

module_conn12 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/12-limma-TRUE/connectivity.rds")
module_conn13 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/13-limma-TRUE/connectivity.rds")
module_conn14 <- readRDS("work/wgcna/2020-12-15/Beta.cell.RNA/14-limma-TRUE/connectivity.rds")

go_terms_of_interest <- c(
  "HPA:0350471",
  "HPA:0350000",
  "HPA:0350472",
  "HPA:035481",
  "KEGG:04950",
  "KEGG:01230",
  "KEGG:01100",
  "REAC:R-HSA-210745",
  "REAC:R-HSA-186712",
  "KEGG:01230",
  "HPA:0350481"
)

tmp_df <- merge(
  go_module_res13$ENSG00000185002$result,
  go_module_res14$ENSG00000185002$result,
  by="term_id",
  all=T
)

tmp_df[is.na(tmp_df)] <- 1
tmp_df$term_name = ifelse(tmp_df$term_name.x == 1, tmp_df$term_name.y, tmp_df$term_name.x)

tmp_plot <- ggplot(tmp_df, aes(
    -log10(p_value.x),
    -log10(p_value.y),
    label = term_name)
  ) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept = 0), linetype = 'dashed', col = 'red', size = 1) +
  gghighlight(term_id %in% go_terms_of_interest, label_key = term_name) +
  labs(x = "-log10(p.adj) (k=13)", y = "-log10(p.adj) (k=14)")

tmp_plot

tmp_upset_df <- merge(
  merge(
    go_module_res12$ENSG00000185002$result %>% dplyr::select(term_name) %>% dplyr::rename(k12 = term_name ),
    go_module_res13$ENSG00000185002$result %>% dplyr::select(term_name) %>% dplyr::rename(k13 = term_name ),
    all=T
  ),
  go_module_res14$ENSG00000185002$result %>% dplyr::select(term_name) %>% dplyr::rename(k14 = term_name ),
  all=T
)


all_terms <- rbind(
  go_module_res12$ENSG00000185002$result[, c("term_id", "term_name")],
  go_module_res13$ENSG00000185002$result[, c("term_id", "term_name")],
  go_module_res14$ENSG00000185002$result[, c("term_id", "term_name")]
) %>% unique()

tmp_heat_df <- all_terms %>%
  left_join(., go_module_res12$ENSG00000185002$result %>% dplyr::select(term_id, p_value) %>% dplyr::rename(k12 = p_value)) %>%
  left_join(., go_module_res13$ENSG00000185002$result %>% dplyr::select(term_id, p_value) %>% dplyr::rename(k13 = p_value)) %>%
  left_join(., go_module_res14$ENSG00000185002$result %>% dplyr::select(term_id, p_value) %>% dplyr::rename(k14 = p_value)) %>%
  unique %>%
  dplyr::filter(grepl("KEGG|REAC|GO", term_id) | grepl("pancrea", term_name))
  #filter(term_id %in% go_terms_of_interest)

tmp_heat_df[is.na(tmp_heat_df)] <- 1

head(tmp_heat_df)
heat_mat <- -log10(as.matrix(tmp_heat_df[, -c(1,2)]))
rownames(heat_mat) <- tmp_heat_df$term_name
ha = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(
    at = which(tmp_heat_df$term_id %in% go_terms_of_interest),
    labels = tmp_heat_df$term_name[which(tmp_heat_df$term_id %in% go_terms_of_interest)]
))

ComplexHeatmap::Heatmap(
  heat_mat,
  right_annotation = ha,
  row_names_side = "left",
  row_split = str_split(tmp_heat_df$term_id, ":", simplify = T)[, 1],
  cluster_rows = F,
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 8),
  col = circlize::colorRamp2(c(0, 5, 10), colors = c("white", "red", "darkred"))
)

tmp_upset_df <- data.frame(
  gene_ids = unique(c(
    module_res12 %>% filter(module == 1) %>% pull(ensembl),
    module_res13 %>% filter(module == 4) %>% pull(ensembl),
    module_res14 %>% filter(module == 3) %>% pull(ensembl)
  ))
)

tmp_upset_df$k12 <- tmp_upset_df$gene_ids %in% dplyr::pull(filter(module_res12, module == 1), ensembl)
tmp_upset_df$k13 <- tmp_upset_df$gene_ids %in% dplyr::pull(filter(module_res13, module == 4), ensembl)
tmp_upset_df$k14 <- tmp_upset_df$gene_ids %in% dplyr::pull(filter(module_res14, module == 3), ensembl)

ComplexUpset::upset(
  tmp_upset_df,
  c("k12", "k13", "k14"),
  name = "Overlapping module genes"
)

nrow(module_conn12) - rank(module_conn12$kWithin)[which(module_conn12$symbol == "RFX6")]
nrow(module_conn13) - rank(module_conn13$kWithin)[which(module_conn13$symbol == "RFX6")]
nrow(module_conn14) - rank(module_conn14$kWithin)[which(module_conn14$symbol == "RFX6")]
