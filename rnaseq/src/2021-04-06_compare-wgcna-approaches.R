###### Approach 1: Comparing different runs using mergeCutHeight

pull_cutoff_data <- function(filename) {
  lapply(celltypes, function(x) {
    print(x)
    lapply(c("0.2", "0.3", "0.15", "0.25"), function(height) {
      print(glue("work/wgcna-explore/2021-04-07/{x}/power-80_cut-{height}/{filename}.rds"))
      df <- readRDS(glue("work/wgcna-explore/2021-04-07/{x}/power-80_cut-{height}/{filename}.rds"))
      if (!class(df) %in% c("data.frame")) {
        df <- as.data.frame(df)
      }
      df$Type <- x
      df$Cutoff <- as.numeric(height)
      df
    }) %>% bind_rows()
  }
  ) %>% bind_rows()
}

mod_assign_df <- pull_cutoff_data("module-assignment")

ggplot(mod_assign_df %>%
         group_by(Type, Cutoff) %>%
         dplyr::summarize(Modules = length(unique(module))),
      aes(Cutoff, Modules, group = Type)) +
  geom_line(aes(col = Type), size = 1.5) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = Modules)) +
  labs(title = "Number of modules vs cutoff for each celltype")


tmp <- recutBlockwiseTrees(
  celltype_norm, goodSamples = rep(T, nrow(celltype_norm)),
  goodGenes = rep(T, ncol(celltype_norm)),
  blocks = rep(1, ncol(celltype_norm)),
  TOMFiles = "work/wgcna-explore/2021-04-07/Beta.cell.RNA/power-80_cut-0.15/blockwiseTOM-block.1.RData",
  corType = "bicor",
  networkType = "signed hybrid"
)

rfx6_modules <- mod_assign_df %>% filter(symbol == "RFX6")
rfx6_modules

filter_for_rfx6_mod <- function(df, rfx6_modules) {
  df %>% rowwise() %>%
  dplyr::filter(any(
    module == rfx6_modules$module & Type == rfx6_modules$Type & Cutoff == rfx6_modules$Cutoff)
  )
}

rfx6_df <- filter_for_rfx6_mod(mod_assign_df, rfx6_modules)

ggplot(rfx6_df %>%
         group_by(Type, Cutoff) %>%
         dplyr::summarize(Genes = n()),
      aes(Cutoff, Genes, group = Type)) +
  geom_line(aes(col = Type), size = 1.5) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = Genes)) +
  labs(title = "Number of Genes in RFX6 modules per cutoff for each celltype")


kegg_assign_df <- pull_cutoff_data("mod_KEGG-summary")  %>% rowwise() %>%
  dplyr::filter(any(
    Cluster == rfx6_modules$module & Type == rfx6_modules$Type & Cutoff == rfx6_modules$Cutoff)
  ) %>%
  filter(p.adjust < 0.01)

go_assign_df <- pull_cutoff_data("mod_GOBP-summary")  %>% rowwise() %>%
  dplyr::filter(any(
    Cluster == rfx6_modules$module & Type == rfx6_modules$Type & Cutoff == rfx6_modules$Cutoff)
  ) %>%
  filter(p.adjust < 0.01)

ggplot(go_assign_df %>% filter(Type == "Beta.cell.RNA"), aes(Description, -log10(pvalue))) +
  ggalt::geom_lollipop(size = 1.25, aes(col = -log10(pvalue))) + coord_flip() +
  facet_grid(Type ~ Cutoff, scales = "free") +
  ggpubr::grids("y") +
  scale_y_continuous(expand = c(0, 0.1))

ggplot(kegg_assign_df, aes(Description, -log10(pvalue))) +
  ggalt::geom_lollipop(size = 1.25, aes(col = -log10(pvalue))) + coord_flip() +
  facet_grid(Type ~ Cutoff, scales = "free") +
  ggpubr::grids("y") +
  scale_y_continuous(expand = c(0, 0.1))


connectivity_df <- pull_cutoff_data("connectivity") %>%
  group_by(Type, Cutoff) %>%
  mutate(Rank_kTotal = rank(-kTotal), Rank_kWithin = rank(-kWithin))

connectivity_df %>% filter(symbol == "RFX6") %>%
  dplyr::select(Type, Cutoff, Rank_kTotal, Rank_kWithin) %>%
  kable()


###### Approach 2: Exploring modules using mergeCloseModules

mod_merge_explore <- lapply(celltypes, function(x) {
  print(x)
  exprData <- readRDS(glue("work/wgcna-explore/2021-04-07/{x}/power-80_cut-0.15/celltype-norm-counts.rds"))
  net <- readRDS(glue("work/wgcna-explore/2021-04-07/{x}/power-80_cut-0.15/blockwiseModules.rds"))
  print(length(unique(net$colors)))
  res <- lapply(seq(.1, .5, by = .05), function(cut) {
    print(cut)
    WGCNA::mergeCloseModules(
      exprData,
      net$colors,
      MEs = net$MEs,
      corFnc = "bicor",
      cutHeight = cut,
      relabel = TRUE
    )
  })
  names(res) <- seq(.1, .5, by = 0.05)
  res
}); names(mod_merge_explore) <- celltypes

mod_merge_summary <- lapply(celltypes, function(x) {
  lapply(names(mod_merge_explore[[x]]), function(cut) {
    data.frame(num_modules = length(unique(mod_merge_explore[[x]][[cut]]$colors)),
               rfx6_gene_count = sum(
                 mod_merge_explore[[x]][[cut]]$colors == mod_merge_explore[[x]][[cut]]$colors[which(names(mod_merge_explore[[x]][[cut]]$colors) == genes_of_interest)]
               ),
               cutHeight = cut, Type = x)
  }) %>% bind_rows()
}) %>% bind_rows()

ggplot(mod_merge_summary, aes(cutHeight, num_modules, col = Type)) +
  geom_line(aes(group = Type)) +
  geom_label(aes(label = num_modules)) +
ggplot(mod_merge_summary, aes(cutHeight, rfx6_gene_count, col = Type)) +
  geom_line(aes(group = Type)) +
  ggrepel::geom_label_repel(aes(label = rfx6_gene_count))


###### Approach 3: Comparing modules using the WGCNA parameter exploration

create_wgcna_explore_df <- function(file) {
  df <- readRDS(file)
  as.data.frame(t(do.call(cbind, lapply(df, function(x) {
    t(as.data.frame(c(
      x[1:7],
      num_modules = length(unique(x$colors)),
      rfx6_mod_size = sum(x$colors == x$colors[which(names(x$colors) == genes_of_interest)]))))
  }))))
}

beta_explore <- create_wgcna_explore_df("work/wgcna-explore/2021-04-07/Beta.cell.RNA/wgcnatest.rds")
alpha_explore <- create_wgcna_explore_df("work/wgcna-explore/2021-04-07/Alpha.cell.RNA/wgcnatest.rds")
islet_explore <- create_wgcna_explore_df("work/wgcna-explore/2021-04-07/Whole.Islet.RNA/wgcnatest.rds")


ggplot(islet_explore, aes(mergeCutHeight, num_modules)) +
  geom_line() +
  geom_label(aes(label = num_modules)) +
  facet_grid(deepSplit ~ minModSize) +
  ggpubr::grids()

ggplot(beta_explore, aes(mergeCutHeight, rfx6_mod_size)) +
  geom_line() +
  geom_label(aes(label = rfx6_mod_size)) +
  facet_grid(deepSplit ~ minModSize) +
  ggpubr::grids()


##### Compare GWAS enrichment results

## Check what is the number of bed regions in different annotations

gwas_res_path <- "~/analyses/2020-01_vanderbilt_rna/work/module-peaks_garfield/2021-04-09"

assign_significance <- function(x, padj) {
  if (x >= 0.05) {
    return("Not significant")
  } else if ((x < 0.05) & (x >= padj)) {
    return("Nominally significant (P<0.05)")
  } else if (x < padj) {
    return("Significant (Bonferroni corrected\nfor effective Annots)")
  }
}

gwas_df <- lapply(
  list.files(gwas_res_path), function(comb) {
    tss_defs <- list.files(glue("{gwas_res_path}/{comb}"))
    lapply(tss_defs, function(tss) {
      trait_enrich_files <- list.files(
        glue("{gwas_res_path}/{comb}/{tss}/results"),
        "enrichment",
        full.names = T
      )

      df <- lapply(trait_enrich_files, function(x) {
        df <- read_delim(x, delim = " ", col_types = cols())
        trait <- str_match(x, ".*\\.enrichment\\.(.*)\\.out")[2]
        mdf <- read_delim(paste0(dirname(x), "/garfield.Meff.", trait, ".out"), delim = "\t", col_types = cols(), col_names = F)

        df$Trait <- trait
        df$Significance <- sapply(df$Pvalue, function(x) assign_significance(x, mdf[2, 2]))
        df
      }) %>% bind_rows()

      df$tss_def <- tss
      df$wgcna_type <- comb
      df
    }) %>% bind_rows()
  }
) %>% bind_rows()

gwas_df$celltype <- ifelse(
  str_split(gwas_df$wgcna_type, "-", simplify = T)[, 1] == "alpha", "Alpha.cell.RNA", "Beta.cell.RNA"
)
gwas_df$merge_cut_height <- str_split(gwas_df$wgcna_type, "-", simplify = T)[, 2]
gwas_df$Annotation <- as.numeric(str_match(gwas_df$Annotation, "\\d+"))
gwas_df$tss_def <- factor(
  gwas_df$tss_def,
  levels = c("hg19.tss.5u1d_ext", "hg19.tss.1kb_ext", "hg19.tss.5kb_ext", "hg19.tss.10kb_ext"),
  ordered = T
)

plot_gwas_grid_all_tss <- function(df, merge_cut, celltype, pthresh, highlight_anno, alt.colorscheme = F) {
  df_subset <- df %>%
    filter(PThresh == pthresh) %>%
    filter(celltype == celltype) %>%
    filter(merge_cut_height == merge_cut) #%>%
    # filter(Annotation != 0)

  half_split <- length(unique(df_subset$Annotation))
  half_split <- ifelse(half_split %% 2, (half_split + 1) / 2, half_split / 2)

  p <- ggplot(df_subset, aes(tss_def, Trait, shape = Significance)) +
      geom_tile(aes(fill = Beta/SE)) +
      geom_point() +
      scale_shape_manual(values = c(96, 96, 8)) +
      ggpubr::rotate_x_text() +
      geom_rect(data = subset(df_subset, Annotation == highlight_anno),
                fill = NA, colour = "red", xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf) +
      facet_wrap(~ Annotation, ncol = half_split) +
    labs(x = "", y = "")

  if (!alt.colorscheme) {
    p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
  }
  p
}

plot_gwas_grid_all_tss(gwas_df, 0.15, "beta", 5e-8, 13)
plot_gwas_grid_all_tss(gwas_df, 0.15, "beta", 5e-8, 13, alt.colorscheme = T)
plot_gwas_grid_all_tss(gwas_df, 0.20, "beta", 5e-8, 1)
plot_gwas_grid_all_tss(gwas_df, 0.20, "beta", 5e-8, 1, T)
plot_gwas_grid_all_tss(gwas_df, 0.25, "beta", 5e-8, 2)
plot_gwas_grid_all_tss(gwas_df, 0.25, "beta", 5e-8, 2, T)
plot_gwas_grid_all_tss(gwas_df, 0.20, "alpha", 5e-8, 9)
plot_gwas_grid_all_tss(gwas_df, 0.3, "beta", 5e-8, 2)

rfx6_gwas_subset <- gwas_df %>% rowwise() %>%
  dplyr::filter(any(
    Annotation == rfx6_modules$module & celltype == rfx6_modules$Type & merge_cut_height == rfx6_modules$Cutoff)
  ) %>%
  filter(PThresh == 5e-8)

ggplot(rfx6_gwas_subset, aes(tss_def, Trait, shape = Significance)) +
  geom_tile(aes(fill = Beta/SE)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  geom_point() +
  scale_shape_manual(values = c(96, 96, 8)) +
  facet_grid(merge_cut_height ~ celltype) +
  labs(x = "", y = "") +
  ggpubr::theme_pubr(base_size = 11, legend = "right") +
  ggpubr::rotate_x_text(30)
