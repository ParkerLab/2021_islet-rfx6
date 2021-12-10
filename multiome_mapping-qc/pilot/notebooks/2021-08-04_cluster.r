#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-c", "--components"), type = "integer", action = "store", help = "Number of principal components to use"),
  make_option(c("-r", "--res"), type = "double", action = "store", dest = "res", help = "Resolution to use"),
  make_option(c("-o", "--output"), action = "store", dest = "output", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$components) || is.null(opt$res) || is.null(opt$output)) {
  stop("arguments misspecified. please see usage.")
}

#print(opt)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(rhdf5)
  library(dplyr)
  library(glue)
})

qc_df <- read.table("/lab/work/vivekrai/2021-06_pilot-rfx6/work/nuclei-qc/metrics_allcols-barcodes_post-qc.txt",
                    sep = " ", header = T)

protein_coding_genes <- rtracklayer::readGFF("/lab/work/vivekrai/2021-06_pilot-rfx6/data/reference/hg19-mCherry-mKate2_protein-coding_genes.gtf") %>%
    filter(!seqid %in% c("chrY", "chrM")) %>%
    mutate(seqid = forcats::fct_drop(seqid))

RNA_HDFS <- Sys.glob("/lab/work/vivekrai/2021-06_pilot-rfx6/work/nuclei-qc/rna/*.hdf5")

load_hdf <- function(HDF) {
    tmp <- h5ls(HDF)
    GROUP <- paste("/", tmp$name[tmp$group=="/"], sep="")
    df <- h5read(HDF, GROUP, native = F, compoundAsDataFrame=T)
    counts <- t(df$block0_values)
    rownames(counts) <- as.character(df$axis1)
    colnames(counts) <- as.character(df$axis0)
    counts <- counts[, colnames(counts) %in% protein_coding_genes$gene_name]
    print("INS" %in% colnames(counts))
    return(counts)
}

rna_objects <- list()

message("1/n Reading data..")
for (f in RNA_HDFS) {
    rna_counts <- load_hdf(f)

    library <- unlist(lapply(strsplit(rownames(rna_counts), '-'), function(x){x[1]}))

    metadata <- data.frame(library = library, nucleus = rownames(rna_counts))
    rownames(metadata) <- metadata$nucleus

    rna <- CreateSeuratObject(
        counts = t(rna_counts),
        min.cells = 1,
        assay = "RNA",
        project = as.character(unique(metadata$library)[1]),
        metadata = metadata
    )
    rna$tech <- "rna"
    rna$library <- metadata$library
    rna_objects[[length(rna_objects) + 1]] <- rna
}

if (length(rna_objects) > 1) {
    rna_additional <- c()
    for (i in 2:length(rna_objects)) {
     rna_additional <- c(rna_additional, rna_objects[[i]])
    }
    rna <- merge(rna_objects[[1]], y = rna_additional, project = "RNA")
} else {
    rna <- rna_objects[[1]]
}

message("2/n Processing..")

rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures = 2000)
rna <- ScaleData(rna)
rna <- RunPCA(rna, npcs = 100, verbose = F)

dir.create(opt$output, showWarnings = F, recursive = T)

out_path <- function(name) {
  paste0(opt$output, "/", name)
}

pdf(out_path("elbow-plot.pdf"), height = 5, width = 5)
ElbowPlot(rna, ndims = 50)
dev.off()

pdf(out_path("pc-heatmap.pdf"), height = 10, width = 10)
PCHeatmap(object = rna, dims = 1:6)
dev.off()

message("2/n UMAP..and clustering..")

rna <- FindNeighbors(rna, dims = 1:opt$components, k.param = 20)
rna <- RunUMAP(rna, reduction = 'pca', dims = 1:opt$components, verbose = F)

pdf(out_path("pc-feature-plot.pdf"), height = 10, width = 10)
FeaturePlot(rna, c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6"))
dev.off()

rna <- FindClusters(rna, resolution = opt$res, n.start = 100)

pdf(out_path("clusters.pdf"), height = 6, width = 7)
DimPlot(rna, reduction = "umap", label = TRUE)
dev.off()

n_libraries <- length(unique(rna$library))

pdf(out_path("clusters-by-library.pdf"), height = 6, width = 6 * n_libraries)
DimPlot(rna, reduction = "umap", label = TRUE, split.by = "library")
dev.off()

pdf(out_path("marker-plots.pdf"), height = 10, width = 10)
FeaturePlot(rna, c("INS-IGF2", "GCG", "SST", "RFX6", "IAPP"))
dev.off()


pdf(out_path("rfx6-by-library.pdf"), height = 6, width = 6 * n_libraries)
FeaturePlot(rna, c("RFX6"), split.by = "library")
dev.off()


pdf(out_path("marker-violin.pdf"), height = 6, width = 6)
VlnPlot(rna, c("INS-IGF2", "RFX6"))
dev.off()

message("3/n Markers and writing tables..")

cluster_markers <- FindAllMarkers(rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.table(cluster_markers,
            file = out_path("cluster-markers.txt"),
            append = F, quote = F, sep = "\t", row.names  = F, col.names = T)

clusters <- as.data.frame(rna@active.ident)
colnames(clusters) <- c("cluster")
clusters$nucleus <- rownames(clusters)
clusters$barcode <- gsub(".*-(.*)", "\\1", clusters$nucleus)
write.table(clusters[, c("nucleus", "cluster")], file = out_path("clusters.txt"),
            append = F, quote = F, sep = "\t", row.names = F, col.names = F)
