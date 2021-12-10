#! /usr/bin/env Rscript

library(optparse)

seed <- 2020

option_list <- list(
  make_option(c("--counts"), action = "store", type = "character", default = NA, help = "[Required] Path to counts RDS file"),
  make_option(c("--outfile"), action = "store", type = "character", default = NA, help = "[Required] Output result"),
  make_option(c("--cores"), action = "store", type = "numeric", default = 2, help = "[Optional] Number of threads to use")
)

option_parser <- OptionParser(
  usage = "usage: Rscript $prog [options]",
  option_list = option_list,
  add_help_option = T
)

opts <- parse_args(option_parser)

if (any(is.na(c(opts$counts, opts$outfile)))) {
  stop("fatal: see help.")
}

suppressPackageStartupMessages({
  library(glue)
  library(WGCNA)
})

root.dir <- "/home/vivekrai/analyses/2020-01_vanderbilt_rna"

# Setup base parameters to be used throughout
network_type <- "signed hybrid"
tom_type <- "signed"
cor_type <- "bicor"

# Part of the code borrowed from the following repository.[1]
# [1]: https://github.com/liukf10/DDPNA/blob/2238eb20c6735a6e797dc7eb356a0e67d265f792/R/Moduleconstruction.R#L32
wgcnatest <- function(data, power = NULL, TOMType = "unsigned",
                      detectCutHeight = NULL, maxBlockSize = 15000,
                      deepSplit = TRUE, minModSize = TRUE,
                      pamRespectsDendro = FALSE,
                      minKMEtoStay = TRUE,
                      minCoreKME = FALSE,
                      reassignThreshold = FALSE,
                      mergeCutHeight = FALSE,
                      maxModNum = 30,
                      minModNum = 8,
                      MaxMod0ratio = 0.3,
                      corType = "bicor",
                      outPrefix = ".",
                      ...) {

  enableWGCNAThreads(opts$cores)

  if( ncol(data) < 50) stop("The number of proteins is too less to fit the function.")
  if( maxBlockSize > 45000) stop("maxBlockSize is too larger to fit the function.")
  if(!is.numeric(power) & length(power) > 1)
    stop("power is not a correct value.")
  if(!is.numeric(detectCutHeight) & length(detectCutHeight) > 1)
    stop("detectCutHeight is not a correct value.")
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop ("WGCNA package is needed. Please install it.",
          call. = FALSE)}

  cor <- WGCNA::cor;
  block <- ceiling(ncol(data) / maxBlockSize);
  .parameter <- function(x, num1, num2, name){
    if(is.null(x) || all(is.na(x)))
      x = num1 else if (isFALSE(x))
        x = num1 else if (isTRUE(x))
          x <- num2 else if (!is.numeric(x))
            stop(paste("wrong", name));
        x
  }

  deepSplit <- .parameter(deepSplit, 2, c(0,1,2,3,4), "deepSplit");
  minModSize <- .parameter(minModSize, 20, c(15, 20, 30, 50), "minModSize");
  minKMEtoStay <- .parameter(minKMEtoStay, 0.3, c(0.1, 0.2, 0.3), "minKMEtoStay");
  reassignThreshold <- .parameter(reassignThreshold, 1e-6, c(0.01, 0.05), "reassignThreshold");
  mergeCutHeight <- .parameter(mergeCutHeight, 0.15, c(0.15, 0.3, 0.45), "mergeCutHeight");
  minCoreKME <- .parameter(minCoreKME, 0.5, c(0.4, 0.5), "minCoreKME");

  if(is.null(power)) {
    message("..estimating soft thresholding power")
    sft <- WGCNA::pickSoftThreshold(
      data,
      blockSize = maxBlockSize,
      corFnc = corType,
      networkType = network_type,
      RsquaredCut = 0.80
    )
    power <- sft$powerEstimate
    message(glue("..OK. Choosing power={power} for at least 80% fit."))
  }
  #r when p<0.05 detectcutheight
  if(is.null(detectCutHeight)){
    num <- nrow(data); p = 0; r <- 0.6
    while (p < 0.05) {
      r <- r-0.01
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p > 0.05) {
      r <- r+0.001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p < 0.05) {
      r <- r-0.0001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    while (p > 0.05) {
      r <- r+0.00001
      p <- pt(-r*sqrt((num-2)/(1-r^2)),num-2)*2
    }
    detectCutHeight = 1-r^power
    if(detectCutHeight < 0.995) detectCutHeight = 0.995
    rm(r,p,num)
  }

  message("..do first WGCNA run with default params")
  net <- try(WGCNA::blockwiseModules(data, power = power, maxBlockSize = maxBlockSize,
                                     TOMType = TOMType, deepSplit = 2, minModuleSize = 20,
                                     reassignThreshold = 0.05, mergeCutHeight = 0.15,
                                     numericLabels = TRUE, pamRespectsDendro = FALSE,
                                     corType = corType, nThreads = opts$cores, networkType = network_type,
                                     saveTOMs = TRUE, loadTOM = TRUE,
                                     saveTOMFileBase = glue("{outPrefix}/blockwisetom"),
                                     verbose = 2, ...), silent = T)

  if (class(net) == "try-error") #190607
    stop("Run require(\"WGCNA\") first, if error again, please check the data.")

  message("..success. Starting exploration..")

  iter <- 1
  results <- list()

  for (i1 in deepSplit) {
    for (i2 in minModSize) {
      for (i3 in minKMEtoStay) {
        for(i4 in reassignThreshold){
          for (i5 in mergeCutHeight) {
            for (i6 in minCoreKME) {
              message(glue("Running: deepSplit {i1}, minModSize {i2}, minKME {i3}, reassign {i4}, mergeCut {i5}, minCore {i6}"))
              net <- WGCNA::blockwiseModules(data, loadTOM = TRUE, power = power,
                                             maxBlockSize = maxBlockSize, TOMType = TOMType,
                                             deepSplit = i1, minModuleSize = i2,
                                             detectCutHeight = detectCutHeight,
                                             minKMEtoStay = i3, minCoreKME = i6,
                                             reassignThreshold = i4, mergeCutHeight = i5,
                                             numericLabels = TRUE,
                                             pamRespectsDendro = pamRespectsDendro,
                                             nThreads = opts$cores,
                                             corType = corType, networkType = network_type,
                                             saveTOMFileBase = glue("{outPrefix}/blockwisetom"),
                                             saveTOMs = FALSE,
                                             verbose = 2, ...)

              iter_res <- list(
                  deepSplit = i1,
                  minModSize = i2,
                  minKMEtoStay = i3,
                  reassignThreshold = i4,
                  mergeCutHeight = i5,
                  minCoreKME = i6,
                  iter = iter,
                  colors = net$colors
              )

              results <- append(results, list(iter_res))
              iter <- iter + 1
            }
          }
        }
      }
    }
  }

  # filetoremove <- paste(glue("{outPrefix}/blockwisetom-block."), 1:block, ".RData", sep="")
  # file.remove(filetoremove)
  results
}

message("(1/2) Reading counts..")
counts <- readRDS(opts$counts)

message("(2/2) Exploring parameter space.. will take a while")
dir.create(dirname(opts$outfile), recursive = T, showWarnings = F)

WGCNAadjust <- wgcnatest(
  counts,
  TOMType = tom_type,
  deepSplit = c(0, 2, 4),
  minModSize = c(20, 30, 50),
  mergeCutHeight = c(0.15, 0.20, 0.25, 0.30, .4),
  maxBlockSize = 15500,
  minCoreKME = 0.5,
  reassignThreshold = 1e-6,
  minKMEtoStay = 0.3,
  corType = cor_type,
  outPrefix = dirname(opts$outfile))

tryCatch({
  saveRDS(WGCNAadjust, file = opts$outfile)
}, error = {
  saveRDS(WGCNAadjust, file = basename(opts$outfile))
})
