#!/usr/bin/env Rscript

library(ggplot2)
library(ggpointdensity)
library(viridis)

args = commandArgs(TRUE)

d = read.csv(args[1], sep='\t')

makeplot = function(d, xstring, ystring, xlab, ylab){
    p = ggplot(d) +
        geom_pointdensity(aes_string(xstring, ystring), size = 0.5, shape = 16) +
        scale_color_viridis() +
        facet_grid(library~out) +
        labs(y = ylab, x = xlab) +
        theme_bw()

    return(p)
}

png(args[2], height=12, width=4, units="in", res=150)
print(makeplot(d, 'log_hqaa', 'log_mito', 'log10(HQAA)', 'log10(Mitochondrial fraction)'))
dev.off()

png(args[3], height=12, width=4, units="in", res=150)
print(makeplot(d, 'log_hqaa', 'log_tss', 'log10(HQAA)', 'log10(TSS enrichment)'))
dev.off()

