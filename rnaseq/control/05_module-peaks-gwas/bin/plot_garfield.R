library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(glue)
options(repr.plot.width=3, repr.plot.height=6, repr.plot.res=150)

savep <- function(p, prefix, h = 10, w = 4) {
    png(glue("{prefix}.png"), height = h, width = w, units = "in", res = 300)
    print(p)
    dev.off()

    pdf(glue("{prefix}.pdf"), height = h, width = w)
    print(p)
    dev.off()
}

makePointPlot_3Pvals <- function(d) {
  p <- ggplot(d, aes(x = reorder(linkID, Beta), y=Beta)) +
    geom_point(data = d, aes(colour=significance), shape=16, size=1.5) +
    geom_hline(yintercept = 0, color="black", size=0.5)  +
    geom_errorbar(data = subset(d, Pvalue < 0.05), aes(ymin = CI95_lower, ymax = CI95_upper), size=0.2, height=0.1) +
    theme(text=element_text(size=8),
          axis.text.x=element_text(size=8),
          strip.text.x=element_text(size=7),
          panel.background = element_rect(fill="white", colour="black"),
          legend.key.size=unit(3, "mm"),
          panel.grid=element_blank(),
          legend.position="right",
          legend.text=element_text(size=7),
          panel.grid.major.y=element_line(colour="grey", size=0.2, linetype="dashed"),
          panel.grid.major.x=element_line(colour="grey", size=0.2, linetype="dashed")) +
    labs(y="Beta (log OR)", x="Annotations") +
    coord_flip() + 
    facet_wrap(~PThresh) +
    guides(col = guide_legend(nrow = 3)) +
    scale_color_manual(values=c("red", "gray", "darkred"))

    return(p)
}

assign_significance <- function(x){
  if (x >= 0.05){
    return("Not significant")
  } else if ((x < 0.05) & (x >= padj)){
    return("Nominally significant (P<0.05)")
  } else if (x < padj){
    return("Significant (Bonferroni corrected\nfor effective Annots)")
  }
}

args <- commandArgs(TRUE)

data = args[1]
prefix = args[2]
padj = as.numeric(args[3])
d <- read.csv(data, header=T, sep=' ')

d$significance <- sapply(d$Pvalue, assign_significance)

## ## SUBSET Point PLOT
lapply(unique(d$PThresh), function(i) {
    f <- glue("{prefix}.{i}")
    p <- makePointPlot_3Pvals(d[d$PThresh == i, ])
    savep(p, f, 6, 6)
})