library(argparse)
library(cowplot)
library(ggplot2)
library(ggrepel)
parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='data file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
args <- parser$parse_args()

res <- readRDS(args$i)
volcano <- ggplot(res, aes(avg_logFC, -log10(p_val_adj))) + 
           geom_point() +
           geom_point(data = subset(res, (myAUC > 0.7)),
                      color = 'orange', alpha = 1.0) + 
           geom_point(data = subset(res, (myAUC < 0.3)),
                      color = 'blue', alpha = 1.0) +
           geom_hline(yintercept = 2, linetype = 'dotted') + 
           geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), 
                               label = ifelse((myAUC > 0.7 | myAUC < 0.3), res$Row.names, ""))) +
           ylim(NA, 400)
  
ggsave(args$o, plot=volcano)

#with(subset(res, padj<.05 ), points(avg_, -log10(pvalue), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
#with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
