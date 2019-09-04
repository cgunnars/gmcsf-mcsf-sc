library(argparse)
library(cowplot)
library(ggplot2)
library(ggrepel)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[n]
}

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
                      color = gg_color_hue(1), alpha = 1.0) + 
           geom_point(data = subset(res, (myAUC < 0.3)),
                      color = gg_color_hue(2), alpha = 1.0) +
           #geom_hline(yintercept = 2, linetype = 'dotted') + 
           geom_text_repel(aes(x = avg_logFC, y = -log10(p_val_adj), 
                               label = ifelse((myAUC > 0.7 | myAUC < 0.3), res$Row.names, ""))) +
           ylim(NA, 400)
  
ggsave(args$o, plot=volcano)
