library(argparse)
library(cowplot)
library(ggplot2)
source('src/plotutils.R')
My_Theme <- get_default_theme()
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
           geom_point(size = 0.1) +
           geom_point(data = subset(res, (myAUC > 0.7)), size = 0.1,
                      color = gg_color_hue(1), alpha = 1.0) + 
           geom_point(data = subset(res, (myAUC < 0.3)), size = 0.1,
                      color = gg_color_hue(2), alpha = 1.0) +
           theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
           #geom_vline(xintercept = 0.25, linetype = 'dotted') + 
           #geom_vline(xintercept = -0.25, linetype = 'dotted') +
           geom_text_repel(size = (6 * 5 / 14), min.segment.length = 0.5, segment.size = 0.5 * 5 / 14, 
                           box.padding = unit(3, 'pt'),
                           force = 25, 
                           aes(x = avg_logFC, y = -log10(p_val_adj), 
                           label = ifelse((myAUC > 0.75 | myAUC < 0.25), res$Row.names, ""))) +
           My_Theme + 
           xlab(expression('Average log'[2]*' fold change')) +
           ylab(expression('-log'[10]*'(Bonferroni-corrected p value)')) +
           ylim(NA, 400)
  
ggsave(args$o, plot=volcano, height = 1.7, width = 1.7, units = 'in')
