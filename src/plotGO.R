library(ggplot2)
library(cowplot)
library(argparse)
library(forcats)
library(scales)
library(ggpubr)
library(stringr)

source('src/plotutils.R')
My_Theme <- get_default_theme()

frac2dec <- function(frac) {
  return(sapply(frac, function(x) eval(parse(text=as.character(x)))))
}

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs='+',
                    help='data file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
args <- parser$parse_args()

plots <- lapply(1:length(args$i), c)

max <- 0
for (i in seq_along(args$i)) {
  go_df <- read.csv(file=args$i[i], header=TRUE, sep=",")

	go_df$GeneRatio <- frac2dec(go_df$GeneRatio)
	go_df$BgRatio <- frac2dec(go_df$BgRatio)
	go_df$Enrichment <- go_df$GeneRatio / go_df$BgRatio

	go_df <- go_df[go_df$p.adjust < 0.05,]
	
  go_df$Color <- gg_color_hue(i)#toString(i)
  
  max <- max(nrow(go_df), max)
# from Tommy's code
	p <- ggplot(go_df, aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust)))) + 
     	     geom_point(aes(size = Enrichment), color = gg_color_hue(i)) +
     	     scale_x_continuous(limits = c(0, ceiling(max(-log10(go_df$pvalue)))), breaks = pretty_breaks()) +
	         scale_size_area(max_size = 4, guide = guide_legend(title.position = "top")) + 
	         scale_y_discrete(labels = function(x) str_wrap(x, width = 15)) +
	         theme_bw() + 
	         theme(legend.title.align = 0.5, panel.border = element_blank(), panel.grid.major = element_blank(),
	               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
	         xlab(expression('-log'[10]*'(FDR)')) +
	         ylab(NULL) +
	         My_Theme 
  plots[[i]] <- p
}

plots <- ggarrange(plotlist = plots, ncol = length(plots), align = 'hv', common.legend = TRUE, legend = 'bottom') 
ggsave(args$o, plot=plots, height = max / 3, width = 3.4, units = 'in')
