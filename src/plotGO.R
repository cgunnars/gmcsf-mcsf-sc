library(ggplot2)
library(cowplot)
library(argparse)
library(forcats)
library(scales)
library(ggpubr)

frac2dec <- function(frac) {
  return(sapply(frac, function(x) eval(parse(text=as.character(x)))))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[n]
}

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs='+',
                    help='data file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
args <- parser$parse_args()

plots <- lapply(1:length(args$i), c)

for (i in seq_along(args$i)) {
  go_df <- read.csv(file=args$i[i], header=TRUE, sep=",")

	go_df$GeneRatio <- frac2dec(go_df$GeneRatio)
	go_df$BgRatio <- frac2dec(go_df$BgRatio)
	go_df$Enrichment <- go_df$GeneRatio / go_df$BgRatio

	go_df <- go_df[go_df$p.adjust < 0.05,]
  go_df$Color <- i
# from Tommy's code
	p <- ggplot(go_df, aes(x = -log10(p.adjust), y = fct_reorder(Description, -log10(p.adjust)))) + 
     	     geom_point(aes(size = Enrichment), color = i) +
     	     theme_bw(base_size = 14) +
     	     scale_x_continuous(limits = c(0, NA), breaks = pretty_breaks()) +
     	     #scale_colour_manual(values = gg_color_hue(i)) +
    	     ylab(NULL) 
  plots[[i]] <- p
}

plots <- ggarrange(plotlist = plots, ncol = length(plots), align = 'hv', common.legend = TRUE, legend = 'right')
ggsave(args$o, plot=plots, width = 12)
