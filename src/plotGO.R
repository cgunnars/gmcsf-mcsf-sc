library(ggplot2)
library(cowplot)
library(argparse)
library(forcats)

frac2dec <- function(frac) {
  return(sapply(frac, function(x) eval(parse(text=as.character(x)))))
}

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='data file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
args <- parser$parse_args()
go_df <- read.csv(file=args$i, header=TRUE, sep=",")

go_df$GeneRatio <- frac2dec(go_df$GeneRatio)
go_df$BgRatio <- frac2dec(go_df$BgRatio)
go_df$Enrichment <- go_df$GeneRatio / go_df$BgRatio

go_df <- go_df[go_df$p.adjust < 0.05,]

# from Tommy's code
p <- ggplot(go_df, aes(x = Count, y = fct_reorder(Description, Count))) + 
     geom_point(aes(size = -log10(p.adjust), color = Enrichment)) +
     theme_bw(base_size = 14) +
     scale_colour_gradient() +
     ylab(NULL) +
     ggtitle("GO pathway enrichment")

ggsave(args$o, plot=p)