library(argparse)
library(Seurat)
library(cowplot)
library(ggplot2)

source('src/plotutils.R')
My_Theme <- get_default_theme()


parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='data file for input')
parser$add_argument('-m', type='character', nargs=1,
                    help='txt file containing list of markers to plot')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
parser$add_argument('--full', action='store_true', default=FALSE,
                    help='plot all cells on heatmap')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
con <- file(args$m, open="r")
markers <- readLines(con)
close(con)

mergeruns$stimdonor <- paste0(mergeruns$stim, mergeruns$donor)
Idents(mergeruns) <- 'stimdonor'
levels(mergeruns) <- c('G1', 'G2', 'M1', 'M2')

colorbar <- c(gg_color_hue(1), gg_color_hue(1), gg_color_hue(2), gg_color_hue(2))
if (args$full) {
  cells = mergeruns
} else {
  cells = subset(mergeruns, downsample = 250)
}
heatmap <- DoHeatmap(cells, features = markers, group.by = 'stimdonor', size = 6 * 5 / 14, 
                     group.bar = TRUE, group.colors = colorbar, label = TRUE, group.bar.height = 1 / length(markers),
                     draw.line = TRUE, lines.width = 15, angle = 0) + 
           guides(colour = 'none') +
           guides(fill = guide_colourbar(title.position = "top")) +
           My_Theme +
           theme(legend.position="bottom", legend.title.align = 0.5, legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(-18, 0, 0, 0), plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = 'in'),
                 axis.text.x=element_blank())


ggsave(args$o, plot=heatmap, height=length(markers) / 15 + 1, width = 2, units = 'in')
