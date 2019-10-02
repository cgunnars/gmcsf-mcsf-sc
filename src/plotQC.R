library(argparse)
library(Seurat)
library(cowplot)
library(ggplot2)
library(scales)
source('src/plotutils.R')
My_Theme <- get_default_theme()

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='data file for input')
parser$add_argument('-m', type='character', nargs=1,
                    help='txt file containing list of markers to plot')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
parser$add_argument('-ncol', type='integer', nargs=1,
                    default = 4, help='number of columns for plotting')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
con <- file(args$m, open="r")
markers <- readLines(con)
close(con)

mergeruns$stimdonor <- paste0(mergeruns$stim, mergeruns$donor)
Idents(mergeruns) <- 'stimdonor'
levels(mergeruns) <- c('G1', 'G2', 'M1', 'M2')
colors = c(gg_color_hue(1), gg_color_hue(1), gg_color_hue(2), gg_color_hue(2))
plot.stim <- VlnPlot(mergeruns, features = markers, cols = colors, group.by = 'stimdonor', pt.size = 0, combine = FALSE, return.plotlist = TRUE)
ylabels = c('Number of genes', 'Total counts', 'Mitochondrial reads (%)')
max = lapply(1:length(plot.stim), c)
left = 0
for (i in seq_along(plot.stim)){
  max[[i]] <- max(FetchData(mergeruns, vars = markers[i]))
  p <- plot.stim[[i]]
  p$layers[[1]]$aes_params$size = 0.5 * 5 / 14
  if (i == length(plot.stim)) {left = 0}
  plot.stim[[i]] <- p + 
                    geom_jitter(size = 0.2, stroke = 0) +
                    geom_boxplot(width = 0.1, fill = 'white', lwd = 0.5 * 5 / 14, outlier.shape = NA) +
                    My_Theme +
                    theme(plot.margin = margin(0.1, 0.1, 0.1, left, unit = 'in'), plot.title = element_text(vjust = 1, margin = margin(0))) + 
                    xlab(NULL) +
		    ylab(ylabels[i]) +
  
  if (max[[i]] < 5) {
    plot.stim[[i]] <- plot.stim[[i]] + scale_y_continuous(breaks=seq(0, round(max[[i]]), by=1))
  }
}
plot.stim <- CombinePlots(plots = plot.stim, ncol = args$ncol, legend = 'none')
ggsave(paste(args$o, '-stim.svg', sep=''), plot=plot.stim, height = ceiling(length(markers) / args$ncol) * 1.7, width = args$ncol * 1.7)

