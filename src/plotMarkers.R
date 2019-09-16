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

plot.stim <- VlnPlot(mergeruns, features = markers, group.by = 'stim', pt.size = 0, combine = FALSE, return.plotlist = TRUE)
max = lapply(1:length(plot.stim), c)
for (i in seq_along(plot.stim)){
  max[[i]] <- max(FetchData(mergeruns, vars = markers[i]))
  p <- plot.stim[[i]]
  p$layers[[1]]$aes_params$size = 0.5 * 5 / 14
  plot.stim[[i]] <- p + 
                    geom_jitter(size = 0.2, stroke = 0) +
                    geom_boxplot(width = 0.1, fill = 'white', lwd = 0.5 * 5 / 14, outlier.shape = NA) +
                    My_Theme +
                    theme(plot.margin = margin(0, 0, 0.05, 0, unit = 'in'), plot.title = element_text(vjust = 1, margin = margin(0))) + 
                    xlab(NULL)
  
  if ((i %% args$ncol) == 1) {
    plot.stim[[i]] <- plot.stim[[i]] + ylab('ln(TPM+1)')
  } else {
    plot.stim[[i]] <- plot.stim[[i]] + ylab(NULL)
  }
  if (max[[i]] < 5) {
    plot.stim[[i]] <- plot.stim[[i]] + scale_y_continuous(breaks=seq(0, round(max[[i]]), by=1))
  }
}
plot.stim <- CombinePlots(plots = plot.stim, ncol = args$ncol, legend = 'none')
ggsave(paste(args$o, '-stim.svg', sep=''), plot=plot.stim, height = ceiling(length(markers) / args$ncol) * 0.85, width = args$ncol * 0.85)

plot.donor <- VlnPlot(mergeruns, features = markers, split.by='stim', group.by='donor', pt.size = 0, combine = FALSE, return.plotlist = TRUE)
for (i in seq_along(plot.donor)){
  p <- plot.donor[[i]]
  p$layers[[1]]$aes_params$size = 0.5 * 5 / 14
  plot.donor[[i]] <- p + 
                     My_Theme + 
                     theme(plot.margin = margin(0, 0, 0.05, 0, unit = 'in'), plot.title = element_text(vjust = 1, margin = margin(0))) + 
                     #geom_boxplot(width=0.1, fill = 'white', outlier.shape = NA) +
                     xlab(NULL) 
  if ((i %% args$ncol) == 1) {
    plot.donor[[i]] <- plot.donor[[i]] + ylab('ln(TPM+1)')
  } else {
    plot.donor[[i]] <- plot.donor[[i]] + ylab(NULL)
  }
  if (max[[i]] < 5) {
    plot.donor[[i]] <- plot.donor[[i]] + scale_y_continuous(breaks=seq(0, round(max[[i]]), by=1))
  }
}
plot.donor <- CombinePlots(plot = plot.donor, ncol = args$ncol, legend = 'right')
print(ceiling(length(markers) / args$ncol) * 0.85)
ggsave(paste(args$o, '-donor.svg', sep=''), plot=plot.donor, 
       height = ceiling(length(markers) / args$ncol) * 0.85 + 0.1, width = args$ncol * 0.85, units = 'in')


