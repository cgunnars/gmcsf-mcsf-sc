library(argparse)
library(Seurat)
library(cowplot)
library(ggplot2)
library(scales)

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

plot.stim <- VlnPlot(mergeruns, features = markers, group.by ='stim', pt.size = 0.01, combine = FALSE, return.plotlist = TRUE)
max = lapply(1:length(plot.stim), c)
for (i in seq_along(plot.stim)){
  max[[i]] <- max(FetchData(mergeruns, vars=markers[i]))
  plot.stim[[i]] <- plot.stim[[i]] + 
                    geom_boxplot(width=0.1, fill = 'white', outlier.shape = NA) +
                    xlab('Stimulation') +
                    ylab('ln(TPM+1)')
  if (max[[i]] < 5) {
    plot.stim[[i]] <- plot.stim[[i]] + scale_y_continuous(breaks=seq(0, round(max[[i]]), by=1))
  }
}
plot.stim <- CombinePlots(plots = plot.stim, ncol = args$ncol, legend = 'none')
ggsave(paste(args$o, '-stim.svg', sep=''), plot=plot.stim, height = length(markers) / args$ncol * 3, width = args$ncol * 3)

plot.donor <- VlnPlot(mergeruns, features = markers, split.by='stim', group.by='donor', pt.size = 0, combine = FALSE, return.plotlist = TRUE)
for (i in seq_along(plot.donor)){
  plot.donor[[i]] <- plot.donor[[i]] + 
                     #geom_boxplot(width=0.1, fill = 'white', outlier.shape = NA) +
                     xlab('Donor') +
                     ylab('ln(TPM+1)')
  if (max[[i]] < 5) {
    plot.donor[[i]] <- plot.donor[[i]] + scale_y_continuous(breaks=seq(0, round(max[[i]]), by=1))
  }
}
plot.donor <- CombinePlots(plot = plot.donor, ncol = args$ncol, legend = 'right')
ggsave(paste(args$o, '-donor.svg', sep=''), plot=plot.donor, height = length(markers) / args$ncol * 3, width = args$ncol * 3)