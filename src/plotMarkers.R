library(argparse)
library(Seurat)
library(cowplot)
library(ggplot2)

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='data file for input')
parser$add_argument('-m', type='character', nargs=1,
                    help='txt file containing list of markers to plot')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for figure output')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
con <- file(args$m, open="r")
markers <- readLines(con)
close(con)

plot.stim <- VlnPlot(mergeruns, features = markers, group.by ='stim', pt.size = 0.01, combine = TRUE, ncol = 4)
ggsave(paste(args$o, '-stim.svg', sep=''), plot=plot.stim, height = length(markers) / 4 * 3, width = 12)

plot.donor <- VlnPlot(mergeruns, features = markers, split.by='stim', group.by='donor', pt.size = 0, combine=TRUE, ncol = 4)
ggsave(paste(args$o, '-donor.svg', sep=''), plot=plot.donor, height = length(markers) / 4 * 3, width = 12)