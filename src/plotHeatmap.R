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
mergeruns$stimdonor <- paste0(mergeruns$stim, mergeruns$donor)
Idents(mergeruns) <- 'stimdonor'
heatmap <- DoHeatmap(mergeruns, features = markers)
ggsave(args$o, plot=heatmap)
