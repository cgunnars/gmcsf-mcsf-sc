markers <- function(data){
  Idents(data) <- 'stim'
  roc.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='roc', min.pct = 0, logfc.threshold = 0, verbose = FALSE)
  wilcox.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='wilcox', min.pct = 0, logfc.threshold = 0, verbose=FALSE)
  g.response = merge(roc.markers, wilcox.markers, by=0, all=TRUE)
  return(g.response)
}

library(Seurat)
library(argparse)
set.seed(666)

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='cleaned file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='filestem for marker output')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
mergeruns.markers <- markers(mergeruns)
markers <- mergeruns.markers[mergeruns.markers$myAUC > 0.70, ]
markers <- c(markers$Row.names)

save(mergeruns.markers, file = paste(args$o, '.rds', sep = ''))
markerfile<-file(paste(args$o, '.txt', sep = ''))
writeLines(markers, markerfile)
close(markerfile)