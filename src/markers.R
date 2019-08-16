markers <- function(data, quick = FALSE){
  Idents(data) <- 'stim'
  if (quick) {
    fc.thresh = 0.25
  } else {
    fc.thresh = 0
  }
  roc.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='roc', min.pct = 0, logfc.threshold = fc.thresh, verbose = FALSE)
  wilcox.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='wilcox', min.pct = 0, logfc.threshold = fc.thresh, verbose=FALSE)
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
parser$add_argument('--quick', action='store_false', default=FALSE,
                    help='omit full data analysis')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
mergeruns.markers <- markers(mergeruns, quick=args$quick)
markers <- mergeruns.markers[mergeruns.markers$myAUC > 0.70, ]
markers <- c(markers$Row.names)

saveRDS(mergeruns.markers, paste('./data/', args$o, '.rds', sep = ''))
markerfile<-file(paste('./data/', args$o, '.txt', sep = ''))
writeLines(markers, markerfile)
close(markerfile)
