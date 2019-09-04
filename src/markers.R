markers <- function(data, quick = FALSE){
  Idents(data) <- 'stim'
  if (quick) {
    fc.thresh = 0.25
  } else {
    fc.thresh = 0
  }
  roc.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='roc', min.pct = 0, logfc.threshold = fc.thresh, verbose = FALSE)
  wilcox.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='wilcox', min.pct = 0, logfc.threshold = fc.thresh, verbose=FALSE)
  print(head(roc.markers))
  print(head(wilcox.markers))
  g.response <- merge(x = roc.markers, y = wilcox.markers, by = 'row.names')
  print(head(g.response))
  return(g.response)
}

writeMarkers <- function(data, outstem, cond) {
  print(head(data))
  if (cond == 'g') {
    markers_all <- data[(data$p_val_adj < 0.05 & data$avg_diff > 0.25), ]
    markers <- data[(data$myAUC > 0.70), ]
  } else {
    markers_all <- data[(data$p_val_adj < 0.05 & data$avg_diff < -0.25), ]
    markers <- data[(data$myAUC < 0.30), ]
  }
  print(head(markers))
  print(head(markers_all))
  markers_all <- c(markers_all$Row.names)
  markers <- c(markers$Row.names)
  
  markers_all_file <- file(paste(outstem, '-all-' , cond, '.txt', sep = ''))
  writeLines(markers_all, markers_all_file)
  close(markers_all_file)

  markers_file <- file(paste(outstem, '-', cond, '.txt', sep=''))
  writeLines(markers, markers_file)
  close(markers_file)
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
parser$add_argument('--split', action='store_false', default=FALSE,
                    help='split by donor for analysis')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)

if(args$split) {
  mergeruns <- SplitObject(mergeruns, split.by='donor')
  markers.donor <- lapply(1:length(mergeruns), c)
  print(markers.donor)
  for (i in seq_along(mergeruns)) {
    markers.donor[[i]] <- markers(mergeruns[[i]], quick=args$quick)
    print(head(markers.donor[[i]]))
    writeMarkers(markers.donor[[i]], paste0(args$o, '-', i), 'g')
    writeMarkers(markers.donor[[i]], paste0(args$o, '-', i), 'm')
  }
  mergeruns.markers <- Reduce(function(x, y) merge(x = x, y = y, by = 'Row.names'),
                              markers.donor)
} else {
  mergeruns.markers <- markers(mergeruns, quick=args$quick)
  writeMarkers(mergeruns.markers, args$o, 'g')
  writeMarkers(mergeruns.markers, args$o, 'm')
}

saveRDS(mergeruns.markers, paste(args$o, '.rds', sep = ''))


