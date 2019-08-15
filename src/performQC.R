#Takes in Seurat Object and returns updated Object with percent mitochondria content in slot "percent.mt"
calcMito <- function(raw) {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = raw)), value = TRUE)
  nmito <- Matrix::colSums(GetAssayData(object = raw, slot = "counts")[mito.genes, ])
  ncounts <- Matrix::colSums(GetAssayData(object = raw, slot = "counts"))
  raw[['percent.mt']] <- nmito / ncounts * 100
  return(raw)
}

performQC <- function(raw, min_count = 2500, max_count = 75000, mt = 20){
  #Calculate percent mitochondrial reads
  raw <- calcMito(raw)
  #Filter
  raw <- subset(raw, subset = nCount_RNA > min_count & nCount_RNA < max_count & percent.mt < mt)
  raw <- NormalizeData(raw, verbose = FALSE)
  raw <- FindVariableFeatures(raw, selection.method = "vst", nfeatures = 4000)
  return(raw)

  }


library(Matrix)
library(Seurat)
library(argparse)
set.seed(666)

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='raw data file for input')
parser$add_argument('-o', type='character', nargs=1,
                    help='qc data file')
args <- parser$parse_args()
mergeruns <- readRDS(args$i)
for (i in seq_along(mergeruns)) {
  run <- performQC(mergeruns[[i]])
  mergeruns[[i]] <- ScaleData(run, vars.to.regress = c("nUMI", "percent.mt"), verbose = FALSE)
}
mergeruns <- merge(mergeruns[[1]], y = mergeruns[2:4])
save(mergeruns, file=args$o)