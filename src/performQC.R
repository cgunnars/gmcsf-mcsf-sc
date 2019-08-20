#Takes in Seurat Object and returns updated Object with percent mitochondria content in slot "percent.mt"
calcMito <- function(raw) {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = raw)), value = TRUE)
  nmito <- Matrix::colSums(GetAssayData(object = raw, slot = "counts")[mito.genes, ])
  ncounts <- Matrix::colSums(GetAssayData(object = raw, slot = "counts"))
  raw[['percent.mt']] <- nmito / ncounts * 100
  return(raw)
}

performQC <- function(raw, minct = 2500, maxct = 75000, mt = 20){
  #Calculate percent mitochondrial reads
  raw <- calcMito(raw)
  #Filter based on counts and mitochondrial reads
  expr <- FetchData(object = raw, vars = c('nCount_RNA'))
  raw <- raw[, which(x = (expr$nCount_RNA > minct & expr$nCount_RNA < maxct))]
  expr <- FetchData(object = raw, vars = c('percent.mt'))
  raw <- raw[, which(x = expr$percent.mt < mt)]
  return(raw)
}

backgroundSet <- function(mergeruns) {
  background <- Reduce(intersect, 
                       lapply(mergeruns, function(x) GetAssayData(object = x)@Dimnames[[1]]))
  background <- Matrix::rowSums(GetAssayData(merge(mergeruns[[1]], mergeruns[2:length(mergeruns)]), slot = 'counts')[background, ])
  
  background <- names(background[background > 100])
  
  outfile<-file('data/background-genes.txt')
  writeLines(background, outfile)
  close(outfile)
}



library(Matrix)
library(Seurat)
library(argparse)
set.seed(666)

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='raw data file for input')
parser$add_argument('--minct', type='integer', nargs=1,
                    default = 2500,
                    help='maximum counts/cell')
parser$add_argument('--maxct', type='integer', nargs=1,
                    default = 75000,
                    help='maximum counts/cell')
parser$add_argument('--mt', type='integer', nargs=1,
                    default = 20, help='mitochondrial content')
parser$add_argument('--regnUMI', action="store_true", default=TRUE,
                    help='Regress out nUMI [default]')
parser$add_argument('--regmito', action='store_true', default=TRUE,
                    help='Regress out mitochondrial content [default]')
parser$add_argument('-o', type='character', nargs=1,
                    help='qc data file')
args <- parser$parse_args()

mergeruns <- readRDS(args$i)

regress_out = vector()
if (args$regnUMI) {
  regress_out = append(regress_out, 'nUMI')
}
if (args$regmito) {
  regress_out = append(regress_out, 'percent.mt') 
}

for (i in seq_along(mergeruns)) {
  mergeruns[[i]] <- performQC(mergeruns[[i]], minct = args$minct, maxct = args$maxct, mt = args$mt)
}
backgroundSet(mergeruns)

mergeruns <- merge(mergeruns[[1]], y = mergeruns[2:length(mergeruns)])

#Normalize data
mergeruns <- NormalizeData(mergeruns, verbose = TRUE)
mergeruns <- ScaleData(mergeruns, vars.to.regress = regress_out, verbose = TRUE)
#mergeruns <- FindVariableFeatures(mergeruns, selection.method = "vst", nfeatures = 2000)

saveRDS(mergeruns, args$o)