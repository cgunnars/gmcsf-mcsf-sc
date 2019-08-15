#Returns Seurat Object with donor information and stimulation condition as metadata
loadData <- function(filename, stimname, donorname) {
  raw_counts <- read.table(file = filename, header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
  raw <- CreateSeuratObject(counts = raw_counts, project = "mcsf-gmcsf", min.cells = 3, min.features = 200)
  raw$stim <- stimname
  raw$donor <- donorname
  Idents(raw) <- 'stim'
  return(raw)
}

loadAll <- function(files, mergeruns, stimnames, donornames, outfile) {
  # Load the datasets
  alldata <- lapply(1:length(files), function(i) i)
  ctr = 0
  for (i in seq_along(files)){
    alldata[[i]] <- loadData(files[i], stimnames[i],  donornames[i])
    # merge runs with the same donor and stimulation 
    if (i == 1 || donornames[i - 1] != donornames[i] || stimnames[i - 1] != stimnames[i]) {
      ctr = ctr + 1
      mergeruns[[ctr]] <- alldata[[i]]
    } else {
      mergeruns[[ctr]] <- merge(mergeruns[[ctr]], alldata[[i]])
    }
    
    save(mergeruns, file=outfile)
    return(mergeruns)
  }
}

library(dplyr)
library(Matrix)
library(Seurat)
library(argparse)
set.seed(666)

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('--mfiles', type="character", nargs='+',
                     help='MCSF data files')
parser$add_argument('--gfiles', type='character', nargs='+',
                     help='GMSCF data files')
parser$add_argument('--mdonors', type='character', nargs='+',
                     help='MCSF donor names')
parser$add_argument('--gdonors', type='character', nargs='+',
                     help='GMCSF donor names')
parser$add_argument('-o', '--out', type='character', nargs=1,
                     default='./data/raw-mergeruns.rds',
                     help='output file for Seurat objects')
args <- parser$parse_args()
#Set up file names and conditions
m_files = c('./data/reads.Sample1_MCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt', './data/mcsf_day6_1.txt', './data/mcsf_day6_2.txt')
m_donors = c('2', '1', '1')
g_files = c('./data/gmcsf_day6_1.txt', './data/gmcsf_day6_2.txt', './data/reads.Sample2_gMCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt')
g_donors = args$gdonors#g_donors = c('1', '1', '2')
files = c(m_files, g_files)

outfile = args$out#outfile = './data/raw-mergeruns.rds'

donornames = c(m_donors, g_donors)
stimnames = c(rep('M', length(m_donors)), rep('G', length(g_donors)))
mergeruns <- lapply(1:(length(unique(g_donors)) + length(unique(m_donors))), function(i) i)

loadAll(files, mergeruns, stimnames, donornames, outfile)
