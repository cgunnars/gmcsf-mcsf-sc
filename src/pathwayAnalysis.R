library(clusterProfiler)
library(argparse)
library(org.Hs.eg.db)
library(dplyr)
library(msigdbr)
set.seed(666)

gs2eg <- function(genes) {
  EG_IDs <- mget(genes, revmap(org.Hs.egSYMBOL), ifnotfound=NA)
  EG_IDs <- unlist(EG_IDs, use.names=F)
  return(EG_IDs)
}

parser <- ArgumentParser(description='File i/o for single cell RNAseq')
parser$add_argument('-i', type="character", nargs=1,
                    help='scRNAseq data')
#parser$add_argument('-d', type='character', nargs=1,
#                    help = 'differential expression analysis results')
parser$add_argument('-m', type='character', nargs=1,
                    help='txt file containing list of significant markers')
parser$add_argument('-o', type='character', nargs=1,
                    help='file stem for output results')
args <- parser$parse_args()

con1 <- file(args$m, open="r")
markers <- readLines(con1)
close(con1)

con2 <- file(args$i, open='r')
background <- readLines(con2)
close(con2)

mkegg_df <- (msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") 
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mreact_df <- (msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
              %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mgomf_df <- (msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mgobp_df <- (msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())

print('KEGG started')
msigk <- setReadable(enrichKEGG(gene = gs2eg(markers), universe = gs2eg(background)), OrgDb = org.Hs.eg.db, keytype="ENTREZID")
write.csv(msigk, file=paste(args$o, '-kegg.csv', sep=''))
print('MF started')
msigmf <- setReadable(enrichGO(gs2eg(markers), org.Hs.eg.db, ont = 'MF', universe = gs2eg(background)), OrgDb = org.Hs.eg.db)
msigmf_simp <- simplify(msigmf, cutoff=0.6)
write.csv(msigmf@result, file=paste(args$o, '-mf.csv', sep=''))
write.csv(msigmf_simp@result, file=paste(args$o, '-mf-simp.csv', sep=''))
print('BP started')
msigbp <- setReadable(enrichGO(gs2eg(markers), org.Hs.eg.db, ont = 'BP', universe = gs2eg(background)), OrgDb = org.Hs.eg.db)
msigbp_simp <- simplify(msigbp, cutoff=0.6)
write.csv(msigbp@result, file=paste(args$o, '-bp.csv', sep=''))
write.csv(msigbp_simp@result, file=paste(args$o, '-bp-simp.csv', sep=''))
print('Reactome started')
mreact <- enricher(gene = markers, TERM2GENE = mreact_df, universe = background)
write.csv(mreact, file=paste(args$o, '-reactome.csv', sep=''))
