#Returns Seurat Object with donor information and stimulation condition as metadata
loadData <- function(filename, stimname, donorname) {
  raw_counts <- read.table(file = filename, header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)
  raw <- CreateSeuratObject(counts = raw_counts, project = "mcsf-gmcsf", min.cells = 3, min.features = 200)
  raw$stim <- stimname
  raw$donor <- donorname
  Idents(raw) <- 'stim'
  return(raw)
}

#Takes in Seurat Object and returns updated Object with percent mitochondria content in slot "percent.mt"
calcMito <- function(raw) {
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = raw)), value = TRUE)
  nmito <- Matrix::colSums(GetAssayData(object = raw, slot = "counts")[mito.genes, ])
  ncounts <- Matrix::colSums(GetAssayData(object = raw, slot = "counts"))
  raw[['percent.mt']] <- nmito / ncounts * 100
  return(raw)
}

#Takes in Seurat Object and 
performQC <- function(raw, plot = FALSE){
  #Calculate percent mitochondrial reads
  raw <- calcMito(raw)
                         
  if (plot) {
    scatter <- FeatureScatter(object = raw, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by = 'percent.mt')
    print(scatter)
    
    #Count QC
    counts = raw$nCount_RNA
    cts <- qplot(counts, geom = 'histogram', bins = 100)
    cts_lo <- qplot(counts[counts < 5000], geom = 'histogram', bins = 100)
    cts_hi <- qplot(counts[counts > 10000], geom = 'histogram', bins = 100)
    cts_grid <- plot_grid(cts, cts_lo, cts_hi)
    print(cts_grid)
    
    #Gene QC
    genes = raw$nFeature_RNA
    gen <- qplot(genes, geom = 'histogram', bins = 100)
    gen_lo <- qplot(genes[genes < 1000], geom = 'histogram', bins = 100)
    gen_grid <- plot_grid(gen, gen_lo)
    print(gen_grid)
    
    vln <- VlnPlot(object = raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'stim', ncol = 3)
    print(vln)
  }
  #Filter
  raw <- subset(raw, subset = nCount_RNA > 2500 & nCount_RNA < 75000 & percent.mt < 20)
  raw <- NormalizeData(raw, verbose = FALSE)
  raw <- FindVariableFeatures(raw, selection.method = "vst", nfeatures = 4000)
  return(raw)
}


scaleData <- function(raw) {
  raw <- ScaleData(raw, vars.to.regress = c("nUMI", "percent.mt"), verbose = FALSE)
  raw <- RunPCA(raw, features = VariableFeatures(object = raw), npcs = 100, verbose = FALSE)
  p2 <- DimPlot(raw)
  return(raw)
}

#Performs differential expression analysis to find markers of stimulation condition,
#returns table of results for each gene
markers <- function(data){
  #DefaultAssay(data) <- 'RNA'
  Idents(data) <- 'stim'
  avg.exp <- log1p(AverageExpression(data, verbose = FALSE)$RNA)
  avg.exp$gene <- rownames(avg.exp)
  p1 <- ggplot(avg.exp, aes(M, G)) + geom_point() + ggtitle("Mph")
  #print(p1)
  roc.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='roc', min.pct = 0, logfc.threshold = 0, verbose = FALSE)
  wilcox.markers <- FindMarkers(data, ident.1 = "G", ident.2 = "M", test.use='wilcox', min.pct = 0, logfc.threshold = 0, verbose=FALSE)
  
  g.response = merge(roc.markers, wilcox.markers, by=0, all=TRUE)
  print(head(g.response, n = 15))
  
  return(g.response)
}


library(dplyr)
library(Matrix)
library(Seurat)
library(cowplot)
library(ggplot2)
library(reticulate)
library(KEGGprofile)
library(topGO)
library(clusterProfiler)
library(enrichR)
library(org.Hs.eg.db)
library(msigdbr)
set.seed(666)

#Set up file names and conditions
m_files = c('./data/reads.Sample1_MCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt', './data/mcsf_day6_1.txt', './data/mcsf_day6_2.txt')
m_donors = c('2', '1', '1')
g_files = c('./data/gmcsf_day6_1.txt', './data/gmcsf_day6_2.txt', './data/reads.Sample2_gMCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt')
g_donors = c('1', '1', '2')
files = c(m_files, g_files)
donornames = c(m_donors, g_donors)
stimnames = c(rep('M', length(m_donors)), rep('G', length(g_donors)))

# Load the datasets
alldata <- lapply(1:length(files), function(i) i)
mergeruns <- lapply(1:(length(unique(g_donors)) + length(unique(m_donors))), function(i) i)
ctr = 0
for (i in seq_along(files)){
  alldata[[i]] <- loadData(files[i], stimnames[i],  donornames[i])
  # merge runs with the same donor and stimulation 
  if (i == 1 || donornames[i - 1] != donornames[i] || stimnames[i - 1] != stimnames[i]) {
    if (ctr > 0) { 
      #cellCycle(mergeruns[[ctr]])
      mergeruns[[ctr]] <- performQC(mergeruns[[ctr]], plot=TRUE) 
      mergeruns[[ctr]] <- scaleData(mergeruns[[ctr]])
    }
    ctr = ctr + 1
    mergeruns[[ctr]] <- alldata[[i]]
  } else {
    mergeruns[[ctr]] <- merge(mergeruns[[ctr]], alldata[[i]])
    
  }
}
mergeruns[[length(mergeruns)]] <- performQC(mergeruns[[length(mergeruns)]], plot=TRUE) 
scaleData(mergeruns[[length(mergeruns)]])
mergeruns <- merge(mergeruns[[1]], y = mergeruns[2:4])
#DefaultAssay(mergeruns) <- 'RNA'
mergeruns.markers <- markers(mergeruns)
markers <- mergeruns.markers[mergeruns.markers$myAUC > 0.70, ]
markers <- c(markers$Row.names)
print(markers)

plot.stim <- VlnPlot(mergeruns, features = markers, group.by ='stim', pt.size = 0.01, combine = TRUE)
print(plot.stim)
plot.donor <- VlnPlot(mergeruns, features = markers, split.by='stim', group.by='donor', pt.size = 0, combine=TRUE)
print(plot.donor)
plot.clust <- VlnPlot(mergeruns, features = markers, split.by='stim', pt.size = 0, combine = TRUE)
print(plot.clust)

genes = tibble::deframe(mergeruns.markers[])
allGenes = tibble::deframe(mergeruns.markers['myAUC'])
names(allGenes) <- mergeruns.markers$Row.names
allGenes <- sort(allGenes, decreasing = TRUE)

geneSel = allGenes[markers]
sampleGOdata <- new("topGOdata", description = "Simple session", ontology = "BP", 
                    allGenes = allGenes, geneSel = (function(x) x > 0.70),
                    nodeSize = 10, annot = annFUN.org, mapping = 'org.Hs.eg.db', ID="symbol")
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultTable <- GenTable(sampleGOdata, 
                        classicFisher = resultFisher, 
                        classicKS = resultKS, 
                        topNodes = 50, 
                        ranksOf = 'classicFisher'
)


EG_IDs = mget(names(allGenes), revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs <- unlist(EG_IDs, use.names=F)
gmtfile <- system.file("extdata", "c5.mf.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

KEGGgenes <- allGenes
names(KEGGgenes) <- EG_IDs

markerKEGG <- mget(markers, revmap(org.Hs.egSYMBOL), ifnotfound=NA)
kk <- enrichKEGG(gene = unlist(markerKEGG),
                 organism = 'hsa',
                 universe = unlist(EG_IDs),
                 pvalueCutoff = 0.05)

kk <- find_enriched_pathway(markerKEGG, species = 'hsa')
kk <- gseKEGG(gene         = KEGGgenes,
              organism     = 'hsa',
              nPerm        = 1000,
              minGSSize    = 50,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
egmt <- enricher(EG_IDs, TERM2GENE=c5)

mkegg_df <- (msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") 
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mreact_df <- (msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
              %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mbio_df <- (msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
            %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mgomf_df <- (msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())
mgobp_df <- (msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
             %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame())

background <- Reduce(intersect, 
                     lapply(alldata, function(x) GetAssayData(object = x)@Dimnames[[1]]))
background <- Matrix::rowSums(GetAssayData(mergeruns, slot = 'counts')[background, ])
background <- names(background[background > 100])

msigk <- enricher(gene = markers, TERM2GENE = mkegg_df, universe = background)
msigmf <- enricher(gene = markers, TERM2GENE = mgomf_df, universe = background)
msigbp <- enricher(gene = markers, TERM2GENE = mgobp_df, universe = background)
mreact <- enricher(gene = markers, TERM2GENE = mreact_df, universe = background)
mbio <- enricher(gene = markers, TERM2GENE = mbio_df, universe = background)

gseak <- GSEA(allGenes, TERM2GENE = mkegg_df)
gseamf <- GSEA(allGenes, TERM2GENE = mgomf_df)
gseabp <- GSEA(allGenes, TERM2GENE = mgobp_df)
