library(argparse)
library(dplyr)
source('src/plotutils.R')

args <- plot_io()
res <- readRDS(args$i)


dropnames <- c('avg_logFC', 'pct.1.y', 'pct.2.y')
res <- res[, -which(names(res) %in% dropnames)]
print(head(res))
res <- rename(res, gene = Row.names, AUC = myAUC, avg_log2FC = avg_diff, pct.g = pct.1.x, pct.m = pct.2.x)

print(head(res))

write.csv(res, file=args$o) 
