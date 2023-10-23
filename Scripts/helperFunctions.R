# helperFunctions
library(cowplot)
library(ggplot2)
geneIDtoSymbol <- function(symb = exp$symb, gid = "WBGene00197333"){
  symb$gene_name[match(gid, symb$gene_id)]
}
geneSymboltoID <- function(symb = exp$symb, gs = "cTel3X.2"){
  symb$gene_id[match(toupper(gs), toupper(symb$gene_name))]
}

plotGene <- function(exp, gs){
  e <- data.matrix(exp$voom$E) 
  dt <- exp$samples
  dt$expressionLevel <- exp$voom$E[which(rownames(e) == geneSymboltoID(exp$symb, gs = gs)),]
  ggplot(dt, aes(y = expressionLevel, x = Genotype, color = Treatment)) + geom_point() + theme_cowplot() + ggtitle(gs)
}