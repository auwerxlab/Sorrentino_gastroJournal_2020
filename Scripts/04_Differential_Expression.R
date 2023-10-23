library("refGenome")
library("reshape2")
library("FactoMineR")
library("ggplot2")
library("edgeR")
library("ggrepel")
library("UpSetR")
library("parallel")
library("grid")
library("clusterProfiler")
library("enrichplot")
library("org.Mm.eg.db")
library("cowplot")
library("DOSE")
library("openxlsx")
library("RUVSeq")
library("EDASeq")
library("RColorBrewer")

LOGFC = 1
PVALUE = 0.05

# Step 1: Load Metadata
# -----------------------
metadata <- read.csv(file = "Data/metadata/20190306_RNAseq sample_info.txt")
colnames_metadata <- c("SampleID","condition","Batch_ID")

# Step 2: Read Counts From Star
# ------------------------------
# Get a list of all read counts in the mapped directory
read_count <- list.files(path = "Data/mapped")
read_count <- read_count[grepl("ReadsPerGene.out.tab",read_count)]
list_count <- lapply(read_count, function(x){
  count_1_condition <- read.delim(file = paste0("Data/mapped/",x),header = FALSE)
  count_1_condition <- count_1_condition[,c(1,2)]
  colnames(count_1_condition) <- c("gene_id","count")
  count_1_condition$SampleID <- gsub("ReadsPerGene.out.tab","",x)
  return(count_1_condition)
})
count_star <- do.call(rbind, list_count)
count_star <- dcast(count_star, gene_id ~ SampleID, value.var = "count")
count_mat <- data.matrix(count_star[,-1])
rownames(count_mat) <- count_star$gene_id

# Step 3: read annotation from gtf file
# --------------------------------------
ens <- ensemblGenome()
read.gtf(ens, "Data/GRCm38/release95/Mus_musculus.GRCm38.95.gtf")
symb <- getGenePositions(ens)
ens_2_symbol <- symb[,c('gene_id','gene_name')]

# Step 4: Defining Function for Filtering And Normalization
# ----------------------------------------------------------
dgList_normalisation_filter <- function(count_mat, group_data){
  # DGEList is the EdgeR data structure
  dgList <- DGEList(
    counts = count_mat, 
    genes = rownames(count_mat),
    group = group_data
  )
  # Counts per million
  countsPerMillion <- cpm(dgList)
  # Which gene have more than 1 cpm
  countCheck <- countsPerMillion > 1
  # Keep genes with more than 1 sample with at least 1 cpm
  # least case it would be in 2 samples with excatly 1 cpm 
  # expression. this is important for better normalization
  # (refer to library composition problem)
  keep <- which(rowSums(countCheck) >= 2)
  dgList <- dgList[keep,]
  # TMM normalisation
  dgList <- calcNormFactors(dgList, method="TMM")
  return(dgList)
}

# Step 5: Defining Function for Multi dimentional scaling
# --------------------------------------------------------
MDS_custom <- function(dgList, metadata, analysis, colnames_meta = colnames_metadata){
  # be to examine the samples for outliers and for other relationships. 
  # The function plotMDS produces a plot in which distances between samples
  # correspond to leading biological coefficient of variation (BCV) between those samples
  MDS <- plotMDS(dgList)
  # take the scaling factors in each dimension and rename of default column names
  MDS$df <- data.frame(MDS$cmdscale.out)
  colnames(MDS$df) <- c("Dim_1","Dim_2")
  # add row names same as samples
  MDS$df$ID <- rownames(MDS$df)
  # match the each sample to it's metadata using sample name
  row_index <- match(as.character(MDS$df$ID),metadata$Sample_Name)
  # get the columns that are in the list colnames_meta
  col_index <- colnames(metadata) %in% colnames_meta
  # join the the metadata of these columns to the df
  MDS$df <- cbind(
    MDS$df,
    metadata[row_index,col_index]
  )
  # geom_text(aes(label=MDS$df$ID),hjust=0, vjust=0, size=1.8)+
  g1 <- ggplot(MDS$df, aes(x = Dim_1, y = Dim_2, col = condition))+
    geom_point()+
    geom_label_repel(aes(label = MDS$df$ID),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', size = 3)+
    theme_classic()+
    scale_color_brewer(palette = "Set1")+
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(face = "plain")
    )+
    ggsave(paste0("Plots/count_qc/MDS_",analysis, ".pdf"),width = 5,height = 4)+
    ggsave(paste0("Plots/count_qc/MDS_",analysis, ".png"),width = 5,height = 4)
  return(MDS)
}


# Step 6: Defining Function for Principle Component Analysis
# -----------------------------------------------------------
PCA_custom <- function(dgList, metadata, analysis, colnames_meta = colnames_metadata){
  # generate the PCA for the data
  pca <- PCA(t(dgList$counts))
  # get the coordinates and add them to pca object then rename rows
  pca$df <- data.frame(pca$ind$coord)
  pca$df$ID <- rownames(pca$df)
  
  pca$df <- cbind(
    pca$df,
    metadata[match(as.character(pca$df$ID),metadata$Sample_Name),
             colnames(metadata) %in% colnames_meta]
  )
  
  pca$df$ID = c("1", "2", "3", "1", "2", "3")
  pca$df$condition = gsub("INT", "INT-777", pca$df$condition)
  g <- ggplot(pca$df, aes(x = Dim.1, y = Dim.2))+
    geom_point(aes(colour=condition), size=5)+
    geom_point(shape = 1, size = 5, colour = "black")+
    geom_label_repel(aes(label = pca$df$ID),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50', 
                     size = 4, 
                     show.legend = FALSE) +
    theme_cowplot()+
    ggtitle("Raw Data (No Batch Correction) \nPrincipal Component Analysis") +
    xlab(paste0("PC1: ",round(pca$eig[1,2],digits = 2),"%"))+
    ylab(paste0("PC2: ",round(pca$eig[2,2],digits = 2),"%"))+
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(face = "plain"),
      legend.direction = "horizontal",
      legend.justification="center" ,
      legend.box.just = "bottom",
      legend.title = element_blank()
    )+
    geom_rug('outside' = TRUE) +
    scale_colour_manual(values=c("#32f207", "#fffff0"))+
    ggsave(paste0("Plots/count_qc/pca_tmm_",analysis, ".pdf"),width = 7.5,height = 4)+
    ggsave(paste0("Plots/count_qc/pca_tmm_",analysis, ".png"),width = 7.5,height = 4)
  
  pca$ggplot <- g
  return(pca)
}


# Step 7: Exact test DEGs
# -----------------------------------------------------------
exact_DEG <- function(dgList, analysis, case){
  design <- model.matrix(~case)
  
  dgList <- estimateDisp(dgList, design)
  edgeR_result <- exactTest(dgList)
  edgeR_result$table$adj_p <- p.adjust(edgeR_result$table$PValue)
  edgeR_result$table$gene <- rownames(edgeR_result$table)
  edgeR_result$table$symbol <- ens_2_symbol$gene_name[match(edgeR_result$table$gene,
                                                            ens_2_symbol$gene_id)]
  edgeR_result$table <- edgeR_result$table[order(edgeR_result$table$adj_p),]
  
  edgeR_result$table$Signif <- abs(edgeR_result$table$logFC)>LOGFC & edgeR_result$table$adj_p<PVALUE
  
  edgeR_result$ggplot <- ggplot(edgeR_result$table, aes(x = logFC, y = -log10(adj_p), col = Signif, Genes = symbol))+
    geom_point()+
    scale_color_manual(values = c("grey","red"))+
    ggsave(paste0("Plots/DEGs/Deg_",analysis, ".pdf"),width = 4,height = 4)+
    ggsave(paste0("Plots/DEGs/Deg_",analysis, ".png"),width = 4,height = 4)
  
  png(filename=paste0("Plots/DEGs/BCV_",analysis, ".png"))
  plotBCV(dgList)
  dev.off()
  
  fit <- glmFit(dgList, design)
  lrt <- glmLRT(fit)
  o <- order(lrt$table$PValue)
  png(filename=paste0("Plots/DEGs/LRT_",analysis, ".png"))
  plotMD(lrt)
  dev.off()
  
  return(edgeR_result)
}

# Step 8: Defining Function for Gene Enrichment Analysis
#-------------------------------------------------------
set.seed(123)
performGSEA <- function(result, analysis){
  gs <- setNames(result$table$logFC, nm = rownames(result$table))
  gs <- sort(gs, decreasing = T)
  
  # generate the gene set enrichment analysis
  gsego <- gseGO(geneList = gs,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 keyType = "ENSEMBL",
                 nPerm = 100000)
  
  # generate dot plot of the GSEA top 30
  png(filename=paste0("Plots/GSEA/Dot_Plot_",analysis, ".png"), height = 1500, width = 1500, res = 100)
  plot(dotplot(gsego, showCategory=30))
  dev.off()
  
  # generate network graph of GSEA result
  png(filename=paste0("Plots/GSEA/Graph_Plot_",analysis, ".png"), height = 3500, width = 3500, res = 400)
  plotGOgraph(gsego)
  dev.off()
  
  
  return(gsego)
}





# Do analysis on all samples
# -------------------------
# specify the groups 
group_data <- c(1,1,1,2,2,2)
# generate the dgList, filter and normalise
dgList <- dgList_normalisation_filter(count_mat, group_data)
# generate mds 
mds <- MDS_custom(dgList, metadata, "All_Data")
# generate pca
pca <- PCA_custom(dgList, metadata, "All_Data")
# perform exact test and generate plots
case <- factor(c(1,1,1,2,2,2))
de <- exact_DEG(dgList, "All_Data", case)
# perform gene set enrichment analysis
gsea <- performGSEA(de, "All_Data")

# compute cpm, remove useless genes & TMM normalize
dgList <- dgList_normalisation_filter(count_mat, metadata$condition)
# define design matrix
design <- model.matrix( ~ case)
# estimate common & tagwise dispursion
dgList <- estimateGLMCommonDisp(dgList, design)
dgList <- estimateGLMTagwiseDisp(dgList, design)
# fit linear model & train it
fit <- glmFit(dgList, design)
lrt_all_samples <- glmLRT(fit, coef=2)
lrt_all_samples$table$Signif <- abs(lrt_all_samples$table$logFC) > LOGFC & lrt_all_samples$table$PValue < PVALUE
lrt_all_samples$table$symbol <- ens_2_symbol$gene_name[match(rownames(lrt_all_samples$table),
                                                        ens_2_symbol$gene_id)]

topTags(lrt_all_samples)

# Do analysis with batch 3 removed
# --------------------------------
# specify the groups 
group_data2 <- c(1,1,2,2)
# generate the dgList, filter and normalise
dgList2 <- dgList_normalisation_filter(count_mat[,c(1,2,4,5)], group_data2)
# generate mds 
mds2 <- MDS_custom(dgList2, metadata[c(1,2,4,5),], "NoOutlier")
# generate pca
pca2 <- PCA_custom(dgList2, metadata[c(1,2,4,5),], "NoOutlier")
# perform exact test and generate plots
case2 <- factor(c(1,1,2,2))
de2 <- exact_DEG(dgList2, "NoOutlier", case2)
# perform gene set enrichment analysis
gsea2 <- performGSEA(de2, "NoOutlier")

# compute cpm, remove useless genes & TMM normalize
dgList <- dgList_normalisation_filter(count_mat[,c(1,2,4,5)], group_data2)
# define design matrix
design <- model.matrix( ~ case2)
# estimate common & tagwise dispursion
dgList <- calcNormFactors(dgList, method="upperquartile")
dgList <- estimateGLMCommonDisp(dgList, design)
dgList <- estimateGLMTagwiseDisp(dgList, design)
# fit linear model & train it
fit <- glmFit(dgList, design)
lrt_no_out <- glmLRT(fit, coef=2)
lrt_no_out$table$Signif <- abs(lrt_no_out$table$logFC) > LOGFC & lrt_no_out$table$PValue < PVALUE
lrt_no_out$table$symbol <- ens_2_symbol$gene_name[match(rownames(lrt_no_out$table),
                                                 ens_2_symbol$gene_id)]



topTags(lrt_no_out)


# Test: Use RUVseq to deal with batch effect - using empirical control genes
#---------------------------------------------------------------------------

# compute cpm, remove useless genes & TMM normalize
dgList <- dgList_normalisation_filter(count_mat, metadata$condition)
# define design matrix
design <- model.matrix( ~ 0 + metadata$condition)
# estimate common & tagwise dispursion
dgList <- estimateGLMCommonDisp(dgList, design)
dgList <- estimateGLMTagwiseDisp(dgList, design)
# fit linear model & train it
fit <- glmFit(dgList, design)
lrt <- glmLRT(fit, coef=2)

set <- newSeqExpressionSet(as.matrix(dgList$counts),
                           phenoData = data.frame(metadata$condition, row.names=colnames(dgList$counts)))

colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-1, 1), col=colors[metadata$condition])
plotPCA(set, col=colors[metadata$condition], cex=1.2)
set <- betweenLaneNormalization(set, which="upper")

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set, empirical, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-1, 1), col=colors[metadata$condition])
plotPCA(set2, col=colors[metadata$condition], cex=1.2)

pca <- PCA(t(assayData(set2)$normalizedCounts))
# get the coordinates and add them to pca object then rename rows
pca$df <- data.frame(pca$ind$coord)
pca$df$ID <- rownames(pca$df)

pca$df <- cbind(
  pca$df,
  metadata[match(as.character(pca$df$ID),metadata$Sample_Name),
           colnames(metadata) %in% colnames_metadata]
)

pca$df$ID = c("1", "2", "3", "1", "2", "3")
pca$df$condition = gsub("INT", "INT-777", pca$df$condition)
g <- ggplot(pca$df, aes(x = Dim.1, y = Dim.2))+
  geom_point(aes(colour=condition), size=5)+
  geom_point(shape = 1, size = 5, colour = "black")+
  stat_ellipse(aes(colour=condition, group=condition), type = "norm")+
  geom_label_repel(aes(label = pca$df$ID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   size = 4, 
                   show.legend = FALSE) +
  theme_cowplot()+
  ggtitle("RUV Data (Batch Correction) \nPrincipal Component Analysis") +
  xlab(paste0("PC1: ",round(pca$eig[1,2],digits = 2),"%"))+
  ylab(paste0("PC2: ",round(pca$eig[2,2],digits = 2),"%"))+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(face = "plain"),
    legend.direction = "horizontal",
    legend.justification="center" ,
    legend.box.just = "bottom",
    legend.title = element_blank()
  )+
  geom_rug('outside' = TRUE) +
  scale_colour_manual(values=c("#fffff0", "#32f207"))+
  ggsave(paste0("Plots/count_qc/pca_tmm_","RUV", ".pdf"),width = 7.5,height = 4)+
  ggsave(paste0("Plots/count_qc/pca_tmm_","RUV", ".png"),width = 7.5,height = 4)


# generate good plot grid for both PCAs for publication
# get the legend
lg <- get_legend(g)
# remove the legend from the other plots
pca_no_lg <- pca$ggplot + theme(legend.position = 'none');
g_no_lg <- g + theme(legend.position = 'none');

pdf("./Plots/count_qc/PCA_Final.pdf", width = 15, height = 5, useDingbats = FALSE, pointsize = 18)
pca_combined_plot <- plot_grid(pca_no_lg, g_no_lg, rel_widths = c(1,1), align = "h", nrow = 1)
pca_combined_plot <- plot_grid(pca_combined_plot, lg, rel_heights = c(1,0.1), nrow = 2, align = "h")
pca_combined_plot
dev.off()  

design <- model.matrix( ~ metadata.condition + W_1, data = pData(set2))
y <- DGEList(counts=counts(set2), group=metadata$condition)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

lrt$table$Signif <- abs(lrt$table$logFC) > LOGFC & lrt$table$PValue < PVALUE
lrt$table$symbol <- ens_2_symbol$gene_name[match(rownames(lrt$table),
                                                          ens_2_symbol$gene_id)]

topTags(lrt)






# ---------------------------------------------------
# Volcano Plots for before and after normalization
# ---------------------------------------------------
Pval_cutoff= 0.05
logFC_cutoff=1
labels = c("") # c("Sh2b2", "Nid1", "Cdh16", "Gjb4", "Lrp2", "Ntf5", "Pear1", "Syt8", "Clcf1")

# get data for this contrast
raw_df <- lrt_all_samples$table
RUV_df <- lrt$table
  
# add level tag to each data frame
raw_df$level <- "Raw Data \n(No Batch Correction)"
RUV_df$level <- "RUV Corrected"
  
# get same columns and merge data
merged_df <- rbind(raw_df, RUV_df)

# Set color for significant hits
merged_df$color <- ifelse(abs(merged_df$logFC) > logFC_cutoff & merged_df$PValue < Pval_cutoff, "#1453a6", "#4c4d4f")
merged_df$signf <- abs(merged_df$logFC) > logFC_cutoff & merged_df$PValue < Pval_cutoff
  
# keep only labels for genes in labels
merged_df$label <- merged_df$symbol
merged_df$label[!merged_df$label %in% labels] <- NA
  
# Number of significantly upregulated and downregulated genes (to be displayed as labels on the graph)
merged_df$up <- merged_df$signf & (merged_df$logFC > logFC_cutoff)
merged_df$down <- merged_df$signf & (merged_df$logFC < -logFC_cutoff)
  

data <- merged_df
data$level <- factor(data$level, levels =  c("Raw Data \n(No Batch Correction)", "RUV Corrected"), ordered = TRUE)

# get limits for all plots
Upper_limit <- ceiling(max(-log10(data$PValue))+1)
FC_limit <- ceiling(max(abs(data$logFC)))

# get the counts for each case
Up_all <- aggregate( up ~ level, data = data, FUN = sum)
Down_all <- aggregate( down ~ level, data = data, FUN = sum)

#data[is.na(data$label),]$label <- ""

volcano <-  ggplot(data = data, aes(x = logFC, y = -log10(PValue), color = color)) +
  # alpha based on significance
  geom_point(size = 2, data = data, aes(x = logFC, y = -log10(PValue), color = color)) +
  #geom_text_repel(aes(label = label), na.rm=TRUE, size = 4, segment.colour = "black", min.segment.length = 0.2,box.padding=0.25,segment.alpha = 0.5) +
  #geom_label_repel(aes(label = label), box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50', size = 4, show.legend = FALSE) +
  #facet_grid(~level) +
  labs("Volcano of gene expression") +
  xlab("log2 fold change") + 
  ylab("-log10 Adj. P-value") +
  scale_y_continuous(limits = c(0, Upper_limit)) +
  scale_x_continuous(limits = c(-FC_limit, FC_limit)) +
  scale_color_manual(values=c("#1453a6", "#4c4d4f")) +
  
  # Cutoff lines
  geom_vline(xintercept = logFC_cutoff, lwd = 0.4, lty = 2, col = "azure4") +
  geom_vline(xintercept = -logFC_cutoff, lwd = 0.4, lty = 2, col = "azure4") +
  geom_hline(yintercept = -log10(Pval_cutoff), lwd = 0.4, lty = 2, col = "azure4") +
  
  # Numbers of up/down as labels on the top right and top left
  geom_label(data = Up_all, mapping = aes(x = FC_limit, y = Upper_limit * 0.96, label = paste0(up," up")), hjust = 1, col = "black", size = 2.5 ) + 
  geom_label(data = Down_all, mapping = aes(x = -FC_limit, y = Upper_limit * 0.96, label = paste0(down," down")), hjust = 0, col = "black", size =2.5 ) +
  
  # Theme
  theme_bw() +
  theme(line = element_line(colour = 'black', size = 0.3), 
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "gray88" , colour="black", size = 0.3),  
        strip.text = element_text(size=10), strip.text.x = element_text(size = 6),
        legend.title = element_text(size = 17, face = "bold"),   legend.text = element_text(size = 12),
        panel.spacing = unit(0, "lines"), panel.border=element_rect(colour = 'black', size = 0.3))
  
  #stat_density2d(geom="tile", data = data[data$color=="#4c4d4f",], n = 300,inherit.aes = F, mapping = aes(x = logFC, y =-log10(adj.P.Val), fill=..density..^0.25, alpha=ifelse(..density..^0.25<0.25,0,1)), contour=FALSE) + 
  #scale_fill_gradientn(colours = colorRampPalette(c("white", "gray24"))(256))

ggsave(volcano, filename = paste0("./test.pdf"), height = 4.3, width = 6.5, useDingbats=FALSE )


# export DE results
wb <- createWorkbook()
addWorksheet(wb, sheetName = "All samples")
writeData(wb, sheet = "All samples", de$table)
addWorksheet(wb, sheetName = "No outlier")
writeData(wb, sheet = "No outlier", de2$table)
addWorksheet(wb, sheetName = "RUV_DE")
writeData(wb, sheet = "RUV_DE", lrt$table)
addWorksheet(wb, sheetName = "lrt_all_samples")
writeData(wb, sheet = "lrt_all_samples", lrt_all_samples$table)
addWorksheet(wb, sheetName = "lrt_no_out")
writeData(wb, sheet = "lrt_no_out", lrt_no_out$table)
saveWorkbook(wb, "./Data/DE.xlsx", overwrite = TRUE)

# export GSEA results
wb2 <- createWorkbook()
addWorksheet(wb2, sheetName = "All samples")
writeData(wb2, sheet = "All samples", gsea@result)
addWorksheet(wb2, sheetName = "No outlier")
writeData(wb2, sheet = "No outlier", gsea2@result)
saveWorkbook(wb2, "./Data/GSEA.xlsx", overwrite = TRUE)






































tmp = data.frame('NES' = gsea@result[-log(gsea@result[,'qvalues'])>4.5,'NES'], 
                 'minus_log_q_Value' = -log(gsea@result[-log(gsea@result[,'qvalues'])>4.5,'qvalues']), 
                 'label' = gsea@result[-log(gsea@result[,'qvalues'])>4.5,'Description'])

png(filename=paste0("Plots/GSEA/GSEA_","All_Data", ".png"), height = 3000, width = 3000, res = 300)
ggplot(tmp, aes(x=NES, y=minus_log_q_Value, label = label))+ 
  geom_point()+ 
  theme_cowplot()+ 
  geom_label_repel(aes(label = tmp[,'label']), box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50', size = 3)+
  theme_classic()
dev.off()



tmp = data.frame('NES' = gsea2@result[-log(gsea@result[,'qvalues'])>4.5,'NES'], 
                 'minus_log_q_Value' = -log(gsea2@result[-log(gsea@result[,'qvalues'])>4.5,'qvalues']), 
                 'label' = gsea2@result[-log(gsea@result[,'qvalues'])>4.5,'Description'])

png(filename=paste0("Plots/GSEA/GSEA_","NoOutlier", ".png"), height = 3000, width = 3000, res = 300)
ggplot(tmp, aes(x=NES, y=minus_log_q_Value, label = label))+ 
  geom_point()+ 
  theme_cowplot()+ 
  geom_label_repel(aes(label = tmp[,'label']), box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50', size = 3)+
  theme_classic()
dev.off()


tmp = data.frame('NES' = gsea@result[,'NES'], 
                 'minus_log_q_Value' = -log(gsea@result[,'qvalues']), 
                 'label' = gsea@result[,'Description'])
tmp[tmp[,'minus_log_q_Value'] < 4.7,][,'label'] = NA

png(filename=paste0("Plots/GSEA/GSEA_full_","All_Data", ".png"), height = 3000, width = 3000, res = 300)
ggplot(tmp, aes(x=NES, y=minus_log_q_Value, label = label))+ 
  geom_point()+ 
  geom_label_repel(aes(label = tmp[,'label']), box.padding = 0.35, point.padding = 0.5,segment.color = 'grey50', size = 1.5)+
  theme_classic()
dev.off()


tmp = data.frame('NES' = gsea2@result[,'NES'], 
                 'minus_log_q_Value' = -log(gsea2@result[,'qvalues']), 
                 'label' = gsea2@result[,'Description'])
tmp[tmp[,'minus_log_q_Value'] < 4.7,][,'label'] = NA

png(filename=paste0("Plots/GSEA/GSEA_full_","NoOutlier", ".png"), height = 3000, width = 3000, res = 300)
ggplot(tmp, aes(x=NES, y=minus_log_q_Value, label = label))+ 
  geom_point()+ 
  geom_label_repel(aes(label = tmp[,'label']), box.padding = 0.35, point.padding = 0.5,segment.color = 'grey50', size = 1.5)+
  theme_classic()
dev.off()



x = de$table[de$table[,'Signif'],'symbol']
y = de2$table[de2$table[,'Signif'],'symbol']
comp = intersect(x, y)
write(comp, file = "Reports/DE_Comp_Res.txt")

x2 = gsea[gsea$qvalues < 0.01,'Description']
y2 = gsea2[gsea2$qvalues < 0.01,'Description']
comp2 = intersect(x2, y2)
write(comp2, file = "Reports/GSEA_Comp_Res.txt")