library("RcisTarget")
library("openxlsx")
library("readxl")
library("plyr")

featherURL <- "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-5kb-7species.mc9nr.feather" 
download.file(featherURL, destfile=basename(featherURL))

# Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_mgi)
motifRankings <- importRankings("./mm9-tss-centered-5kb-7species.mc9nr.feather")

# Without outliers 
# ----------------
no_out <- read_excel("./Data/DE.xlsx", sheet = "No outlier")
geneSet_no_out <- no_out[no_out$Signif & no_out$logFC >= 1,]$symbol
geneLists_no_out <- list(geneSetName = geneSet_no_out)

# Motif enrichment analysis:
motifEnrichmentTable_no_out <- cisTarget(geneLists_no_out, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations_mgi)


# All Samples
# -----------
out <- read_excel("./Data/DE.xlsx", sheet = "All samples")
geneSet_all_samples <- out[out$Signif & out$logFC >= 1,]$symbol
geneLists_all_samples <- list(geneSetName = geneSet_all_samples)

# Motif enrichment analysis:
motifEnrichmentTable_all_samples <- cisTarget(geneLists_all_samples, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations_mgi)


# RUV normalized
# --------------
RUV_out <- read_excel("./Data/DE.xlsx", sheet = "RUV_DE")
geneSet_all_samples <- RUV_out[RUV_out$Signif & RUV_out$logFC >= 1,]$symbol
geneLists_all_samples <- list(geneSetName = geneSet_all_samples)

# Motif enrichment analysis:
motifEnrichmentTable_RUV <- cisTarget(geneLists_all_samples, 
                                              motifRankings,
                                              motifAnnot=motifAnnotations_mgi)


# LRT all samples
# ---------------
RUV_out <- read_excel("./Data/DE.xlsx", sheet = "lrt_all_samples")
geneSet_all_samples <- RUV_out[RUV_out$Signif & RUV_out$logFC >= 1,]$symbol
geneLists_all_samples <- list(geneSetName = geneSet_all_samples)

# Motif enrichment analysis:
motifEnrichmentTable_all_samples <- cisTarget(geneLists_all_samples, 
                                      motifRankings,
                                      motifAnnot=motifAnnotations_mgi)

# LRT no outliers
# ---------------
RUV_out <- read_excel("./Data/DE.xlsx", sheet = "lrt_no_out")
geneSet_all_samples <- RUV_out[RUV_out$Signif & RUV_out$logFC >= 1,]$symbol
geneLists_all_samples <- list(geneSetName = geneSet_all_samples)

# Motif enrichment analysis:
motifEnrichmentTable_no_out <- cisTarget(geneLists_all_samples, 
                                      motifRankings,
                                      motifAnnot=motifAnnotations_mgi)


# export the results
wb <- createWorkbook()
addWorksheet(wb, sheetName = "no_outlier")
writeData(wb, sheet = "no_outlier", motifEnrichmentTable_no_out)
addWorksheet(wb, sheetName = "all_samples")
writeData(wb, sheet = "all_samples", motifEnrichmentTable_all_samples)
addWorksheet(wb, sheetName = "RUV")
writeData(wb, sheet = "RUV", motifEnrichmentTable_RUV)
addWorksheet(wb, sheetName = "lrt_all_samples")
writeData(wb, sheet = "lrt_all_samples", motifEnrichmentTable_all_samples)
addWorksheet(wb, sheetName = "lrt_no_out")
writeData(wb, sheet = "lrt_no_out", motifEnrichmentTable_no_out)
saveWorkbook(wb, paste0("./Reports/Morif_Enrichment.xlsx"), overwrite = TRUE)



# consider only ones with high confidence
data <- motifEnrichmentTable_RUV[!motifEnrichmentTable_RUV$TF_highConf == "",]

# labels minimal
data$label <- gsub("(.*)\\(.*\\).*","\\1", data$TF_highConf)
data$label <- gsub("1700080O16Rik; 3830417A13Rik; Magea1; Magea10; Magea2; Magea3; Magea4; Magea5; Magea6; Magea8", 
                   "Magea1; Magea10;", data$label)
data$label[24] <- "Tead4 (2nd Occurence)"
data$label[25] <- "Tead3 (2nd Occurence)"
data$label[10] <- "Egr1; Egr2; Egr3; Egr4; Wt (2nd Occurence)"
data$label <- factor(data$label, levels = unique(data$label))

# Sort by dose and supp
data_sorted <- data[order(data[,3],decreasing=TRUE),]

top_25_barplot <- ggplot(data = data[1:25,], aes(x=label, y= NES, fill = NES)) +
  geom_bar(stat="identity", width=0.9)+
  scale_fill_gradient(low = "green", high = "darkgreen")+
  geom_text(aes(label=NES), hjust=-0.3, size=3.5) +
  ggtitle("Top 25 Enriched Motifs") +
  xlab("Transcription factors annotated to the motif (High Confidence)") + 
  # Transcription factors annotated to the motif according to ‘motifAnnot_highConfCat’
  ylab("Normalized enrichment score") + # Normalized enrichment score of the motif in the gene-set
  theme(line = element_line(colour = 'black', size = 0.3), 
        panel.grid = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "gray88" , colour="black", size = 0.3),  
        strip.text = element_text(size=10), strip.text.x = element_text(size = 6),
        legend.title = element_text(size = 17, face = "bold"),   legend.text = element_text(size = 12),
        panel.spacing = unit(0, "lines"), panel.border=element_rect(colour = 'black', size = 0.3) )+ 
  coord_flip()

ggsave(top_25_barplot, filename = paste0("./Plots/top_25_Enriched_motifs.pdf"), height = 9, width = 10, useDingbats=FALSE )
