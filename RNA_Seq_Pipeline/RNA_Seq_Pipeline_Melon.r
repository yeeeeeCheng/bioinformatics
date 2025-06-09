#Step 1. Prepare counts data and metadata of the samples and transform it into the format that can be processed by DESeq2.

library(tidyverse)

#Import the count matrix. Because the count matrix contains some basic information of that run in the first row that are not required for the following analysis, the attributes "skip = 1" can remove the first row whilte the attribute "header = FALSE" set the values in the second row as the column names.

counts_data <- read.csv('/Users/yeecheng/Bioinformatics/RNA_Seq_Pipeline/results/Quantification/USDA_Cmelo_AY_1.0/GCF_025177605.1_USDA_Cmelo_AY_1.0_count_matrix.csv', skip = 1, header = TRUE)

#Preview the first few rows in the dataframe counts_data

head(counts_data)

#Rename each column name. Only the codes of the samples were left, others are removed with the rename function. The function "remove_rownames" can delete the defalut row names (e.g. 1, 2, 3, 4, 5, 6, ... in this case). Once the default row names were removed, the values in the first column was set as the new row names while its column name would be removed.

counts_data_modified <- counts_data %>%
	rename(CK21 = "CK21_FR.bam",
		   CK22_1 = "CK22.1_FR.bam",
		   CK22_2 = "CK22.2_FR.bam",
		   CK23 = "CK23_FR.bam",
		   T21 = "T21_FR.bam",
		   T22 = "T22_FR.bam",
		   T23 = "T23_FR.bam",
		   T24 = "T24_FR.bam") %>%
	remove_rownames %>%
	column_to_rownames(var = "Geneid")
head(counts_data_modified)

#Confirm the column names of the dataframe counts_data_modified. The order of the sample names must be in accordance with the sample names provided in the metadata.

colnames(counts_data_modified)

#Import the sample metadata. The table can be created manually in the terminal.

col_data <- read.csv('/Users/yeecheng/Bioinformatics/RNA_Seq_Pipeline/results/Quantification/USDA_Cmelo_AY_1.0/Melon_sample_data.csv')
col_data

#Adjust the appearance of the metadata table as described above.

col_data_modified <- col_data %>%
	remove_rownames %>%
	column_to_rownames(var = "sample")
col_data_modified

#To confirm the level in the factor "treatment", the function "factor" can help us to declare each level in the factor.

col_data_modified$treatment <- factor(col_data_modified$treatment, levels = c("Treated", "Control"))
col_data_modified$treatment

#Make sure the row names in the col_data matches to the column names in counts_data_modified, if both the return values are TRUE, then we can proceed to Step 2.

all(rownames(col_data_modified) %in% colnames(counts_data_modified))
all(rownames(col_data_modified) == colnames(counts_data_modified))

#-----------------------------------------------------------------------#

#Step 2. Perform Differentially Expressed Gene Analysis with DESeq2.

library(DESeq2)
library(ggtext)
library(scales)
library(cowplot)

dds <- DESeqDataSetFromMatrix(countData = counts_data_modified,
							  colData = col_data_modified,
							  design = ~ treatment)
dds
#It's normal that the console might return a warning message saying that "some variables in design formula are characters, converting to factors". Since it's just a warning sign, it just tell us that this thing happened.

#It's important to tell DESeq2 which group serves as the control. Here, the function "relevel" was used to get this done. The line with the function "factor" tries to do the same thing to announce the levels in the factor "treatment". Yet, since we have already done this in Step. 1, maybe this line can be omitted. I don't remember any question arise when this line was omitted.

dds$treatment <- factor(dds$treatment, levels = c("Control", "Treated"))
dds$treatment <- relevel(dds$treatment, ref = "Control")

#Pre-filter the genes that have counts fewer than 10, which is a recommended cut-off threshold proposed by the official.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#In our case of study, the dimension of the table decrease almost 10000 of the rows

#Run DESeq2

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

#baseMean is the normalized result of each gene is each sample. log2FoldChange is the difference between the treated group and the untreated group expressed in the fold change. The positive value means the gene was upregulated in the treated group, while the negative value means the gene was downregulated. lfcSE means the standard estimate (or error?) of log2FoldChange. The column stat contain the result after performing Wald test. pvalue means the significance of the difference between the two groups, while padj means the corrected one. While using the result of padj to tell whether a gene was upregulated or downregulated after treatment, it can avoid the case of calling the false positive results.

#Set the threshold of padj (or the so-called "q value", the default correctionn method was FDR) to 0.05. the result table remain the same as the one generated with the command "res" as above, but the summary would show the number of the genes that met this criteria. since the summary table remain the same, maybe it's better to export the table and manipulate it in Excel to filter out genes which got NA in any of one column. here, we can export the result table of the object res into a csv file (optional), later on we can use excel or R to view the result table.
res05 <- results(dds, alpha = 0.05)
res05
summary(res05)
res05_df <- as.data.frame(res05)
#write.csv(res05_df, "/Users/yeecheng/desktop/res05_df.csv")

#Tag the genes according to the conditions of log2FC and padj. First, add a column of NAs, then give the data that met certain criteria the corresponding tags for the subsequent MA plot and volcano plot.

#Tags for the volcano plot. In the column created by the following functions, p-adj should less than 0.05 to be considered as the differentially expressed genes (DEG). If the log2FC is higher or lower than 1 or -1, the gene was considered as upregulated or downregulated DEG.
res05_df$diffexpressed <- "Non-DEGs"
#If log2Foldchange > 1 (FoldChange > 2) and padj < 0.05, set as "UP" 
res05_df$diffexpressed[res05_df$log2FoldChange > 1 & res05_df$padj < 0.05] <- "Upregulated"
# if log2Foldchange (FoldChange < 1/2) < -1 and padj < 0.05, set as "DOWN"
res05_df$diffexpressed[res05_df$log2FoldChange < -1 & res05_df$padj < 0.05] <- "Downregulated"

#Tags for the MA plot, yet the significance level of MA plot was set as 0.1. Whether are we going to use the value 0.1 as the significance level should be confirmed with PI.
res05_df$MASigs <- "MA-Non_Sigs"
res05_df$MASigs[res05_df$padj < 0.05] <- "MA_Sigs"

#Assign the colors for volcano plot to the corresponding data according to the tags
mycolors_Vol <- c("darkgoldenrod1", "black", "deepskyblue2")
names(mycolors_Vol) <- c("Upregulated", "Non-DEGs", "Downregulated")

# Re-plot but this time color the points with "diffexpressed"
res05_df %>%
	ggplot(aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+
	geom_point(alpha = 0.5)+
	geom_vline(xintercept = c(-1, 1), col="gray", linetype = "dashed")+
    geom_hline(yintercept = 1.301, col="gray", linetype = "dashed")+
    xlim(-7, 7)+
    ylim(0, 38)+
    xlab("log<sub>2</sub>FC")+
    ylab("-log<sub>10</sub>p-adj")+
    labs(color = "Types")+
    scale_color_manual(values = mycolors_Vol)+
	theme_classic()+
	theme(axis.title.x = element_markdown(),
		  axis.title.y = element_markdown(),
		  axis.text.x = element_text(hjust = 1, vjust = 1),
		  axis.text = element_text(color = "#000000", face = 2),
		  axis.line = element_line(linewidth = 0.8),
		  axis.ticks = element_line(linewidth = 0.8),
		  axis.title = element_text(face = 2),
		  legend.position = "inside",
		  legend.position.inside = c(0.15, 0.85),
		  legend.title = element_text(face = 2),
		  legend.text = element_text(face = 2))

#When it comes to outliers, just add the footnote to mention that there are still some outliers outside the range of the plot (but does it matter?). In this case, genes that are marked as outliers have NA in both the columns pvalue and padj. Genes marked as low counts have NA in the column padj.

#create the MA plot with the native DESeq2 function. x axis - counts; y axis MA plot - log2FoldChange
plotMA(res05, ylim = c(-10, 10))

#Generate the MA plot with ggplot2.
mycolors_MA <- c("black", "deepskyblue2")
names(mycolors_MA) <- c("MA-Non_Sigs", "MA_Sigs")

#MA plot powered by ggplot2
res05_df %>%
	ggplot(aes(x = baseMean, y = log2FoldChange, col = MASigs))+
	geom_point(alpha = 0.5)+
    geom_hline(yintercept = 0, col="darkgray", linetype = "dashed")+
    ylim(-7, 7)+
    xlab("Mean of normalizaed counts")+
    ylab("log<sub>2</sub>FC")+
    labs(color = "Types")+
    scale_color_manual(values = mycolors_MA)+
	scale_x_continuous(trans = "log10", n.breaks = 6)+
	theme_classic()+
	theme(axis.title.y = element_markdown(),
		  axis.text = element_text(color = "#000000", face = 2),
		  axis.line = element_line(linewidth = 0.8),
		  axis.ticks = element_line(linewidth = 0.8),
		  axis.title = element_text(face = 2),
		  legend.position = "right",
		  legend.title = element_text(face = 2),
		  legend.text = element_text(face = 2))

#create the PCA plot with the native DESeq2 function
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
plotPCA(vsd, intgroup=c("treatment"))

#PCA plot powered by ggplot2
pcaData <- plotPCA(vsd, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = treatment))+
	geom_point(size = 3)+
	geom_hline(yintercept = 0, col="black", linetype = "dashed", linewidth = 0.8)+
	geom_vline(xintercept = 0, col="black", linetype = "dashed", linewidth = 0.8)+
	xlab(paste0("PC1: ",percentVar[1],"% variance"))+
	ylab(paste0("PC2: ",percentVar[2],"% variance"))+
	labs(color = "Treatment")+
	theme(axis.text.x = element_text(hjust = 1, vjust = 1),
		  axis.text = element_text(color = "#000000", face = 2),
		  axis.ticks = element_line(linewidth = 0.8),
		  axis.title = element_text(face = 2),
		  legend.position = "bottom",
		  legend.title = element_text(face = 2),
		  legend.text = element_text(face = 2),
		  legend.key = element_rect(fill = "white", colour = "white"),
		  panel.background = element_rect(fill = "white", colour = "black"),
		  panel.border = element_rect(linetype = "solid", fill = NA, linewidth = 1))+
	coord_fixed()

#filter the DEGs, combine the condition of padj less than 0.05 and the absolute value of log2FoldChange larger than 1 with the function subset(), and extract the gene ids into a character vector.

resSig <- subset(res05, padj < 0.05 & abs(res05$log2FoldChange) > 1)
resSig <- na.omit(resSig)
genes_to_test <- rownames(resSig)
genes_to_test
resSig_df <- as.data.frame(resSig)

#-----------------------------------------------------------------------#

#Step 3. Perform GO enrichment analysis

library(clusterProfiler)
library(AnnotationDbi)
library(org.Cmelo.eg.db)

GO_mapped <- read.csv("/Users/yeecheng/desktop/GO.csv")

GO_mapped_tidy <- GO_mapped %>%
	gather(key = 'samples', value = 'GO_IDs', -From) %>%
	subset(GO_IDs != "") %>%
	dplyr::select(-samples)
GO_mapped_tidy

genes_to_test_2 <- append(genes_to_test, GO_mapped_tidy$GO_IDs)
genes_to_test_2_new <- genes_to_test_2[! genes_to_test_2 %in% genes_to_test]
genes_to_test_2_new

#GO terms in the category "Biological Process"
enrichGO_BP <- enrichGO(gene = genes_to_test_2_new,
					 OrgDb = org.Cmelo.eg.db,
					 ont = "BP",
					 keyType = "GO",
					 pvalueCutoff = 0.05,
					 pAdjustMethod = "BH",
					 qvalueCutoff = 0.05)
head(enrichGO_BP)
as.data.frame(enrichGO_BP)
write.csv(enrichGO_BP, "/Users/yeecheng/desktop/enrichGO_BP.csv")
enrichGO_BP_vis <- plot(barplot(enrichGO_BP, showCategory = 20))

png("enrichGO_BP.png", res = 250, width = 2400, height = 2400)
print(enrichGO_BP_vis)
dev.off()

enrichGO_BP_vis
goplot(enrichGO_BP)

#GO terms in the category "Cellular Component"
enrichGO_CC <- enrichGO(gene = genes_to_test_2_new,
					 OrgDb = org.Cmelo.eg.db,
					 ont = "CC",
					 keyType = "GO",
					 pvalueCutoff = 0.05,
					 pAdjustMethod = "BH",
					 qvalueCutoff = 0.05)
head(enrichGO_CC)
as.data.frame(enrichGO_CC)
write.csv(enrichGO_CC, "/Users/yeecheng/desktop/enrichGO_BP.csv")
enrichGO_CC_vis <- plot(barplot(enrichGO_CC, showCategory = 20))

png("enrichGO_CC.png", res = 250, width = 2400, height = 2400)
print(enrichGO_CC_vis)
dev.off()

enrichGO_CC_vis
goplot(enrichGO_CC)

#GO terms in the category "Molecular Function"
enrichGO_MF <- enrichGO(gene = genes_to_test_2_new,
					 OrgDb = org.Cmelo.eg.db,
					 ont = "MF",
					 keyType = "GO",
					 pvalueCutoff = 0.05,
					 pAdjustMethod = "BH",
					 qvalueCutoff = 0.05)
head(enrichGO_MF)
as.data.frame(enrichGO_MF)
write.csv(enrichGO_MF, "/Users/yeecheng/desktop/enrichGO_BP.csv")
enrichGO_MF_vis <- plot(barplot(enrichGO_MF, showCategory = 20))

png("enrichGO_MF.png", res = 250, width = 2400, height = 2400)
print(enrichGO_MF_vis)
dev.off()

enrichGO_MF_vis
goplot(enrichGO_MF)

#-----------------------------------------------------------------------#
#Step 4. Perform KEGG enrichment analysis

#Performing this step in the official KEGG website is more convenient than doing it in R, therefore the overview of the procedures will be presented in the slides.

#-----------------------------------------------------------------------#