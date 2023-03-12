
sessionInfo()
rm(list=ls())

#################################################################################################

# Activating appropriate libraries
library(DESeq2)
library(tidyverse)
library(janitor)
library(xlsx)
library(data.table)
library(pheatmap)
library(ggplotify)
library(EnhancedVolcano)

# Set working directory
setwd("Experiment_1")
getwd()


#################################################################################################
#################################################################################################
#### Loading count data and metadata (sample information data)
#### Organizing and filtering to keep samples that are present all files
#################################################################################################
#################################################################################################


# Read in count data and restructure names to match sample info
countdata <- read.csv(file = 'input/220304_UNC51-VH00562_34_AAATNNVM5_2022-03-24_miRNA.csv')
countdata <- countdata %>% column_to_rownames("miRNA")
old_col_names <- colnames(countdata) #list of sample IDs with extended name
new_col_names <- str_extract(old_col_names, "\\d_\\d_\\d_\\d+") #use a regular expression to extract the sample ID portion from the extended name that matches the sample info IDs in the coldata file below
colnames(countdata) <- new_col_names #use the shortened sample IDs for the countdata such that the sample IDs match between datasets


# Check for duplicate miRNA IDs in the countdata 
summary(duplicated(rownames(countdata)))


# Read in metadata (sample information file)
coldata <- read.xlsx(file = "input/PFAS_miRNA_tracseqmanifest_RAGER.xlsx", sheetIndex = 1)


#################################################################################################
#################################################################################################
#### Create dataframes that are formatted for proceeding code, as well as DESeq2 functions
#################################################################################################
#################################################################################################


# Filter coldata to only include sample names from count data. This is already the case, so no changes.
samples <- colnames(countdata)
coldata <- subset(coldata, `Sample.ID..Required.` %in% all_of(samples))

# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$Sample.ID..Required.))

# Double checking that the same sample IDs appear between the two dataframes
setequal(as.character(coldata$Sample.ID..Required.), colnames(countdata))
identical(as.character(coldata$Sample.ID..Required.), colnames(countdata))


###################################################################################################
###################################################################################################
#### Background Filter
###################################################################################################
###################################################################################################

# Total number of samples = 12
nsamp <- ncol(countdata)

# Median expression level across all genes and all samples = 0
total_median <- median(as.matrix(countdata))

# Add back in the miRNA column
countdata <- countdata %>% rownames_to_column("miRNA")

# Get list of miRNAs that have an expression greater than the total median in at least 20% of the samples
miRNAs_above_background <- countdata %>% 
  pivot_longer(cols=!miRNA, names_to = "sampleID", values_to="expression") %>% 
  mutate(above_median=ifelse(expression>total_median,1,0)) %>% 
  group_by(miRNA) %>% 
  summarize(total_above_median=sum(above_median)) %>% 
  filter(total_above_median>.2*nsamp) %>% 
  select(miRNA)

# filter countdata for only the miRNAs above background
countdata <- left_join(miRNAs_above_background, countdata, by="miRNA")


# Write out filtered countdata
write_csv(countdata, paste0("output/all_samples/bg_20/RawCounts_AboveBack_bg20.csv"))

#################################################################################################
#################################################################################################
#### Subject Filter
#################################################################################################
#################################################################################################

# Transpose filtered countdata and confirm all samples have expression
countdata_T <- countdata %>% 
  pivot_longer(cols=!miRNA, names_to="sampleID",values_to="expression") %>% 
  pivot_wider(names_from=miRNA, values_from=expression) %>% 
  adorn_totals(where="col", name="rowsum")
 
# all samples have some expression so no need to remove anything
nrow(countdata_T %>% filter(rowsum==0))



#################################################################################################
#################################################################################################
#### RNASeq QA/QC on raw count data to identify potential outlier samples 
#################################################################################################
#################################################################################################

# make the miRNAs the rownames so we can perform PCA
countdata <- countdata %>% column_to_rownames("miRNA")

# Calculate principal components using transposed count data for both the unscaled and scaled data
pca <- prcomp(t(countdata))
screeplot(pca)


pca_scale <- prcomp(scale(t(countdata)))
screeplot(pca_scale)


# Make dataframe for PCA plot generation using first two components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample=colnames(countdata))
pca_df <- merge(pca_df, coldata %>% select(Sample.ID..Required., External.Code), by.x="Sample", by.y="Sample.ID..Required.")


pca_scale_df <- data.frame(PC1 = pca_scale$x[,1], PC2 = pca_scale$x[,2], Sample=colnames(countdata))
pca_scale_df <- merge(pca_scale_df, coldata %>% select(Sample.ID..Required., External.Code), by.x="Sample", by.y="Sample.ID..Required.")



# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
pca_scale_percent <- round(100*pca_scale$sdev^2/sum(pca_scale$sdev^2),1)


# Generating PCA plot annotated by Exposure Status
ggplot(pca_df, aes(PC1,PC2, color = External.Code))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
ggsave(paste0("figures/all_samples/bg_20/PCA_exposure_bg20.png"), width = 8, height = 5)
dev.off()


ggplot(pca_scale_df, aes(PC1,PC2, color = External.Code))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_scale_percent[1],"%)"), y=paste0("PC2 (",pca_scale_percent[2],"%)"))
ggsave(paste0("figures/all_samples/bg_20/PCA_exposure_scaled_bg20.png"), width = 8, height = 5)
dev.off()



# Generating PCA plot annotated by SampleID
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
ggsave(paste0("figures/all_samples/bg_20/PCA_IDs_bg20.png"), width = 10, height = 5)
dev.off()

ggplot(pca_scale_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_scale_percent[1],"%)"), y=paste0("PC2 (",pca_scale_percent[2],"%)"))
ggsave(paste0("figures/all_samples/bg_20/PCA_IDs_scaled_bg20.png"), width = 10, height = 5)
dev.off()

# run hierarchical clustering
p <- pheatmap(scale(t(countdata)), main="Hierarchical Clustering",
         cluster_rows=TRUE, cluster_cols = FALSE,
         fontsize_col = 7, treeheight_row = 60, show_colnames = FALSE)

p <- as.ggplot(p)
ggsave(paste0("figures/all_samples/bg_20/Hierarchical_Clustering_bg20.pdf"), width = 8, height = 5)

#Upon review, decided to keep all samples

#########################################################
# Set up DESeq2 experiment
#########################################################

#make sure exposure status, either PFAS or control, is a factor
coldata$External.Code <- as.factor(coldata$External.Code)

#design reflects Case vs. Control
dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~External.Code)

#estimate and check size factors
dds <- estimateSizeFactors(dds)
sizeFactors(dds) 

# normalized counts
normcounts <- counts(dds, normalized=TRUE)
write.csv(normcounts, paste0("output/all_samples/bg_20/NormCounts_bg20.csv"), row.names=TRUE)

# pseudocounts
ps_normcounts <- normcounts + 1
write.csv(ps_normcounts, paste0("output/all_samples/bg_20/NormCounts_ps_bg20.csv"),row.names=TRUE)

# log2 pseudocounts (y=log2(n+1))
log2normcounts <- log2(normcounts+1)
write.csv(log2normcounts, paste0("output/all_samples/bg_20/NormCounts_pslog2_bg20.csv"), row.names=TRUE)

# variance stabilizing matrix
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_matrix <-as.matrix(assay(vsd))
write.csv(vsd_matrix, paste0("output/all_samples/bg_20/VSDCounts_bg20.csv"), row.names=TRUE)


# run experiment
dds <- DESeq(dds, betaPrior=FALSE)

# get results for PFAS vs. Control and order by adjusted p-value
res <- results(dds, pAdjustMethod = "BH")
head(res)
ordered <- as.data.frame(res[order(res$padj),])
ordered <- ordered %>% rownames_to_column("miRNA")

# export results
write_csv(ordered, "output/all_samples/bg_20/PFAS_vs_Ctrl_statistical_results_bg20.csv")



# MA plot
# MA <- data.frame(res)

# pull non-significant genes, and significant up-, and down-regulated genes
# MA_ns <- MA %>% filter(padj>=0.01)
# MA_up <- MA %>% filter(padj<0.01 & log2FoldChange>1)
# MA_down <- MA %>% filter(padj<0.01 & log2FoldChange<(-1))

# Plot data with counts on x-axis and log2 fold change on y-axis
# ggplot(MA_ns, aes(x = baseMean, y = log2FoldChange)) + 
#   geom_point(color="gray75", size = .5) + # Plot non-significant genes
#   geom_point(data = MA_up, color="firebrick", size=1, show.legend = TRUE) + # Plot the up-regulated significant genes
#   geom_point(data = MA_down, color="dodgerblue2", size=1, show.legend = TRUE) + # Plot down-regulated significant genes
#   theme_bw() + # Change theme of plot from gray to black and white
#   scale_x_continuous(trans = "log10", breaks=c(1,10,100, 1000, 10000, 100000, 1000000), # We want to log10 transform x-axis for better visualizations
#                      labels=c("1","10","100", "1000", "10000", "100000", "1000000")) +   
#   xlab("Expression (Normalized Count)") + # Add labels for axes
#   ylab(" log2 Fold Change") + 
#   labs(title="MA Plot PFAS EV") +   # Add title
#   geom_hline(yintercept=0) # Add horizontal line at 0


# ggsave(paste0("figures/all_samples/bg_20/MA_bg20.png"), width = 8, height = 5)
