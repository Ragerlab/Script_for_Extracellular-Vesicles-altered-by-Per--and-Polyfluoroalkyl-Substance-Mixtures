miRNA Differential Expression from EVs Isolated from HepG2 Cells Treated
with a PFAS Mixture
================
Rager Lab
2023-03-12

``` r
# Activating appropriate libraries
library(DESeq2)
library(tidyverse)
library(janitor)
library(xlsx)
library(data.table)
library(pheatmap)
library(ggplotify)
library(EnhancedVolcano)
```

``` r
setwd("Experiment_1")
```

``` r
# Read in count data and restructure names to match sample info
countdata <- read.csv(file = 'input/220304_UNC51-VH00562_34_AAATNNVM5_2022-03-24_miRNA.csv')
countdata <- countdata %>% column_to_rownames("miRNA")

#list of sample IDs with extended name
old_col_names <- colnames(countdata) 

#use a regular expression to extract the sample ID portion from the extended name that matches the sample info IDs in the coldata file below
new_col_names <- str_extract(old_col_names, "\\d_\\d_\\d_\\d+") 

#use the shortened sample IDs for the countdata such that the sample IDs match between datasets
colnames(countdata) <- new_col_names 


# Check for duplicate miRNA IDs in the countdata 
summary(duplicated(rownames(countdata)))
```

    ##    Mode   FALSE 
    ## logical    1124

``` r
# Read in metadata (sample information file)
coldata <- read.xlsx(file = "input/PFAS_miRNA_tracseqmanifest_RAGER.xlsx", sheetIndex = 1)
```

``` r
## Create dataframes that are formatted for proceeding code, as well as DESeq2 functions

# Filter coldata to only include sample names from count data. This is already the case, so no changes.
samples <- colnames(countdata)
coldata <- subset(coldata, `Sample.ID..Required.` %in% all_of(samples))

# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$Sample.ID..Required.))

# Double checking that the same sample IDs appear between the two dataframes
setequal(as.character(coldata$Sample.ID..Required.), colnames(countdata))
```

    ## [1] TRUE

``` r
identical(as.character(coldata$Sample.ID..Required.), colnames(countdata))
```

    ## [1] TRUE

``` r
# Background Filter

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
```

``` r
# Subject Filter

# Transpose filtered countdata and confirm all samples have expression
countdata_T <- countdata %>% 
  pivot_longer(cols=!miRNA, names_to="sampleID",values_to="expression") %>% 
  pivot_wider(names_from=miRNA, values_from=expression) %>% 
  adorn_totals(where="col", name="rowsum")
 
# all samples have some expression so no need to remove anything
nrow(countdata_T %>% filter(rowsum==0))
```

    ## [1] 0

``` r
#RNASeq QA/QC on raw count data to identify potential outlier samples 

# make the miRNAs the rownames so we can perform PCA
countdata <- countdata %>% column_to_rownames("miRNA")

# Calculate principal components using transposed count data for both the unscaled and scaled data
pca <- prcomp(t(countdata))
screeplot(pca)
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-1.png)<!-- -->

``` r
pca_scale <- prcomp(scale(t(countdata)))
screeplot(pca_scale)
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-2.png)<!-- -->

``` r
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
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-3.png)<!-- -->

``` r
ggplot(pca_scale_df, aes(PC1,PC2, color = External.Code))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_scale_percent[1],"%)"), y=paste0("PC2 (",pca_scale_percent[2],"%)"))
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-4.png)<!-- -->

``` r
# Generating PCA plot annotated by SampleID
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-5.png)<!-- -->

``` r
ggplot(pca_scale_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_scale_percent[1],"%)"), y=paste0("PC2 (",pca_scale_percent[2],"%)"))
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-6.png)<!-- -->

``` r
# run hierarchical clustering
p <- pheatmap(scale(t(countdata)), main="Hierarchical Clustering",
         cluster_rows=TRUE, cluster_cols = FALSE,
         fontsize_col = 7, treeheight_row = 60, show_colnames = FALSE)
```

![](DESeq2_miRNA_bg20_files/figure-gfm/outlier%20ID-7.png)<!-- -->

No clear outliers so we will retain all samples.

``` r
# Set up DESeq2 experiment

#make sure exposure status, either PFAS or control, is a factor
coldata$External.Code <- as.factor(coldata$External.Code)

#design reflects Case vs. Control
dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~External.Code)
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

``` r
#estimate and check size factors
dds <- estimateSizeFactors(dds)
```

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

``` r
sizeFactors(dds) 
```

    ##   9_5_3_1   9_5_3_2   9_5_3_3   9_5_3_4   9_5_3_5   9_5_3_6   9_5_3_7   9_5_3_8 
    ## 0.5429027 0.3621633 0.7578586 1.6581753 3.4833112 3.0622006 0.5106434 0.5170499 
    ##   9_5_3_9  9_5_3_10  9_5_3_11  9_5_3_12 
    ## 0.3576277 2.4219811 2.5115672 1.8570315

``` r
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
```

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## final dispersion estimates

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## fitting model and testing

``` r
# get results for PFAS vs. Control and order by adjusted p-value
res <- results(dds, pAdjustMethod = "BH")
head(res)
```

    ## log2 fold change (MLE): External.Code PFAS vs Control. 
    ## Wald test p-value: External.Code PFAS vs Control. 
    ## DataFrame with 6 rows and 6 columns
    ##                         baseMean log2FoldChange     lfcSE      stat      pvalue
    ##                        <numeric>      <numeric> <numeric> <numeric>   <numeric>
    ## hsa-let-7a-1_mature    153.26735       -2.68671  0.870968 -3.084740 2.03730e-03
    ## hsa-let-7a-1_stemloop    2.48717        3.73370  2.836713  1.316207 1.88105e-01
    ## hsa-let-7a-2_mature   1876.59358        2.10311  0.379766  5.537913 3.06097e-08
    ## hsa-let-7a-3_mature     91.81008       -1.80577  0.988804 -1.826212          NA
    ## hsa-let-7b_mature      478.46639        0.32961  0.789660  0.417407 6.76381e-01
    ## hsa-let-7c_mature       58.64673       -1.61249  0.886937 -1.818046 6.90572e-02
    ##                              padj
    ##                         <numeric>
    ## hsa-let-7a-1_mature   1.00910e-02
    ## hsa-let-7a-1_stemloop 3.54936e-01
    ## hsa-let-7a-2_mature   5.70781e-07
    ## hsa-let-7a-3_mature            NA
    ## hsa-let-7b_mature     8.47481e-01
    ## hsa-let-7c_mature     1.79435e-01

``` r
ordered <- as.data.frame(res[order(res$padj),])
ordered <- ordered %>% rownames_to_column("miRNA")

# export results
write_csv(ordered, "output/all_samples/bg_20/PFAS_vs_Ctrl_statistical_results_bg20.csv")
```
