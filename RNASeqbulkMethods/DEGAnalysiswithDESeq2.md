---
title: "RNA-seq Differential Expression (DE) Analysis Using DESeq2"
author: "Jessica Randall"
date: "1/30/2020"
output:
  bookdown::html_document2:
    keep_md: yes
    theme: "paper"
    highlight: "default"
    toc: yes
---

<img src="eicc_logo.png" width="949" />

The original paper from Anders & Huber that introduces the concepts implemented in DESeq2 is available here https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)

We will be using the pasilla package for our example data which is cited here. Wolfgang Huber and Alejandro Reyes (2019). pasilla: Data package with per-exon and per-gene read counts of RNA-seq samples of Pasilla knock-down by Brooks et al., Genome Research 2011.. R package version 1.12.0.

Data for the package is part of this paper: "Conservation of an RNA regulatory map between Drosophila and mammals" by Brooks AN, Yang L, Duff MO, Hansen KD, Park JW, Dudoit S, Brenner SE, Graveley BR, Genome Res. 2011 Feb;21(2):193-202, Epub 2010 Oct 4, PMID: 20921232.

Since we will be using the apeglm shrinkage estimator to create our MA plot, I have also cited those authors
Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

Briefly, DESeq2 uses a Wald test to determine differential gene expression between at least two groups. DESeq2 assumes a negative binomial model which mathematically accounts for the fact that we are assessing gene counts and we are assuming that most genes we are comparing between the groups will not be differentially expressed.

Please note that DESeq2 has a number of capabilities that we will not be covering. You can import data for use in DESeq2 in several different ways, there are options for single cell projects, for incorporating Bayesian statistics, for time-series experiments, for outliers, for obtaining all of the results that DESeq2 functions produce, and even more options for graphing (including an R shiny app we plan to cover in a future walkthrough) so we strongly encourage you to reach out to EICC with questions regarding options available to you with DESeq2.

### Definition of terms {-}

###### Wald test: {-}

DESeq2 offers two options with how it determines differential expression, the Wald test and the Likelihood Ratio Test. Users of edgeR may be familiar with the latter. In DESeq2, the default test for pairwise comparison analysis is called the Wald Test and it is looking at your control and your experimental samples to see if the difference between them is equal to zero or not. Compared to edgeR it has been our experience that DESeq2 is more liberal in its calling of genes as significantly DE. If you are not sure which or how many DE genes you are expecting to find in your experiment, DESeq2 would likely the best analysis tool for your data.

###### Shrinkage estimator: {-}

Shrinkage estimators help us more accurately measure exactly how much the log transformation is condensing the data as it tends to do, especially when there are very small and very large values in one dataset together. Choosing the best estimator to assess this depends on your experimental aims and what you would like to do with your data but picking one is essential for accurately visualizing data analyzed in DESeq2.

###### Unadjusted p values vs adjusted p values/(FDR): {-} 

In DE analysis, a single p-value tells you how likely it is that a single gene is differentially expressed between at least two groups (ex: a control and a treatment group) due to some actual difference between the groups as opposed to random chance. False Discovery Rate (FDR) tells you how likely it is that all genes identified as DE are false positives. A FDR of 5% means that among all genes called DE, an average of 5% of those are truly not DE. DE genes are only considered significantly so if they meet the adjusted p value, not only the unadjusted p value.

### Load our example data {-}

Our very first step is to load the libraries we'll need to acsess the functions required for analysis and graphing.Please see http://bioconductor.org/ for information about initial installation and use of Bioconductor and its packages.

Our example data are from the pasilla package available on Bioconductor. The experiment studied RNAi knockdown of Pasilla, the Drosophila melanogaster ortholog of mammalian NOVA1 and NOVA2, on the transcriptome. Data are provided by NCBI Gene Expression Omnibus under accession numbers GSM461176 to GSM461181.

DESeq2 offers many options for importing count data and data about your samples. Here we will demonstrate importing the count matrix and sample data from the pasilla package.

Please reach out to EICC if you would like to compare 3 or more groups as this is a simplified example. It may also be the case you will need more than 6 samples per experimental group or that you may need to remove genes with average counts greater than 5, 10, 15, or even 20 for sufficient statistical power. Please see our PROPER walkthrough for an example of our of power and sample size analysis.


```r
require ("pacman")
p_load("readr", "dplyr", "knitr", "DESeq2", "vsn", "ggplot2", "pheatmap", "EnhancedVolcano")
theme_set(theme_classic())

countdata <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)

countdata <- as.matrix(read.csv(countdata, sep="\t", row.names="gene_id"))


sampledata <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)

sampledata <- read.csv(sampledata, row.names=1)
sampledata <- sampledata[,c("condition","type")]
```

Our data is almost ready to analyze but first we have to check that the column names of our count data are the same as the row names of our sample data. If this is not true, we cannot proceed. DESeq2 requires that the rownames of the sample data are the same as the column names of the count data since these are both referring to the samples and should be using the same names.


```r
rownames(sampledata) <- sub("fb", "", rownames(sampledata))

countdata <- countdata[, rownames(sampledata)]

all(rownames(sampledata) == colnames(countdata))
```

```
## [1] TRUE
```
Now that we know this is true, we can proceed.

### Preparation for Differential Expression Analysis {-}

In order to preform a pairwise comparison using the default negative binomial Wald test we need to specify some information about our data. In DESeq2 we must create a special object called a DESeqdataset object, here abbreviated as dds. This object takes in the countdata, the sample data, and the variable we would like to compare between the samples as inputs.

Next, we specify that the untreated group is our reference group to which we would like to compare our treated samples.


```r
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = sampledata,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "untreated")
```

If we had any additional data to add about the samples that we wanted to include in our analysis we would add it next but since this is a simplified example, we are only comparing treated and control samples without taking into account any additional information about them.

At this point we generate our first exploratory visualization, the principal components analysis plot. This will show us how your data cluster or how similar each sample is to others of the same group. There are percentages along the axes and the percentage on the x-axis tells us how much the differences between the samples is explained by them being treated or untreated. 

We start by transforming our data using the variance stabilizing transformation available from the vsn library (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010). This is similar to using a log2 transformation with normally distributed data with many very small values. VST accounts for the library size and normalization factors used in the mean-normalizing that DESeq2 uses on raw counts to make those many very small values easier to see on graphs. 


```r
vsd_dds <- vst(dds, blind=FALSE)

(exp_pca <- plotPCA(vsd_dds, intgroup=c("condition")))
```

![](DEGAnalysiswithDESeq2_files/figure-html/exploratory pca-1.png)<!-- -->

We would interpret this as the samples being clustered clearly by group and interpret the percentage on the x-axis as 58% of the variability between these samples is due to them being treated or untreated. The y-axis tells us how much variability between these samples is due to other factors in our model or if we have none, sources of variability we may not have accounted for like sex or ethnicity which are aften leading contributors of variability between samples and should be accounted for in experimental design if you wish to control for their effects. 

### Prefilter genes with low counts {-}

Typically we want to ignore genes that have counts of zero across all samples since these are adding statistical noise. We may also want to be more stringent and remove genes with rows that sum to 10, 20 or even 30 or less since these could also be contributing noise.

DESeq2 will filter genes it deems as low counts automatically based on the sum of the mean-normalized counts in each row. We'll see the criteria chosen when we view the results of our analysis later on. The results will tell you how many genes were removed and how many remain. If you would like to specify your own cut-offs for filtering or if you do not want DESeq2 to do any additional filtering, these are parameters that can be adjusted. 

If you choose to do DE analysis through EICC we typically rely on DESeq2's robust filtering since it tends to increase power to detect DE genes but we would customize this part of the analysis based on your data should you choose to do so.

### Testing {-}

DESeq2's analysis step is quite neat compared to the analysis steps of edgeR and baySeq. It will tell you the steps it is taking with your data and you have the option to ask for additional output and customizations. The default analysis explained here is the use of the Wald test. DESeq2 also offers the option of the Likelihood Ratio test. Both of these tests rely on the assumption that your count data follow a negative binomial distrubition which means that we assume that most counts are very low and that there are more non-DE genes between the groups than there are DE genes. Which test we choose will depend on your experimental design and properties of your data. Our summary presentations for clients typically include information on the model we chose, justification, and the null and alternative hypotheses of that model.

Finally, we create our results with the results function, use the summary function to see a tabular summary of them, and save them as an R dataframe for further manipulation. we check that the dataframe was created successfully by using the informal unit test of dimension with the expected number of rows and columns and telling the program to stop if the file does not have these dimensions. After this runs successfull we would typically export them as a .csv file for you.

The summary gives you information about the total number of genes with non-zero read counts, the adjusted p value it used to determine significantly DE genes (default is 0.1 but this can be changed in the DESeq function above), the exact number and total percentage of the up and down regulated DE genes, the presence of any outliers, and the removal of any additional genes with low counts.


```r
dds <- DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
res <- results(dds)

summary(res)
```

```
## 
## out of 12359 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 521, 4.2%
## LFC < 0 (down)     : 540, 4.4%
## outliers [1]       : 1, 0.0081%
## low counts [2]     : 4035, 33%
## (mean count < 7)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
res_df <- as.data.frame(res) %>%
  mutate(ID = as.factor((row.names(res))))
stopifnot(nrow(res_df) == 14599 & ncol(res_df) == 7)
```

### Visualizations {-}

We now perform additional data visualizations. Typically we provide PCA plots, heatmaps, and volcano plots and we would be happy to work with you to customize these for publication. Please see our Data Visualization menu for more options and examples from previous projects.

For our heatmaps, we typically show the mean normalized counts of the samples. Since we want to use the vsd transformed data from the results and it's currently in special Bioconductor format, we extract those transformed counts from that format and use it in the more common data frame format. 

Next, we put the samples in the order of treated to untreated and select only the 20 genes with the lowest adjusted p-values be displayed. It is possible to keep all DE genes on a heatmap but after 20 it becomes nearly impossible to read each of the individual gene names.

We use the pheatmap function in the pheatmap library to specify that we want to graph only the mean normalized counts of the genes with the 20 lowest adjusted p-values from the dataframe of the VST tranformed data using the sample names from the "condition" variable in our sample data so we can see which samples belong to which group.

```r
vsd_dds_df <- as.data.frame(assay(vsd_dds))

vsd_dds_df <- vsd_dds_df[c("treated1", "treated2", "treated3",
                           "untreated1", "untreated2", "untreated3", "untreated4")]

sel_padj <- order(res_df$padj, decreasing = FALSE)[1:20]

#select vars of interest
annotation <- as.data.frame(colData(dds)["condition"])

(ht <- pheatmap(vsd_dds_df[sel_padj, ],
                cluster_rows=FALSE,
                show_rownames = TRUE,
                cluster_cols = FALSE,
                annotation_col = annotation,
                width = 1))
```

![](DEGAnalysiswithDESeq2_files/figure-html/visualize-1.png)<!-- -->

Here we see the differences in mean-normalized counts between the samples in the treated vs untreated groups in the genes sorted by smalles adjusted p-value. Please note that hese are sorted for convenience but the gene at the top of the list is no more significant than the gene at the bottom of the list. As is the case with nominal p-values, a smaller adjusted p-value does not make a gene more statistically significant than one with a larger adjusted p-value. If the genes are below the threshold, they are all equally statistically significantly differentially expressed. These are sorted for convenience but the gene at the top of the list is no more significant than the gene at the bottom of the list. 


We also typically provide clients with an initial volcano plot created with the EnhancedVolcano R library. Similar to the PCA plot and heatmap, this is a highly customizable graph and we would like to work with you to design graphs which best tell the story of your results. 

A volcano plot is technically a scatter plot where the x-axis has the log2 transformed fold changes between the compared samples and the y axis has the local adjusted p-values for each gene, also called the q-value. Here we have also labelled the genes with FDR < 0.1 as that is where we set our threshold when we generated our results. The points in red are those which meet the threshold for statistical significance with a qvalue  less than or equal to 0.1 and a log2 fold change of 1.0 or greater. Points in green are those with only log2 fold changes >1.0 and those in blue have qvlaues < 0.1. The points in grey are non statistically significant by any measure. All of these parameters can be adjusted based on your cutoffs and thresholds. 


```r
(volc<- EnhancedVolcano(res,
                        lab = rownames(res),
                        x = "log2FoldChange",
                        y = "padj",
                        xlim = c(-6, 6),
                        title= NULL,
                        subtitle= "Log2 Fold Change vs q values",
                        FCcutoff = 1.0,
                        pLabellingCutoff = 0.1,
                        pCutoff = 0.1,
                        legendPosition = "bottom",
                        legend=c("NS", "Log2 fold-change", "adj P-value",
                                  "adj P-value & Log2 fold-change")))
```

![](DEGAnalysiswithDESeq2_files/figure-html/volcano-1.png)<!-- -->

There are many more functions and many more specifications to functions than are used here in order to show a simplified example of one of the tools we use for differential expression analysis. Obtaining specific, actionable, and publication quality results from analysis requires a deeper understanding of your specific data set and we would love the opportunity to discuss these options with you.

While we encourage clients to reach out prior to sequencing so that we can collaborate to design the experiment to answer your specific questions, we look forward to hearing from you at any stage of your RNAseq project. Please find our contact information available here https://www.cores.emory.edu/eicc/about/index.html

### Session information and References {-}


```r
sessionInfo()
```

```
## R version 3.6.2 (2019-12-12)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] EnhancedVolcano_1.4.0       ggrepel_0.8.1              
##  [3] pheatmap_1.0.12             ggplot2_3.2.1              
##  [5] vsn_3.54.0                  DESeq2_1.26.0              
##  [7] SummarizedExperiment_1.16.1 DelayedArray_0.12.2        
##  [9] BiocParallel_1.20.1         matrixStats_0.55.0         
## [11] Biobase_2.46.0              GenomicRanges_1.38.0       
## [13] GenomeInfoDb_1.22.0         IRanges_2.20.2             
## [15] S4Vectors_0.24.3            BiocGenerics_0.32.0        
## [17] dplyr_0.8.3                 readr_1.3.1                
## [19] pacman_0.5.1                knitr_1.27                 
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6           bit64_0.9-7            RColorBrewer_1.1-2    
##  [4] tools_3.6.2            backports_1.1.5        R6_2.4.1              
##  [7] affyio_1.56.0          rpart_4.1-15           Hmisc_4.3-0           
## [10] DBI_1.1.0              lazyeval_0.2.2         colorspace_1.4-1      
## [13] nnet_7.3-12            withr_2.1.2            tidyselect_1.0.0      
## [16] gridExtra_2.3          bit_1.1-15.1           compiler_3.6.2        
## [19] preprocessCore_1.48.0  htmlTable_1.13.3       labeling_0.3          
## [22] bookdown_0.17          scales_1.1.0           checkmate_1.9.4       
## [25] genefilter_1.68.0      affy_1.64.0            stringr_1.4.0         
## [28] digest_0.6.23          foreign_0.8-72         rmarkdown_2.1         
## [31] XVector_0.26.0         base64enc_0.1-3        jpeg_0.1-8.1          
## [34] pkgconfig_2.0.3        htmltools_0.4.0        limma_3.42.0          
## [37] htmlwidgets_1.5.1      rlang_0.4.4            rstudioapi_0.10       
## [40] RSQLite_2.2.0          farver_2.0.3           acepack_1.4.1         
## [43] RCurl_1.98-1.1         magrittr_1.5           GenomeInfoDbData_1.2.2
## [46] Formula_1.2-3          Matrix_1.2-18          Rcpp_1.0.3            
## [49] munsell_0.5.0          lifecycle_0.1.0        stringi_1.4.5         
## [52] yaml_2.2.0             zlibbioc_1.32.0        grid_3.6.2            
## [55] blob_1.2.1             crayon_1.3.4           lattice_0.20-38       
## [58] splines_3.6.2          annotate_1.64.0        hms_0.5.3             
## [61] locfit_1.5-9.1         pillar_1.4.3           geneplotter_1.64.0    
## [64] XML_3.99-0.3           glue_1.3.1             evaluate_0.14         
## [67] latticeExtra_0.6-29    data.table_1.12.8      BiocManager_1.30.10   
## [70] png_0.1-7              vctrs_0.2.2            gtable_0.3.0          
## [73] purrr_0.3.3            assertthat_0.2.1       xfun_0.12             
## [76] xtable_1.8-4           survival_3.1-8         tibble_2.1.3          
## [79] AnnotationDbi_1.48.0   memoise_1.1.0          cluster_2.1.0
```

```r
citation("readr")
```

```
## 
## To cite package 'readr' in publications use:
## 
##   Hadley Wickham, Jim Hester and Romain Francois (2018). readr: Read
##   Rectangular Text Data. R package version 1.3.1.
##   https://CRAN.R-project.org/package=readr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {readr: Read Rectangular Text Data},
##     author = {Hadley Wickham and Jim Hester and Romain Francois},
##     year = {2018},
##     note = {R package version 1.3.1},
##     url = {https://CRAN.R-project.org/package=readr},
##   }
```

```r
citation("dplyr")
```

```
## 
## To cite package 'dplyr' in publications use:
## 
##   Hadley Wickham, Romain François, Lionel Henry and Kirill Müller
##   (2019). dplyr: A Grammar of Data Manipulation. R package version
##   0.8.3. https://CRAN.R-project.org/package=dplyr
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {dplyr: A Grammar of Data Manipulation},
##     author = {Hadley Wickham and Romain François and Lionel Henry and Kirill Müller},
##     year = {2019},
##     note = {R package version 0.8.3},
##     url = {https://CRAN.R-project.org/package=dplyr},
##   }
```

```r
citation("knitr")
```

```
## 
## To cite the 'knitr' package in publications use:
## 
##   Yihui Xie (2020). knitr: A General-Purpose Package for Dynamic Report
##   Generation in R. R package version 1.27.
## 
##   Yihui Xie (2015) Dynamic Documents with R and knitr. 2nd edition.
##   Chapman and Hall/CRC. ISBN 978-1498716963
## 
##   Yihui Xie (2014) knitr: A Comprehensive Tool for Reproducible
##   Research in R. In Victoria Stodden, Friedrich Leisch and Roger D.
##   Peng, editors, Implementing Reproducible Computational Research.
##   Chapman and Hall/CRC. ISBN 978-1466561595
## 
## To see these entries in BibTeX format, use 'print(<citation>,
## bibtex=TRUE)', 'toBibtex(.)', or set
## 'options(citation.bibtex.max=999)'.
```

```r
citation("ggplot2")
```

```
## 
## To cite ggplot2 in publications, please use:
## 
##   H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
##   Springer-Verlag New York, 2016.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Book{,
##     author = {Hadley Wickham},
##     title = {ggplot2: Elegant Graphics for Data Analysis},
##     publisher = {Springer-Verlag New York},
##     year = {2016},
##     isbn = {978-3-319-24277-4},
##     url = {https://ggplot2.tidyverse.org},
##   }
```

```r
citation("DESeq2")
```

```
## 
##   Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change
##   and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
##   (2014)
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
##     author = {Michael I. Love and Wolfgang Huber and Simon Anders},
##     year = {2014},
##     journal = {Genome Biology},
##     doi = {10.1186/s13059-014-0550-8},
##     volume = {15},
##     issue = {12},
##     pages = {550},
##   }
```

```r
citation("vsn")
```

```
## 
## To cite the vsn package in publications use:
## 
##   Wolfgang Huber, Anja von Heydebreck, Holger Sueltmann, Annemarie
##   Poustka and Martin Vingron. Variance Stabilization Applied to
##   Microarray Data Calibration and to the Quantification of Differential
##   Expression. Bioinformatics 18, S96-S104 (2002).
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression},
##     author = {Wolfgang Huber and Anja {von Heydebreck} and Holger Sueltmann and Annemarie Poustka and Martin Vingron},
##     journal = {Bioinformatics},
##     year = {2002},
##     volume = {18 Suppl. 1},
##     pages = {S96-S104},
##   }
```

```r
citation("pheatmap")
```

```
## 
## To cite package 'pheatmap' in publications use:
## 
##   Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version
##   1.0.12. https://CRAN.R-project.org/package=pheatmap
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {pheatmap: Pretty Heatmaps},
##     author = {Raivo Kolde},
##     year = {2019},
##     note = {R package version 1.0.12},
##     url = {https://CRAN.R-project.org/package=pheatmap},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
```

```r
citation("EnhancedVolcano")
```

```
## 
## To cite package 'EnhancedVolcano' in publications use:
## 
##   Kevin Blighe, Sharmila Rana and Myles Lewis (2019). EnhancedVolcano:
##   Publication-ready volcano plots with enhanced colouring and labeling.
##   R package version 1.4.0.
##   https://github.com/kevinblighe/EnhancedVolcano
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {EnhancedVolcano: Publication-ready volcano plots with enhanced colouring and
## labeling},
##     author = {Kevin Blighe and Sharmila Rana and Myles Lewis},
##     year = {2019},
##     note = {R package version 1.4.0},
##     url = {https://github.com/kevinblighe/EnhancedVolcano},
##   }
## 
## ATTENTION: This citation information has been auto-generated from the
## package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
```






