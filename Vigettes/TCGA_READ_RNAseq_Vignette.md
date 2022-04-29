Removing library size and plate effects from the TCGA rectum
adenocarcinoma RNA-seq data using RUV-III-PRPS
================
true
15-02-2020

<style type="text/css">
h1.title {
  font-size: 28px;
  color: DarkRed;
}
h1 { /* Header 1 */
  font-size: 24px;
  color: DarkBlue;
}
h2 { /* Header 2 */
    font-size: 20px;
  color: DarkBlue;
}
h3 { /* Header 3 */
    font-size: 18px;
  color: DarkBlue;
}
h4 { /* Header 3 */
    font-size: 16px;
  color: DarkBlue;
}
</style>
<style>
p.caption {
  font-size: 46em;
  font-style: italic;
  color: black;
}
</style>

# Introduction

Effective removal of unwanted variation is essential to derive
meaningful biological results from RNA-seq data, particularly when the
data comes from large and complex study. We have previously proposed a
new method, removing unwanted variation III (RUV-III) to normalize gene
expression data [(R.Molania, NAR,
2019)](https://academic.oup.com/nar/article/47/12/6073/5494770?login=true).
The RUV-III method requires well-designed technical replicates
(well-distributed across sources of unwanted variation) and negative
control genes to estimate known and unknown sources of unwanted
variation and remove it from the data.  
We propose a novel strategy, pseudo-replicates of pseudo-samples (PRPS)
[R.Molania, bioRxiv,
2021](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1), for
deploying RUV-III to normalize RNA-seq data in situations when technical
replicate are not available or well-designed. Using the RUV-III with
PRPS is straightforward Our approach requires at least one **roughly**
known biologically homogenous subclass of samples presented across
sources of unwanted variation. For example, in a cancer RNA-seq study
where there are normal tissue samples present across all sources of
unwanted variation we can use these samples to create PRPS.  
To create PRPS, we first need to identify the sources of unwanted
variation, which we call batches in the data. Then the gene expression
measurements of biologically homogeneous sets of samples are averaged
within batches, and the results called pseudo-samples. Since the
variation between pseudo-samples in different batches is mainly unwanted
variation, by defining them as pseudo-replicates and used in RUV-III as
replicates, we can easily and effectively remove the unwanted variation.
We refer to our paper for more technical details [R.Molania, bioRxiv,
2021](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1).  

Here, we use the TCGA rectum adenocarcinoma (READ) RNA-seq data as an
example to show how to remove library size and batch effects (plate
effects) from the data. We illustrate the value of our approach by
comparing it to the standard TCGA normalizations on the TCGA READ data.
Further, we demonstrate how unwanted variation can compromise several
downstream analyses that can lead to wrong biological conclusions. We
will also assess the performance of RUV-III with poorly chosen PRPS and
in situations where biological labels are partially known.  
Note that RUV-III with PRPS is not limited to TCGA data: it can be used
for any large genomics project involving multiple labs, technicians,
platforms, …  

## Data preparation

The TCGA consortium aligned RNA sequencing reads to the hg38 reference
genome using the STAR aligner and quantified the results at gene level
using the HTseq and Gencode v22 gene-annotation
[Ref](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).
The TCGA RNA-seq data are publicly available in three formats: raw
counts, FPKM and FPKM with upper-quartile normalization (FPKM.UQ). All
these formats for individual cancer types (33 cancer types, \~ 11000
samples) were downloaded using the
*[TCGAbiolinks](https://bioconductor.org/packages/3.14/TCGAbiolinks)*
R/Bioconductor package (version 2.16.1). The TCGA normalized microarray
gene expression data were downloaded from the Broad GDAC
[Firehose](https://gdac.broadinstitute.org) repository , data version
2016/01/28. Tissue source sites (TSS), and batches of sequencing-plates
were extracted from individual TCGA [patient
barcodes](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/),
and sample processing times were downloaded from the [MD Anderson Cancer
Centre TCGA Batch Effects
website](https://bioinformatics.mdanderson.org/public-software/tcga-batch-effects).
Pathological features of cancer patients were downloaded from the Broad
GDAC Firehose repository (<https://gdac.broadinstitute.org>). The
details of processing the TCGA BRCA RNA-seq samples using two flow cell
chemistries were received by personal communication from Dr. K Hoadley.
The TCGA survival data reported by [Liu et
al.](https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418302290%3Fshowall%3Dtrue)
were used in this paper. The consensus measurement of purity estimation
(CPE) were downloaded from the [Aran et
al](https://www.nature.com/articles/ncomms9971) study.  
We have generated SummarizedExperiment objects for all the TCGA RNA-seq
datasets. These datasets can be found here
[TCGA\_PanCancerRNAseq](https://zenodo.org/record/6326542#.YimR0C8Rquo).
Unwanted variation of all the datasets can be explored using an Rshiny
application published in [(R.Molania, bioRxiv,
2021)](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1.article-metrics).  
All datasets that are required for this vignette can be found here
[link](https://doi.org/10.5281/zenodo.6392171)

# TCGA READ gene expression data

## RNA-seq data

We load the TCGA\_SummarizedExperiment\_HTseq\_READ.rds file. This is a
SummarizedExperiment object that contains:  
**assays:**  
-Raw counts  
-FPKM  
-FPKM.UQ  
**colData:**  
-Batch information  
-Clinical information (collected from different resources)  
**rowData:**  
-Genes’ details (GC, chromosome, …)  
-Several lists of housekeeping genes  

The lists of housekeeping genes might be suitable to use as negative
control genes (NCG) for RUV-III normalization.  
The Libraries\_HelperFunctions\_ForTcgaReadVignette.R containes all
helper functions that are required for this vignette.

``` r
# source('Libraries_HelperFunctions_ForTcgaReadVignette.R')
# read.se <- readRDS(
#   '../TCGA_SummarizedExperiment_HTseq_READ.rds'
#   )
# source('Scripts/Libraries_HelperFunctions_ForTcgaReadVignette.R')
# read.se <- readRDS(
#   'TCGA_SummarizedExperiment_HTseq_READ.rds'
#   )
```

## Microarray gene expression data

We also load the TCGA READ microarray gene expression data. This data
will be used as an orthogonal platform to assess the performance of
different RNA-seq normalizations. The data was downloaded from the TCGA
firehouse repositories. This data contains 17814 gene and 72 samples.

<!-- ```{r tcgaReadMicroArray, error=F, message=F, warning=F} -->
<!-- read.micro.array <- base::readLines( -->
<!--   '../READ.transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt') -->
<!-- read.micro.array <- read.delim( -->
<!--   file = base::textConnection(read.micro.array[-2]),   -->
<!--   row.names = 1, -->
<!--   stringsAsFactors = FALSE -->
<!--   ) -->
<!-- read.micro.array <- read.micro.array[ -->
<!--   complete.cases(read.micro.array) ,  -->
<!--   ] -->
<!-- colnames(read.micro.array) <- gsub( -->
<!--   '\\.',  -->
<!--   '-',  -->
<!--   colnames(read.micro.array) -->
<!--   ) -->
<!-- row.lables <- row.names(read.micro.array) -->
<!-- read.micro.array <- apply( -->
<!--   read.micro.array,  -->
<!--   2,  -->
<!--   as.numeric -->
<!--   ) -->
<!-- row.names(read.micro.array) <- row.lables -->
<!-- ### Sample annotation -->
<!-- index.cancer <- read.se$tissue == 'cancer' -->
<!-- read.micro.array.sampleAnnot <- as.data.frame( -->
<!--   SummarizedExperiment::colData(read.se[ , index.cancer]) -->
<!--   ) -->
<!-- common.samples <- intersect( -->
<!--   read.micro.array.sampleAnnot$Sample, -->
<!--   colnames(read.micro.array) -->
<!--   ) # 67 -->
<!-- read.micro.array.sampleAnnot <- read.micro.array.sampleAnnot[  -->
<!--   common.samples , ] -->
<!-- read.micro.array <- read.micro.array[ , common.samples] # 17812 genes * 67 samples -->
<!-- ``` -->
<!-- ## Removing genes and plates or samples -->
<!-- Here, we explain what kind of genes and plates or samples are removed before any down-stream analysis. -->
<!-- ### Lowly expressed genes -->
<!-- We identify lowly expressed genes in cancer and normal samples separately. Genes with at least 15 raw counts in at least 10% of cancer and normal samples are retained for down-stream analyses. These details can be found in the "keep.cancer" and 'keep.normal' columns of the gene annotation file in the SummarizedExperiment object. -->
<!-- ```{r LowlyExprGene,  warning=F, message=F, error=F, } -->
<!-- normal.tissues <- SummarizedExperiment::colData(read.se)$tissue == 'normal' -->
<!-- cancer.tissues <- SummarizedExperiment::colData(read.se)$tissue == 'cancer' -->
<!-- keep.genes.normal <- apply( -->
<!--   SummarizedExperiment::assay(read.se[ , normal.tissues], 'HTseq_counts'),  -->
<!--   1,  -->
<!--   function(x) length(x[x > 15]) >= round(.1*sum(normal.tissues), digits = 0) -->
<!--   ) -->
<!-- keep.genes.normal <- names(keep.genes.normal[keep.genes.normal == TRUE]) -->
<!-- keep.genes.cancer <- apply( -->
<!--   SummarizedExperiment::assay(read.se[ , cancer.tissues], 'HTseq_counts'),  -->
<!--   1,  -->
<!--   function(x) length(x[x > 15]) >= round(.1*sum(cancer.tissues), digits = 0) -->
<!--   ) -->
<!-- keep.genes.cancer <- names( -->
<!--   keep.genes.cancer[keep.genes.cancer == TRUE]) -->
<!-- ``` -->
<!-- ### Keep protein conding genes -->
<!-- We also keep only protein coding genes. This detail can be found in the "gene_type" column of the gene annotation file. This is a arbitrary filtering, any gene_type of interest can be retained in the data. -->
<!-- ```{r keepProteinCoding, message=FALSE, warning=F} -->
<!-- proteinCoding.index <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.se) -->
<!--   )$gene_type. == 'protein.coding' -->
<!-- keep.proteinCoding <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.se) -->
<!--   )$gene_id.v[proteinCoding.index] -->
<!-- selected.genes <- intersect( -->
<!--   unique(c(keep.genes.normal, keep.genes.cancer)), -->
<!--   keep.proteinCoding) -->
<!-- SummarizedExperiment::rowData(read.se)$selected.genes <- 'remove' -->
<!-- SummarizedExperiment::rowData(read.se)$selected.genes[ -->
<!--   SummarizedExperiment::rowData(read.se)$gene_id.v  -->
<!--   %in% selected.genes] <- 'keep' -->
<!-- ``` -->
<!-- ### Remove genes that have no or duplicated ENTREZ or gene symbol ids -->
<!-- Further, we remove genes without or with duplicated ENTREZ gene id or gene symbol. -->
<!-- ```{r Enterz, message=FALSE, warning=F} -->
<!-- ### entrez ids -->
<!-- SummarizedExperiment::rowData(read.se)$entrezgene.use <- 'keep' -->
<!-- na.duplicated <- -->
<!--   is.na(SummarizedExperiment::rowData(read.se)$entrezgene_id_BioMart) | -->
<!--   duplicated(SummarizedExperiment::rowData(read.se)$entrezgene_id_BioMart) -->
<!-- SummarizedExperiment::rowData(read.se)$entrezgene.use[na.duplicated] <- -->
<!--   'remove' -->
<!-- ### gene symbol -->
<!-- SummarizedExperiment::rowData(read.se)$geneName.use <- 'keep' -->
<!-- na.duplicated <- -->
<!--   is.na(SummarizedExperiment::rowData(read.se)$gene_name.) | -->
<!--   duplicated(SummarizedExperiment::rowData(read.se)$gene_name.) -->
<!-- SummarizedExperiment::rowData(read.se)$geneName.use[na.duplicated] <- -->
<!--   'remove' -->
<!-- keep.genes <- -->
<!--   SummarizedExperiment::rowData(read.se)$entrezgene.use == 'keep' & -->
<!--   SummarizedExperiment::rowData(read.se)$geneName.use == 'keep' -->
<!-- read.se <- read.se[keep.genes,] -->
<!-- ``` -->
<!-- ### Keep genes that are required for the CMS identification -->
<!-- Colorectal cancers are classified into four transcriptomic-based subtypes, consensus molecular subtypes [(CMS)](https://www.nature.com/articles/nm.3967), with distinct features. We also keep all genes that are used by the [CMScaller R package](https://www.nature.com/articles/s41598-017-16747-x) to identify the CMS in the data. -->
<!-- ```{r CMSgenes, message=FALSE, warning=F} -->
<!-- cms.geneSets <- CMScaller::geneSets.CMS -->
<!-- cms.genes <- unname(unlist(cms.geneSets)) -->
<!-- SummarizedExperiment::rowData(read.se)$cms.genes <- 'no' -->
<!-- index <- -->
<!--   SummarizedExperiment::rowData(read.se)$entrezgene_id_BioMart %in%  -->
<!--   unique(cms.genes) -->
<!-- SummarizedExperiment::rowData(read.se)$cms.genes[index] <- 'yes' -->
<!-- keep.genes <- -->
<!--   SummarizedExperiment::rowData(read.se)$selected.genes == 'keep' | -->
<!--   SummarizedExperiment::rowData(read.se)$cms.genes == 'yes' -->
<!-- read.se <- read.se[keep.genes,] -->
<!-- ``` -->
<!-- ### Keep plates with at least 2 samples -->
<!-- The READ RNA-seq study involved 177 assays generated using 14 plates over four years [(R.Molania, bioRxiv, 2021)](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1.article-metrics). We keep plates with at least three samples for down-stream analyses. The plate "A32Y" has only one sample, so this will be excluded from the analysis. -->
<!-- ```{r KeepPlates, message=FALSE, warning=F} -->
<!-- keep.plates <- names(which(table(read.se$plate_RNAseq) > 2)) -->
<!-- keep.plates <- read.se$plate_RNAseq %in% keep.plates -->
<!-- read.se <- read.se[ -->
<!--   SummarizedExperiment::rowData(read.se)$gene_id.v ,  -->
<!--   keep.plates] # 16327 176 -->
<!-- ``` -->
<!-- After the filtering above, the READ SummarizedExperiment object contains 16327 genes and 176 samples.  -->
<!-- ## Library size (sequencing-depth) -->
<!-- After removing genes and samples, we compute library size (total counts) and add this to the SummarizedExperiment object. We will use log2 of library size for all down-stream analyses. -->
<!-- ```{r LibSize, message=FALSE, warning=F} -->
<!-- read.se$libSize <-log2(colSums( -->
<!--     SummarizedExperiment::assay(read.se, 'HTseq_counts'))) -->
<!-- ``` -->
<!-- Figure \@ref(fig:lsPlots)  shows the library size of the TCGA READ RNA-seq data across years. Substantial library size differences between samples profiled in 2010 and the rest of the samples are clearly visible. -->
<!-- ```{r lsPlots, message=FALSE, warning=FALSE, size='small', fig.dim=c(7,3), fig.align = 'center', fig.cap='Library size of the TCGA READ RNA-seq data coloured  by different years.'} -->
<!-- df <- data.frame( -->
<!--   ls = read.se$libSize,  -->
<!--   samples = c(1:ncol(read.se)),  -->
<!--   time = read.se$year_mda -->
<!--   ) -->
<!-- ggplot(df, aes(x = samples, y = ls, color = time)) + -->
<!--   geom_point() + -->
<!--   scale_color_manual(values = years.colors, name = 'Time (Years)') + -->
<!--   ylab(expression(Log[2] ~ 'library size')) + -->
<!--   xlab('Samples') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 10), -->
<!--     axis.text.y = element_text(size = 10), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 18)) -->
<!-- ``` -->
<!-- Further, figure \@ref(fig:lsHk)  shows the library sizes of several housekeeping lists in the TCGA READ RNA-seq data. This result shows that the library size differences across samples are unwanted variation. -->
<!-- ```{r lsHk, message=FALSE, warning=FALSE, fig.dim=c(10,6), fig.cap='Library size of several housekeeping gene lists of the TCGA READ RNA-Seq data coloured by different years'} -->
<!-- hk.list <- c( -->
<!--     "sc.hk" , -->
<!--     "rnaseq.hk", -->
<!--     "array.hk", -->
<!--     "nanostring.hk", -->
<!--     "sinscore.hk") -->
<!-- ls.hk <- lapply( -->
<!--   hk.list,  -->
<!--   function(x){ -->
<!--     hk.index <- SummarizedExperiment::rowData(read.se)[ , x] == 'yes' -->
<!--     log2(colSums( -->
<!--     SummarizedExperiment::assay( -->
<!--       read.se[hk.index , ],  -->
<!--       'HTseq_counts') -->
<!--     )) -->
<!--   }) -->
<!-- names(ls.hk) <- -->
<!--   c( -->
<!--     'scRNASeq HK', -->
<!--     'RNASeq HK', -->
<!--     'Microarray HK', -->
<!--     'Nanostring PanCancer HK', -->
<!--     'singscore PanCancer HK') -->
<!-- ls.hk <- do.call(cbind, ls.hk) %>% -->
<!--   data.frame() %>% -->
<!--   dplyr::mutate(samples = factor( -->
<!--     x = c(1:ncol(read.se)),  -->
<!--     levels = c(1:ncol(read.se)))) %>% -->
<!--   dplyr::mutate(time = read.se$year_mda) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(samples, time),  -->
<!--     values_to = 'ls',  -->
<!--     names_to = 'hk') %>% -->
<!--   dplyr::mutate(hk = gsub('\\.', ' ', hk)) %>%   -->
<!--   dplyr::mutate(hk = factor( -->
<!--     hk,  -->
<!--     levels = c( -->
<!--       'RNASeq HK',  -->
<!--       'Microarray HK',  -->
<!--       'singscore PanCancer HK',  -->
<!--       'Nanostring PanCancer HK',  -->
<!--       'scRNASeq HK'))) %>%  -->
<!--   data.frame() -->
<!-- ggplot(ls.hk, aes(x = samples, y = ls, color = time)) + -->
<!--   geom_point() + -->
<!--   ylab(expression(Log[2]~'library size')) + -->
<!--   scale_color_manual(values = years.colors, name = 'Time (Years)') + -->
<!--   xlab('Samples') + -->
<!--   facet_wrap(~ hk, scale = 'free') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 0), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 12), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 14)) -->
<!-- ``` -->
<!-- ## Major time interval -->
<!-- Based on library size variation, we divide samples into two major time intervals, 2010 and 2011:2014, for down-stream analyses. We explored a range of clinical details of samples from the two major time intervals and did not find any specific markers that are associated with the time intervals. -->
<!-- ```{r MajorTimeInterval, message=FALSE, warning=F} -->
<!-- read.se$time.points <- '2010' -->
<!-- index <- read.se$year_mda != 2010 -->
<!-- read.se$time.points[index] <- '2011:2014' -->
<!-- read.se$time.points <- factor( -->
<!--   read.se$time.points , -->
<!--   levels = c( -->
<!--     '2010', -->
<!--     '2011:2014')) -->
<!-- ``` -->
<!-- ## Major gene expression-based biological populations -->
<!-- We identify major gene expression-based biological populations in order to create pseudo-samples [(R.Molania, bioRxiv, 2021)](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1.article-metrics). In the TCGA READ RNA-seq data, we use the microsatellite instability (MSI) and consensus molecular subtypes (CMS) to create different sets of PRPS.\ -->
<!-- Note that, any biological populations should roughly show distinct gene expression profile, otherwise they are not helpful for PRPS. For example, we did not find distinct gene expression patterns associated with tumor stage in the TCGA READ RNA-seq data, so we do not use tumor stage for PRPS. -->
<!-- ### Microsatellite instability (MSI) -->
<!-- The details of the MSI status can be found in the sample annotation file. There are:\ -->
<!-- msi-h: microsatellite instability high\ -->
<!-- msi-l: microsatellite instability low\ -->
<!-- mss: microsatellite stable\ -->
<!-- indeterminate -->
<!-- Figure \@ref(fig:MsiGroups) shows the numbers of individual MSI status in the TCGA READ RNA-Seq data. -->
<!-- ```{r MsiGroups, cache=T, message=F, warning=F, fig.dim=c(6,4), fig.align = 'center', fig.cap='MSI status in the TCGA READ-RNA-Seq data.'} -->
<!-- colnames(SummarizedExperiment::colData(read.se))[932] <- -->
<!--   'msi.status' -->
<!-- msi.groups <- sort(unique( -->
<!--   SummarizedExperiment::colData( -->
<!--     read.se)$msi.status)) -->
<!-- msi.new.names <- c( -->
<!--   'Indeterminate', -->
<!--   'MSI-H', -->
<!--   'MSI-L', -->
<!--   'MSS') -->
<!-- for(i in 1:4){ -->
<!--   index <- SummarizedExperiment::colData( -->
<!--     read.se)$msi.status == msi.groups[i] -->
<!--   SummarizedExperiment::colData( -->
<!--     read.se)$msi.status[index] <-  msi.new.names[i] -->
<!-- } -->
<!-- index <- read.se$tissue == 'normal' -->
<!-- read.se$msi.status[index] <- 'Adjacent normal' -->
<!-- read.se$msi.status <- factor( -->
<!--   x = read.se$msi.status, -->
<!--   levels = c( -->
<!--     'MSS', -->
<!--     'MSI-L', -->
<!--     'MSI-H', -->
<!--     'Indeterminate', -->
<!--     'Adjacent normal')) -->
<!-- ### plot -->
<!-- ggplot(data = as.data.frame(SummarizedExperiment::colData(read.se)), -->
<!--        aes(x = msi.status)) + -->
<!--   geom_bar() + -->
<!--   xlab('MSI') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 14, angle = 20, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 14)) -->
<!-- ``` -->
<!-- ### Consensus molecular subtypes (CMS) -->
<!-- As we mentioned above, colorectal cancers can be classified into four widely accepted consensus molecular subtypes (CMS) based on their gene expression profiles. This classification provides a framework for stratifying the treatment of patients with colon and rectum cancer. There are CMS1, CMS2, CMS3, CMS4 subtypes.\ -->
<!-- Here, we use the [CMScaller R package](https://www.nature.com/articles/s41598-017-16747-x) to identify the CMS of the TCGA READ RNA-seq samples. The CMScaller provides a classification based on pre-defined cancer-cell intrinsic CMS templates. Figure \@ref(fig:cmsTCGA) shows heatmaps of the relative expression levels of the CMS marker genes (vertical bar) with classifications indicated below (horizontal bar, white indicating prediction confidence p-values). -->
<!-- ```{r cmsTCGA, cache=TRUE, message=FALSE, warning=F, fig.cap='Consensus molecular subtypes (CMS) identification using the CMScaller R package in the TCGA raw counst, FPKM and FPKM.UQ data (from left to right). Heatmaps show the relative expression levels of CMS marker genes (vertical bar) with classifications indicated below (horizontal bar, white indicating prediction confidence p-values).'} -->
<!-- tcga.harmonized <- names(SummarizedExperiment::assays(read.se)) -->
<!-- index.cancer <- read.se$tissue == 'cancer' -->
<!-- set.seed(2010301149) -->
<!-- par(mfrow = c(1,3)) -->
<!-- cms.clusters.cancer.tcga <- lapply( -->
<!--   tcga.harmonized, -->
<!--   function(x){ -->
<!--     expr.matrix <- as.matrix(SummarizedExperiment::assay(read.se[ , index.cancer], x)) -->
<!--     row.names(expr.matrix) <- SummarizedExperiment::rowData(read.se)$entrezgene_id_BioMart -->
<!--     cms.cluster <- CMScaller::CMScaller( -->
<!--       emat = log2(expr.matrix + .5), -->
<!--       RNAseq = FALSE, -->
<!--       FDR = 0.05, -->
<!--       verbose = FALSE) -->
<!--     return(cms.cluster) -->
<!--     }) -->
<!-- ``` -->
<!-- We add the obtained CMS subtypes to the sample annotation object. -->
<!-- ```{r CmsOnCancerSamples, warning=FALSE, message=FALSE} -->
<!-- names(cms.clusters.cancer.tcga) <- tcga.harmonized -->
<!-- col.names <- paste0( -->
<!--   'cms.cancer.', -->
<!--   c('rawCounts', -->
<!--     'fpkm', -->
<!--     'fpkmUq')) -->
<!-- for(i in 1:3){ -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]] <- -->
<!--     'Adjacent normal' -->
<!--   index <- match( -->
<!--     row.names(cms.clusters.cancer.tcga[[i]]), -->
<!--     read.se$Sample -->
<!--   ) -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]][index] <- -->
<!--     as.character(cms.clusters.cancer.tcga[[i]]$prediction) -->
<!--   index <- is.na(SummarizedExperiment::colData( -->
<!--     read.se)[ , col.names[i]] -->
<!--     ) -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]][index] <- -->
<!--     'Not classified' -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]] <- factor( -->
<!--     x = SummarizedExperiment::colData(read.se)[ , col.names[i]], -->
<!--     levels = c( -->
<!--       'CMS1', -->
<!--       'CMS2', -->
<!--       'CMS3', -->
<!--       'CMS4', -->
<!--       'Adjacent normal', -->
<!--       'Not classified')) -->
<!--   } -->
<!-- ``` -->
<!-- We also apply the classifier to samples within the major time intervals (2010 and 2011:2014) using the raw counts, FPKM and FPK.UQ normalized datasets (Figure \@ref(fig:CmsWithinTimes)). The reason for applying the classifier within each key time interval was to assess the effect of large library differences on the CMS classifications. -->
<!-- ```{r CmsWithinTimes, cache=TRUE, message=FALSE, warning=F, fig.dim=c(8,6), fig.cap='Identification of consensus molecular subtypes (CMS) in the TCGA READ RNA-seq studies. Heatmaps show the relative expression levels of CMS marker genes (vertical bar) with classifications indicated below (horizontal bar, white indicating prediction confidence p-values). The CMS classification were performed within each key time intervals (first row is2010 and the second row is 2010-2014) to assess the impact of batch effects on the classification.'} -->
<!-- set.seed(2010301149) -->
<!-- par(mfrow = c(2,3)) -->
<!-- read.cancer.se <- read.se[ , index.cancer] -->
<!-- cms.clusters.time.points.cancer.tcga <- lapply( -->
<!--   levels(read.se$time.points), -->
<!--   function(x){ -->
<!--     index.time.points <- read.cancer.se$time.points == x -->
<!--     cms.clusters.cancer.tcga <- lapply( -->
<!--       tcga.harmonized, -->
<!--       function(y){ -->
<!--         expr.matrix <- as.matrix(SummarizedExperiment::assay( -->
<!--           read.cancer.se[, index.time.points], y)) -->
<!--         row.names(expr.matrix) <- SummarizedExperiment::rowData( -->
<!--           read.cancer.se)$entrezgene_id_BioMart -->
<!--         cms.cluster.time <- CMScaller::CMScaller( -->
<!--           emat = log2(expr.matrix + .5), -->
<!--           RNAseq = FALSE, -->
<!--           FDR = 0.05, -->
<!--           verbose = FALSE) -->
<!--         return(cms.cluster.time) -->
<!--     }) -->
<!--     names(cms.clusters.cancer.tcga) <- tcga.harmonized -->
<!--     return(cms.clusters.cancer.tcga) -->
<!--   }) -->
<!-- names(cms.clusters.time.points.cancer.tcga) <- paste0( -->
<!--   'CMS', -->
<!--   '_', -->
<!--   levels(read.se$time.points)) -->
<!-- ``` -->
<!-- Here we add the CMS details to the sample annotation object. -->
<!-- ```{r  cache=TRUE, message=FALSE, warning=F } -->
<!-- raw.counts <- as.data.frame(rbind( -->
<!--   cms.clusters.time.points.cancer.tcga$CMS_2010[[1]], -->
<!--   cms.clusters.time.points.cancer.tcga$`CMS_2011:2014`[[1]]) -->
<!--   ) -->
<!-- fpkm <- as.data.frame(rbind( -->
<!--   cms.clusters.time.points.cancer.tcga$CMS_2010[[2]], -->
<!--   cms.clusters.time.points.cancer.tcga$`CMS_2011:2014`[[2]]) -->
<!--   ) -->
<!-- fpkm.uq <- as.data.frame(rbind( -->
<!--   cms.clusters.time.points.cancer.tcga$CMS_2010[[3]], -->
<!--   cms.clusters.time.points.cancer.tcga$`CMS_2011:2014`[[3]]) -->
<!--   ) -->
<!-- cms.clusters.time.points.cancer.tcga <- list( -->
<!--   raw.counts = raw.counts, -->
<!--   fpkm = fpkm, -->
<!--   fpkm.uq = fpkm.uq) -->
<!-- col.names <- paste0( -->
<!--   'cms.cancer.time.points.', -->
<!--   c('rawCounts', -->
<!--     'fpkm', -->
<!--     'fpkmUq') -->
<!--   ) -->
<!-- for(i in 1:3){ -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]] <- 'Adjacent normal' -->
<!--   index <- match( -->
<!--     row.names(cms.clusters.time.points.cancer.tcga[[i]]), -->
<!--     read.se$Sample -->
<!--   ) -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]][index] <- -->
<!--     as.character(cms.clusters.time.points.cancer.tcga[[i]]$prediction) -->
<!--   index <- is.na(SummarizedExperiment::colData(read.se)[ , col.names[i]]) -->
<!--   SummarizedExperiment::colData(read.se)[ , col.names[i]][index] <- 'Not classified' -->
<!--  SummarizedExperiment::colData(read.se)[ , col.names[i]] <- factor( -->
<!--     x = SummarizedExperiment::colData(read.se)[ , col.names[i]], -->
<!--     levels = c('CMS1', -->
<!--                'CMS2', -->
<!--                'CMS3', -->
<!--                'CMS4', -->
<!--                'Adjacent normal', -->
<!--                'Not classified')) -->
<!-- } -->
<!-- read.cancer.se <- read.se[ , index.cancer] -->
<!-- ``` -->
<!-- # Study outline -->
<!-- Figure \@ref(fig:studyOutline) shows the outline of the TCGA READ RNA-seq study. This study involved 176 assays generated using 14 plates over four years. -->
<!-- ```{r studyOutline, warning=F, message=F, error=F, fig.cap='Outline of the TCGA rectum adenocarcinoma RNA-seq study. 176 rectum adenocarcinoma and adjacent normal tissues were collected from 13 tissue source sites (TSS) and distributed across 14 sequencing plates for profiling at 14 time points over a span of 4 years. The consensus molecular subtypes were obtained using the R package CMScaller on the FPKM.UQ normalized data. The MSI status were obtained from the pathological reports in the TCGA clinical data. The library sizes are calculated after removing lowly expressed genes and log2 transformed. The tumour purity scores are obtained from 1- stromal&immune scores.'} -->
<!-- selected.columns <- c( -->
<!--   'year_mda', -->
<!--   'plate_RNAseq', -->
<!--   'tss_RNAseq', -->
<!--   'libSize', -->
<!--   'purity_HTseq_FPKM', -->
<!--   'tissue', -->
<!--   'cms.cancer.fpkmUq', -->
<!--   'msi.status' -->
<!--   ) -->
<!-- sample.info <- as.data.frame( -->
<!--   SummarizedExperiment::colData(read.se)[ , selected.columns]) -->
<!-- sample.info$plate_RNAseq <- factor( -->
<!--   sample.info$plate_RNAseq, -->
<!--   levels = unique(sample.info$plate_RNAseq) -->
<!--   ) -->
<!-- sample.info$tss_RNAseq <- factor( -->
<!--   sample.info$tss_RNAseq, -->
<!--   levels = unique(sample.info$tss_RNAseq) -->
<!--   ) -->
<!-- ## Time (years) -->
<!-- H.time <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$year_mda), -->
<!--   cluster_columns  = F, -->
<!--   col =  years.colors, -->
<!--   name = 'Time (years)', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     ncol = 2, -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## Plates -->
<!-- colfunc <- grDevices::colorRampPalette( -->
<!--   RColorBrewer::brewer.pal(n = 11, name = 'PRGn')[-6]) -->
<!-- color.plates <- colfunc(length(unique(sample.info$plate_RNAseq))) -->
<!-- names(color.plates) <- levels(sample.info$plate_RNAseq) -->
<!-- H.plate <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$plate_RNAseq), -->
<!--   cluster_rows = FALSE, -->
<!--   cluster_columns = FALSE, -->
<!--   col = color.plates, -->
<!--   name = 'Plates', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     ncol = 4, -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## TSS -->
<!-- colfunc <- grDevices::colorRampPalette( -->
<!--   RColorBrewer::brewer.pal(n = 11, name = 'BrBG')[-6]) -->
<!-- color.tss <- colfunc(length(unique(sample.info$tss_RNAseq))) -->
<!-- names(color.tss) <- levels(sample.info$tss_RNAseq) -->
<!-- H.tss <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$tss_RNAseq), -->
<!--   cluster_rows = FALSE, -->
<!--   cluster_columns = FALSE, -->
<!--   col = color.tss, -->
<!--   name = 'Tissue source sites', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     ncol = 4, -->
<!--     title_gp = grid::gpar(fontsize = 14) -->
<!--     )) -->
<!-- ## Tissue -->
<!-- H.tissue <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$tissue), -->
<!--   cluster_rows = FALSE, -->
<!--   col = RColorBrewer::brewer.pal(9, 'Greys')[c(8,3)], -->
<!--   name = 'Tissues', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     direction = "vertical", -->
<!--     ncol = 1, -->
<!--     title_gp = grid::gpar(fontsize = 14), -->
<!--     labels = c('Primary tumor', 'Normal tissue'))) -->
<!-- ## CMS -->
<!-- H.cms <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$cms.cancer.fpkmUq), -->
<!--   cluster_rows = FALSE, -->
<!--   col = cms.colors.so, -->
<!--   name = 'CMS', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     direction = "vertical", -->
<!--     ncol = 1, -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## MSI -->
<!-- H.msi <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$msi.status), -->
<!--   cluster_rows = FALSE, -->
<!--   col = msi.colors.so, -->
<!--   name = 'MSI', -->
<!--   heatmap_legend_param = list( -->
<!--     color_bar = "discrete" , -->
<!--     direction = "vertical", -->
<!--     ncol = 1, -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## Purity -->
<!-- H.purity <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$purity_HTseq_FPKM), -->
<!--   cluster_rows = FALSE, -->
<!--   name = 'Tumor purity scores', -->
<!--   col = viridis::plasma(n = 10), -->
<!--   heatmap_legend_param = list( -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## Library size -->
<!-- H.ls <- ComplexHeatmap::Heatmap( -->
<!--   rev(sample.info$libSize), -->
<!--   cluster_rows = FALSE, -->
<!--   name = "Log2 library size", -->
<!--   col = viridis::viridis(n = 10), -->
<!--    heatmap_legend_param = list( -->
<!--     title_gp = grid::gpar(fontsize = 14))) -->
<!-- ## All -->
<!-- ComplexHeatmap::draw( -->
<!--   H.time + -->
<!--     H.plate + -->
<!--     H.tss + -->
<!--     H.tissue + -->
<!--     H.cms + -->
<!--     H.msi + -->
<!--     H.ls + -->
<!--     H.purity, -->
<!--   merge_legends = FALSE, -->
<!--   heatmap_legend_side = 'right') -->
<!-- ``` -->
<!-- # RUV-III normalization -->
<!-- Here, we will explain how to select a suitable set of negative control genes and pseudo-replicates of pseudo-samples to remove library size and plate effects from the TCGA READ RNA-seq data. -->
<!-- ## Selection of negative control genes (NCG) -->
<!-- First, we need to highlight that in our usage, a negative control gene is one that is not expected to change much across biological factors of interests [(R.Molania, NAR, 2019)](https://academic.oup.com/nar/article/47/12/6073/5494770?login=true). Second, in the presence of unwanted variation usually a subset of genes are affected in different ways, so we need to emphasize that our approach to negative controls is pragmatic: if using a given gene in RUV-III as one of the set of negative control genes helps, as indicated by various measures, then whether or not it is an ideal negative control gene is moot: it helped. This does not rule out the fact that we may be able to do better by replacing that gene with a different gene designated a negative control. Third, we point out that theoretical analyses not presented here show that it is not the extent to which individual genes designated as negative controls are ideal or less than ideal negative controls that drives the success or otherwise of RUV-III; that is a property of the full set of negative controls. We will frequently get very good results using the entire set of genes being studied as negative controls, even when many genes are changing. (That this is not unreasonable follow from the theory just mentioned.) Fourth, it is usually the case that using more genes as negative control genes is better than fewer. That is a matter of stability, but as stated in the introduction, there is a bias-variance trade-off here: using too many genes as negative controls may be counter-productive. You must look and see. Fifth, endogenous genes generally make more suitable negative controls than spike-ins. The reason here is the obvious one, namely, that endogenous genes have shared the complete sample experience of the other genes, whereas spike-ins can only reflect unwanted variation in the process from the point at which they were added onward. The best source of endogenous negative control genes are ones that were found to be stable in previous studies similar to the one being analyzed.\ -->
<!-- Here, we explain an approach the we used in our paper [R.Molania, bioRxiv, 2021](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1) to select a suitable set of negative control genes for the TCGA READ RNA-seq data.\ -->
<!-- First, we select samples that they have the same CMS subtypes obtained by applying the CMS classifier within and between the major time interval using the TCGA FPKM.UQ data. We found 118 samples that meet this criteria.\ -->
<!-- ```{r ConCms, message=FALSE, warning=F} -->
<!-- SummarizedExperiment::colData(read.cancer.se)$cms.use <- 'no' -->
<!-- a <- as.character( -->
<!--   SummarizedExperiment::colData(read.cancer.se)$cms.cancer.fpkmUq -->
<!--   ) -->
<!-- b <- as.character( -->
<!--   SummarizedExperiment::colData(read.cancer.se)$cms.cancer.time.points.fpkmUq -->
<!--   ) -->
<!-- SummarizedExperiment::colData(read.cancer.se)$cms.use[a==b] <- 'yes' -->
<!-- ### Data and gene annot -->
<!-- index.sample.use <- SummarizedExperiment::colData( -->
<!--   read.cancer.se)$cms.cancer.fpkmUq!='Not classified' & -->
<!--   SummarizedExperiment::colData(read.cancer.se)$cms.use == 'yes' -->
<!-- ``` -->
<!-- Then, we apply ANOVA on individual genes expression with the CMS being a factor using the TCGA FPKM.UQ data.\ -->
<!-- ```{r NcgAnovaCms, message=FALSE, warning=F} -->
<!-- n.cores <- 5 -->
<!-- ftest.cms <- .Ftest( -->
<!--   data = as.matrix(SummarizedExperiment::assay( -->
<!--     read.cancer.se[ , index.sample.use], -->
<!--     'HTseq_FPKM.UQ') -->
<!--     ), -->
<!--   variable = droplevels( -->
<!--     SummarizedExperiment::colData( -->
<!--       read.cancer.se -->
<!--       )$cms.cancer.fpkmUq[index.sample.use]), -->
<!--   is.log = FALSE, -->
<!--   n.cores = n.cores) -->
<!-- ``` -->
<!-- We rank genes based on their F-statistics and select genes that have the rank below 1000. These genes most likely do not capture the CMS variation. In addition, we remove genes that are on Y chromosome. Finally, we end up with 997 genes as a potential set of negative control genes for RUV-III normalization. -->
<!-- ```{r NcgSelection, message=FALSE, warning=F} -->
<!-- SummarizedExperiment::rowData(read.cancer.se)$cmsDE.genes <- rank( -->
<!--   ftest.cms$FValue) -->
<!-- negative.control.genes <-  SummarizedExperiment::rowData( -->
<!--   read.cancer.se -->
<!--   )$cmsDE.genes < 1000 & -->
<!--   SummarizedExperiment::rowData( -->
<!--     read.cancer.se -->
<!--     )$chromosome_name_BioMart != 'Y' -->
<!-- ``` -->
<!-- ### Assessments of negative control genes -->
<!-- We perform PCA on the raw counts of the TCGA READ RNA-seq data using only the selected negative control genes to assess their performance (Figure \@ref(fig:PcaOnNcg)). Ideally, they should capture the library size differences, but not much of variation related to the CMS. Note that, the RUV-III method is generally robust to negative control genes, but not always. Figure \@ref(fig:PcaOnNcg) shows that the selected negative control genes capture the library size differences and they do not capture variation related to the CMS. -->
<!-- ```{r PcaOnNcg, warning=F, message=F, fig.dim=c(11,6), fig.cap='PCA plots of the raw counts of the TCGA READ RNA-seq data using only the negative control genes (997 genes).'} -->
<!-- pca.ncg <- .pca( -->
<!--   data = as.matrix( -->
<!--     SummarizedExperiment::assay( -->
<!--       read.cancer.se[negative.control.genes , ], -->
<!--       'HTseq_counts') -->
<!--     ), -->
<!--   is.log = FALSE -->
<!--   ) -->
<!-- p1 <- .scatter.density.pc( -->
<!--   pcs = pca.ncg$sing.val$u[, 1:3], -->
<!--   pc.var = pca.ncg$var, -->
<!--   group.name = 'Time (years)', -->
<!--   group = SummarizedExperiment::colData(read.cancer.se)$time.points, -->
<!--   color = major.times.colors, -->
<!--   strokeSize = .2, -->
<!--   pointSize = 3, -->
<!--   strokeColor = 'gray30', -->
<!--   alpha = .5 -->
<!--   ) -->
<!-- p2 <- .scatter.density.pc( -->
<!--   pcs = pca.ncg$sing.val$u[, 1:3], -->
<!--   pc.var = pca.ncg$var, -->
<!--   group.name = 'CMS', -->
<!--   group = SummarizedExperiment::colData(read.cancer.se)$cms.cancer.fpkmUq, -->
<!--   color = cms.colors, -->
<!--   strokeSize = .2, -->
<!--   pointSize = 3, -->
<!--   strokeColor = 'gray30', -->
<!--   alpha = .5 -->
<!--   ) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(p1, p2, ncol = 4)) -->
<!-- ``` -->
<!-- ## Pseudo-replicates of pseudo-samples (PRPS) -->
<!-- As we have mentioned above, the RUV-III method also requires technical replicates to estimate one aspect of unwanted variation.  -->
<!-- Here we propose a new approach, pseudo-replicates of pseudo-samples (PRPS), for deploying RUV-III method to remove unwanted variation from the TCGA READ RNA-seq data. This approach requires at least one roughly known biologically homogeneous subclass of samples shared across the sources of unwanted variation. To create PRPS, we first need to identify the sources of unwanted variation, which we will call batches in the data. Then the gene expression measurements of biologically homogeneous sets of samples are averaged within batches, and the results called pseudo-samples. Since the variation between pseudo-samples in different batches is mainly unwanted variation, by defining them as pseudo-replicates and used in RUV-III as replicates, we can easily and effectively remove the unwanted variation. -->
<!-- ### Selection of biological populations to create PRPS -->
<!-- To select homogeneous biological populations, we consider two major biological factors:\ -->
<!-- 1. The CMS subtypes that were obtained by applying the CMS classifier across and within the major time intervals using the TCGA FPKM.UQ data.\ -->
<!-- 2. The MSI status.\ -->
<!-- The 11 combinations (we do not have CMS4_MSI-H) of the 4 CMS, and the 3 MSI statuses were considered to be homogeneous biological populations for the purpose of creating PRPS.\ -->
<!-- In general, in the TCGA RNA-seq data, plates are completely confounded with times, making it difficult to distinguish plate effects from time effects. We consider plates as batches to create pseudo-samples. Note that, we need to make sure that our pseudo-samples span the major time intervals, otherwise we will not be able to remove the library size differences in the data. -->
<!-- ```{r PrPsGeneration, message=FALSE, warning=FALSE, results=F} -->
<!-- samples.to.use <- -->
<!--   read.cancer.se$cms.cancer.fpkmUq != 'Not classified' & -->
<!--   read.cancer.se$msi.status != 'Indeterminate' & -->
<!--   read.cancer.se$cms.use == 'yes' -->
<!-- sample.info <- droplevels( -->
<!--   as.data.frame(SummarizedExperiment::colData(read.cancer.se[ , samples.to.use])) -->
<!--   ) -->
<!-- raw.counts <- as.data.frame( -->
<!--   SummarizedExperiment::assay(read.cancer.se[ , samples.to.use] , 'HTseq_counts') -->
<!--   ) -->
<!-- read.prps <- -->
<!--   .CreatePseudoSamplesForLsPurityBatch( -->
<!--     expr.data = as.data.frame( -->
<!--       SummarizedExperiment::assay(read.cancer.se[, samples.to.use] , 'HTseq_counts') -->
<!--     ), -->
<!--     sample.info = droplevels(as.data.frame( -->
<!--       SummarizedExperiment::colData(read.cancer.se[, samples.to.use]) -->
<!--     )), -->
<!--     batch = 'PlateId_mda', -->
<!--     biology = c('cms.cancer.fpkmUq', 'msi.status'), -->
<!--     purity = FALSE, -->
<!--     include.ls = FALSE, -->
<!--     include.purity = FALSE, -->
<!--     minSamplesPerBatchPS = 2) -->
<!-- ``` -->
<!-- ### PRPS map -->
<!-- The plot \@ref(fig:PrPsMap) shows the distribution of the homogeneous biological populations across plates in the data. To Create PS, we average gene expression of at least two samples with respect to the biological populations and plates. There 3 plates that we do not have any PS for them. -->
<!-- ```{r PrPsMap, message=FALSE, warning=FALSE, fig.cap='Plot showing the sample sizes of the major biological groups across plates in the TCGA READ RNA-seq data.'} -->
<!-- new.info <- droplevels( -->
<!--   as.data.frame( -->
<!--     SummarizedExperiment::colData( -->
<!--       read.cancer.se[ , samples.to.use])) -->
<!--   ) -->
<!-- new.info$new.batch <- paste0( -->
<!--   new.info$year_mda, -->
<!--   '_', -->
<!--   new.info$PlateId_mda -->
<!--   ) -->
<!-- new.info$biololy <- paste0( -->
<!--   new.info$cms.cancer.fpkmUq, -->
<!--   '_', -->
<!--   new.info$msi.status) -->
<!-- df_count <- new.info %>% -->
<!--     dplyr::count(new.batch, biololy) -->
<!-- df_count$use <- 'Un-selected' -->
<!-- df_count$use[df_count$n > 1] <- 'Selected' -->
<!-- ### Plot -->
<!-- ggplot(df_count, aes(x = new.batch, y = biololy)) + -->
<!--   geom_count(aes(color = use), size = 7) + -->
<!--   geom_text(aes( -->
<!--     label = n, -->
<!--     hjust = 0.5, -->
<!--     vjust = 0.5 -->
<!--   )) + -->
<!--   xlab('Years-plates') + -->
<!--   ylab('Biological groups') + -->
<!--   theme_bw()+ -->
<!--   theme( -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text( -->
<!--       size = 12, -->
<!--       angle = 30, -->
<!--       hjust = 1 -->
<!--     ), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.position = 'none') -->
<!-- ``` -->
<!-- ### Library size of PRPS -->
<!-- Figure \@ref(fig:LsOfPrPs) shows the library sizes of the pseudo-samples of each pseudo-replicate sets across plates. As expected, they capture the large library size differences that we aim to remove from the data. -->
<!-- ```{r LsOfPrPs, message=FALSE, warning=FALSE, fig.align='right', fig.dim=c(6,7), fig.cap='Library sizes of pseudo-samples created in the TCGA READ RNA-Seq data'} -->
<!-- ### Make names out of plate ans years -->
<!-- ps.samples <- base::strsplit( -->
<!--   x = colnames(read.prps$ps.batch), -->
<!--   split = '_') -->
<!-- year <- lapply( -->
<!--   1:length(ps.samples), -->
<!--   function(x) { -->
<!--     index <- which(new.info$plate_RNAseq == ps.samples[[x]][3]) -->
<!--     unique(new.info$year_mda[index]) -->
<!--   }) -->
<!-- year.plate <- sapply( -->
<!--   1:length(ps.samples), -->
<!--   function(x) paste( -->
<!--     year[x], -->
<!--     ps.samples[[x]][3], -->
<!--     sep = '_' -->
<!--     )) -->
<!-- cms.msi <- sapply( -->
<!--   ps.samples, -->
<!--   function(x) paste( -->
<!--     x[1], -->
<!--     x[2], -->
<!--     sep = '_')) -->
<!-- ps <- data.frame( -->
<!--   bio = cms.msi, -->
<!--   year.plate = year.plate, -->
<!--   ls = log2(colSums(read.prps$ps.batch)) -->
<!--   ) -->
<!-- ### Plot -->
<!-- ggplot(ps, aes(x = year.plate, y = ls))  + -->
<!--   geom_point(size = 2) + -->
<!--   geom_line(aes(x = as.numeric(as.factor(year.plate)), y = ls)) + -->
<!--   xlab('Years_plates') + -->
<!--   ylab(expression(Log[2]~'library size')) + -->
<!--   facet_grid(bio ~.) + -->
<!--   theme( -->
<!--     axis.text.x = element_text(size = 8, angle = 30, hjust = 1), -->
<!--     axis.text.y = element_text(size = 8), -->
<!--     axis.title.x = element_text(size = 12), -->
<!--     axis.title.y = element_text(size = 12), -->
<!--     strip.text.y = element_text(size = 10)) -->
<!-- ``` -->
<!-- ### RUV-III-PRPS normalization -->
<!-- Here, we apply the RUV-III normalization with the PRPS and selected negative control genes. We refer to [(R.Molania, bioRxiv, 2021)](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1.article-metrics), and [(R.Molania, NAR, 2019)](https://academic.oup.com/nar/article/47/12/6073/5494770?login=true) for more details about the RUV-III method. -->
<!-- ```{r RuviiiNorm, message=FALSE, warning=FALSE} -->
<!-- ### prps -->
<!-- prps.batch <- read.prps$ps.batch -->
<!-- colnames(prps.batch) <- paste( -->
<!--   lapply( -->
<!--   colnames(prps.batch), -->
<!--   function(x){ -->
<!--     unlist(strsplit(x, '[_]'))[1] -->
<!--   }), -->
<!--   lapply( -->
<!--   colnames(prps.batch), -->
<!--   function(x){ -->
<!--     unlist(strsplit(x, '[_]'))[2] -->
<!--   }), -->
<!--   sep = '_') -->
<!-- ### ruv input data -->
<!-- ruv.data.input <- cbind( -->
<!--   SummarizedExperiment::assay( -->
<!--     read.cancer.se , -->
<!--     'HTseq_counts' -->
<!--     ), -->
<!--   prps.batch) -->
<!-- ### replicate matrix -->
<!-- ruv.rep.matrix <- ruv::replicate.matrix( -->
<!--   colnames(ruv.data.input)) -->
<!-- ruviii.norm <- RUV_III_PRPS( -->
<!--   Y = t(log2(ruv.data.input + 1)), -->
<!--   M = ruv.rep.matrix, -->
<!--   ctl = negative.control.genes, -->
<!--   k = 20, -->
<!--   eta = NULL, -->
<!--   return.info = TRUE -->
<!--   ) -->
<!-- ruviii.prps.norm <- t(ruviii.norm$newY[1:ncol(read.cancer.se), ]) -->
<!-- ``` -->
<!-- Note that, we remove all the PS from the data before any downstream analysis. -->
<!-- ## Consensus molecular subtypes of RUV-III normalized data -->
<!-- We apply the CMS classifier on the RUV-III normalized data. -->
<!-- ```{r CmsOnRuv, message=FALSE, error=FALSE, fig.align='center', fig.dim= c(5,5), fig.cap='Consensus molecular subtypes (CMS) indentification using the CMScaller R package in the RUV-III normalized data of the TCGA READ RNA-seq.'} -->
<!-- set.seed(2010221017) -->
<!-- row.names(ruviii.prps.norm) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$entrezgene_id_BioMart -->
<!-- cms.cluster.cancer.ruv <- CMScaller::CMScaller( -->
<!--   emat =  ruviii.prps.norm, -->
<!--   RNAseq = FALSE, -->
<!--   FDR = 0.05, -->
<!--   verbose = FALSE -->
<!--   ) -->
<!-- ### Addin new cms lables -->
<!-- index <- match( -->
<!--   row.names(cms.cluster.cancer.ruv), -->
<!--   colnames(read.cancer.se) -->
<!--   ) -->
<!-- read.cancer.se$cms.cancer.ruv[index] <- as.character(cms.cluster.cancer.ruv$prediction) -->
<!-- index <- is.na(read.cancer.se$cms.cancer.ruv) -->
<!-- read.cancer.se$cms.cancer.ruv[index] <- 'Not classified' -->
<!-- read.cancer.se$cms.cancer.ruv <- factor( -->
<!--   x = read.cancer.se$cms.cancer.ruv, -->
<!--   levels = c( -->
<!--     'CMS1', -->
<!--     'CMS2', -->
<!--     'CMS3', -->
<!--     'CMS4', -->
<!--     'Not classified')) -->
<!-- ``` -->
<!-- # Performance assessments for normalizations -->
<!-- We make use of both **global** and **gene-level** approaches to assess the performance of different normalization methods as removers of unwanted and preservers of biological variation in the data.\ -->
<!-- Our global approaches involve the use of principal component analysis (PCA) plots, linear regression, vector correlation analyses, silhouette coefficients, adjusted rand indices (ARI), and relative log expression (RLE) plots. Our PCA plots are each of the first three principal components (PC) against each other, coloured by known sources of unwanted variation, e.g. time, or known biology, e.g. cancer subtypes. Linear regression is used to quantify the relationship between the first few PC and continuous sources of unwanted variation such as (log) library size. The R2 calculated from the linear regression analyses indicates how strongly the PC capture unwanted variation in the data, and we do these calculations cumulatively, i.e. continuous source vs all of (PC1,…,PCk), for k = 1,…,5 or 10. Similar to linear regression, we used vector correlation analysis to assess the effect on the data of discrete sources of unwanted variation such as years or year intervals. Silhouette coefficients and adjusted Rand indices (ARI) were used to quantify how well experimental batches are mixed and known biology is separated. Finally, relative log expression (RLE) plots were used to assess the performance of different normalizations in terms of removing unwanted variation from the data.\ -->
<!-- The gene-level approach includes differential expression analyses between experimental batches, looking at p-value histograms and assessing the expression levels of negative control genes, positive control genes (genes whose behaviour we know), Spearman correlation and ANOVA between individual gene expression and sources of unwanted variation. These methods assess and quantify the effects of unwanted variation on individual gene expression levels in the RNA-seq datasets. We refer to the Methods section for more details about the assessment tools.\ -->
<!-- Here, we compare the performance of the RUV-III normalized data with the TCGA FPKM and FPKM.UQ datasets. We create a new SummarizedExperiment object that contains the RUV-III normalized data as well as the TCGA normalized datasets. -->
<!-- ```{r SeOnAllData, message=FALSE, warning=FALSE} -->
<!-- raw.count.data <- SummarizedExperiment::assay( -->
<!--   read.cancer.se,  -->
<!--   'HTseq_counts') -->
<!-- row.names(ruviii.prps.norm) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$gene_id.v -->
<!-- read.cancer.se <- SummarizedExperiment::SummarizedExperiment( -->
<!--   assays = list( -->
<!--     HTseq_counts = log2(SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_counts') + 1), -->
<!--     HTseq_FPKM = log2(SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_FPKM') + 1), -->
<!--     HTseq_FPKM.UQ = log2(SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_FPKM.UQ') + 1), -->
<!--     RUV_III = ruviii.prps.norm -->
<!--     ), -->
<!--   colData = S4Vectors::DataFrame( -->
<!--     SummarizedExperiment::colData(read.cancer.se)), -->
<!--   rowData = as.data.frame( -->
<!--     SummarizedExperiment::rowData(read.cancer.se)) -->
<!--   ) -->
<!-- read.sampleAnnot <- as.data.frame( -->
<!--   SummarizedExperiment::colData(read.cancer.se)) -->
<!-- normalizations <- names( -->
<!--   SummarizedExperiment::assays(read.cancer.se) -->
<!--   ) -->
<!-- normalizations.names <- c( -->
<!--   'Raw counts',  -->
<!--   'FPKM',  -->
<!--   'FPKM.UQ',  -->
<!--   'RUV-III') -->
<!-- ``` -->
<!-- ## Library size effects -->
<!-- Large library size variation between samples profiled in 2010 and the other samples are clearly visible in the PCA plots of the raw count data (figure \@ref(fig:LsEffectsAllPca)). Although the FPKM and FPKM.UQ normalizations reduce the variation caused by library size differences, both methods exhibited shortcomings, e.g. by not fully mixing samples from different times (figure \@ref(fig:LsEffectsAllPca)). PCA plots of the RUV-III normalized data illustrate that this normalization improved upon the FPKM and FPKM.UQ normalizations in removing the library size effects from the data. -->
<!-- ```{r LsEffectsAllPca, message=FALSE, warning=FALSE, fig.align='center', fig.dim=c(12,12), fig.cap='The scatter plots of first three principal components for raw counts, FPKM, FPKM.UQ and RUV-III normalized data coloured by key time points (2010 vs. 2011-2014).'} -->
<!-- pca.all <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     .pca( -->
<!--       data = as.matrix( -->
<!--         SummarizedExperiment::assay(read.cancer.se, x) -->
<!--         ), -->
<!--       is.log = TRUE) -->
<!--   }) -->
<!-- names(pca.all) <- normalizations -->
<!-- pp <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     pcs <- pca.all[[x]] -->
<!--     p <- .scatter.density.pc( -->
<!--       pcs = pcs$sing.val$u[,1:3], -->
<!--       pc.var = pcs$var, -->
<!--       group.name = 'Time (years)', -->
<!--       group = read.cancer.se$time.points, -->
<!--       color = major.times.colors, -->
<!--       strokeSize = .2, -->
<!--       pointSize = 3, -->
<!--       strokeColor = 'gray30', -->
<!--       alpha = .5) -->
<!--     p -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[1]], -->
<!--     pp[[2]], -->
<!--     pp[[3]], -->
<!--     pp[[4]], -->
<!--     ncol = 4)) -->
<!-- ``` -->
<!-- ### Association between PCs and library size -->
<!-- The first 5-10 PCs should have no or weak association with library size in an well- normalized dataset. The linear regression between the first ten PC, taken cumulatively, and library size clearly show that RUV-III outperforms other normalization in removing library size variation from the data (figure \@ref(fig:LsEffectsAllLreg)). -->
<!-- ```{r LsEffectsAllLreg, message=FALSE, warning=FALSE, fig.cap='Left-hand side: A plot showing the R-squared of linear regression between library size and up to the first 10 principal components (taken cumulatively) for different normalization methods. Right-hand side: percentage of PCs variation.'} -->
<!-- lreg.pcs.ls<- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     pcs <- pca.all[[x]]$sing.val$u -->
<!--     tcga.ls.rSquared <- sapply( -->
<!--       1:10, -->
<!--       function(y) { -->
<!--         lm.ls <- summary(lm( -->
<!--           read.sampleAnnot$libSize ~ pcs[, 1:y]) -->
<!--           )$r.squared -->
<!--     }) -->
<!--   }) -->
<!-- names(lreg.pcs.ls) <- normalizations -->
<!-- pcs.ls.lnreg <- as.data.frame(lreg.pcs.ls) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'r.sq') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III'))) -->
<!-- p1 <- ggplot(pcs.ls.lnreg, aes(x = pcs, y = r.sq, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = 1) + -->
<!--   geom_point(aes(color = datasets), size = 3) + -->
<!--   xlab('') +  -->
<!--   ylab (expression("R"^"2")) + -->
<!--   scale_color_manual( -->
<!--     values = c(dataSets.colors), -->
<!--     name = 'Datasets', -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III')) + -->
<!--   scale_x_continuous( -->
<!--     breaks = (1:10), -->
<!--     labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous( -->
<!--     breaks = scales::pretty_breaks(n = 5), -->
<!--     limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 35, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     legend.position = 'none') -->
<!-- # PCA variation -->
<!-- pcs.variation <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     pca.all[[x]]$variation[1:10] -->
<!--   }) -->
<!-- names(pcs.variation) <- normalizations -->
<!-- pcs.variation <- as.data.frame(pcs.variation) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'var') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III'))) %>% -->
<!--   data.frame() -->
<!-- ## plot -->
<!-- p2 <- ggplot(pcs.variation, aes(x = pcs, y = var, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = 1) + -->
<!--   geom_point(aes(color = datasets), size = 3) + -->
<!--   xlab('') + -->
<!--   ylab ('Percentage of explained variance') + -->
<!--   scale_color_manual( -->
<!--     values = dataSets.colors, -->
<!--     labels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III'), name = 'Datasets') + -->
<!--   scale_x_continuous( -->
<!--     breaks = (1:10), -->
<!--     labels = c(paste0('PC', 1:10) )) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 25, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1, -->
<!--   p2, -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- ### DE analysis between sample with low and high library size -->
<!-- Further, we evaluate the effects of library differences on the data using differential expression (DE) analyses between sample with low and high library size (2010 vs. 2011:2014). DE analyses were performed using the Wilcoxon signed-rank test with log2 transformed of the raw counts and normalized datasets.  In the absence of any library size effects, the histogram of the resulting unadjusted p-values should be uniformly distributed (Figure \@ref(fig:LsEffectsAllDe)). -->
<!-- ```{r LsEffectsAllDe, message=FALSE, warning=FALSE, fig.cap='P-value histograms obtained from differential expression analysis between samples with low and high library size.'} -->
<!-- de.ls.high.low <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     data <- SummarizedExperiment::assay(read.cancer.se, x) -->
<!--     de <- .wilcoxon.test( -->
<!--       expr.data = data, -->
<!--       is.log = TRUE, -->
<!--       variable = read.sampleAnnot$time.points, -->
<!--       n.cores = n.cores) -->
<!--     de -->
<!--   }) -->
<!-- names(de.ls.high.low) <- normalizations -->
<!-- pval.de.time.interval <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     de.ls.high.low[[x]]$pvalue -->
<!--   }) -->
<!-- names(pval.de.time.interval) <- normalizations -->
<!-- pval.de.time.interval <- pval.de.time.interval %>% -->
<!--   as.data.frame() %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'p.val') %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III'))) -->
<!-- ### Plot -->
<!-- ggplot(pval.de.time.interval, aes( p.val)) + -->
<!--   geom_histogram(binwidth = .1) + -->
<!--   scale_x_continuous(breaks = c(seq(0, 1, .5))) + -->
<!--   xlab('p_values') + ylab('Frequency') + -->
<!--   facet_wrap( ~ datasets, ncol = 4) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 16), -->
<!--     axis.title.y = element_text(size = 16), -->
<!--     plot.title = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 16)) -->
<!-- ``` -->
<!-- ### Association between gene expression and library size -->
<!-- Ideally, normalized gene expression values should have no significant association with library size in RNA-seq data. The relationship between individual normalized gene expression measurements and library size are assessed using Spearman correlation. Figure \@ref(fig:LsEffectsAllCorr) shows distribution of Spearman correlation coefficients between the gene expression levels and library size. These results show that there reasonable number of genes that have either positive or negative correlation with library size in the TCGA normalized datasets, whereases we did no see this in the RUV-III normalized data. -->
<!-- ```{r LsEffectsAllCorr, message=FALSE, warning=FALSE, fig.dim=c(6, 4), fig.cap='Spearman correlation coefficients between individual gene expression levels and library size.'} -->
<!-- corr.geneLs <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     .corr.gene.variable( -->
<!--       expr.data = as.matrix(SummarizedExperiment::assay(read.cancer.se, x)), -->
<!--       is.log = TRUE, -->
<!--       variable = read.sampleAnnot$libSize, -->
<!--       method = 'spearman', -->
<!--       n.cores = n.cores, -->
<!--       group = 'ls' -->
<!--       ) -->
<!--     }) -->
<!-- names(corr.geneLs) <- normalizations -->
<!-- gene.ls.corr.coeff <- lapply( -->
<!--   normalizations, -->
<!--   function(x) corr.geneLs[[x]]$ls_rho -->
<!--   ) -->
<!-- names(gene.ls.corr.coeff) <- normalizations -->
<!-- gene.ls.corr.coeff <- gene.ls.corr.coeff %>% -->
<!--   as.data.frame() %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'corr.coeff') %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III'))) -->
<!-- # plot -->
<!-- ggplot(gene.ls.corr.coeff, aes(x = datasets, y = corr.coeff, fill = datasets)) + -->
<!--   geom_violin() +  -->
<!--   ylab("Spearman correlation") + -->
<!--   xlab('') + -->
<!--   scale_fill_manual(values = dataSets.colors, guide = 'none') + -->
<!--     theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- ``` -->
<!-- ### Silhouette coefficient and ARI index analyses -->
<!-- We use Silhouette coefficient and ARI index to assess the performance of normalizations in terms of mixing samples from the two major time intervals. The results demonstrates that RUV-III removed the unwanted variation between samples from the major time intervals (Figure \@ref(fig:SilAriTime)). -->
<!-- ```{r SilAriTime, message=F, warning=F, fig.cap='Silhouette coefficient and ARI index analyses shows the performance of different methods in removing the unwanted variation between samples from the major time intervals.'} -->
<!-- silCoef.time <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     .silhouette.coeff( -->
<!--       pcs = pca.all[[x]]$sing.val$u, -->
<!--       variable = read.sampleAnnot$time.points, -->
<!--       nPCs = 3) -->
<!--     }) -->
<!-- names(silCoef.time) <- normalizations -->
<!-- pcs.time.silCoef <- as.data.frame(silCoef.time) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'silCoef') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) -->
<!--     ) -->
<!-- p1 <- ggplot(pcs.time.silCoef, aes(x = datasets, y = silCoef)) + -->
<!--   geom_col() + -->
<!--   ylab('Silhouette coefficient') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 14), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- # ARI -->
<!-- nPCs <- 3 -->
<!-- set.seed(2011110837) -->
<!-- ari.time <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     pcs <- pca.all[[x]]$sing.val$u[,1:nPCs] -->
<!--     BIC <- mclust::mclustBIC(data = pcs) -->
<!--     mod <- mclust::Mclust(data = pcs, x = BIC) -->
<!--     mclust::adjustedRandIndex( -->
<!--       mod$classification, -->
<!--       read.sampleAnnot$time.points -->
<!--       ) -->
<!--     }) -->
<!-- names(ari.time) <- normalizations -->
<!-- pcs.time.ari <- as.data.frame(ari.time) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'ari') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) -->
<!--     ) -->
<!-- # Plot -->
<!-- p2 <- ggplot(pcs.time.ari, aes(x = datasets, y = ari)) + -->
<!--   geom_col() + -->
<!--   ylab('ARI') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 14), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1,  -->
<!--   p2,  -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- ## Plate effects -->
<!-- To examine plate effects and separate this variation from the large library size variation in the data, we perform our evaluation within each major time interval. -->
<!-- ### Association between PCs and plates -->
<!-- Here, we apply PCA within each key time interval. Then, we use the Rozeboom squared vector correlation to quantify the strength of (linear) relationships between two sets of variables such as the first k PCs (i.e. 1≤k≤10) and dummy variables representing time, batches, plates, and biological variables [(R.Molania, bioRxiv, 2021)](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1.article-metrics). Here, we perform the vector correlation between the first ten PCs and plates. Figure \@ref(fig:PlateEffectsAllCaa) shows that the RUV-III performs better compared to the other normalization. -->
<!-- ```{r PlateEffectsAllCaa, message=FALSE, warning=FALSE, fig.cap='A plot showing the vector correlation coefficient between plates and the first 10 principal components within each time interval. The solid and dashed lines represent samples from 2010 and 2011:2014 respectively'} -->
<!-- ## PCA -->
<!-- pca.main.time.all <- lapply( -->
<!--   levels(read.sampleAnnot$time.points), -->
<!--   function(x){ -->
<!--     index <- read.sampleAnnot$time.points == x -->
<!--     pca.times  <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         pcs <- .pca( -->
<!--           data = as.matrix(SummarizedExperiment::assay( -->
<!--             read.cancer.se[ , index], y)), -->
<!--           is.log = TRUE) -->
<!--     }) -->
<!--     names(pca.times) <- normalizations -->
<!--     return(pca.times) -->
<!--   }) -->
<!-- names(pca.main.time.all) <- paste0( -->
<!--   'Time_', -->
<!--   levels(read.sampleAnnot$time.points)) -->
<!-- ## Vector correlation -->
<!-- cca.plates.time.interval <- lapply( -->
<!--   c(1:2), -->
<!--   function(x){ -->
<!--     pca.times <- pca.main.time.all[[x]] -->
<!--     index.time <- read.cancer.se$time.points == levels(read.cancer.se$time.points)[x] -->
<!--     plates.dummies <- fastDummies::dummy_cols(read.cancer.se$plate_RNAseq[index.time]) -->
<!--     plates.dummies <- plates.dummies[, c(2:ncol(plates.dummies))] -->
<!--     cca.plates <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         pcs <- pca.times[[y]]$sing.val$u -->
<!--         sapply( -->
<!--           1:10, -->
<!--           function(z) { -->
<!--             cca.plates <- stats::cancor( -->
<!--               x = pcs[, 1:z, drop = FALSE], -->
<!--               y = plates.dummies) -->
<!--             1 - prod(1 - cca.plates$cor ^ 2) -->
<!--           }) -->
<!--         }) -->
<!--     names(cca.plates) <- normalizations -->
<!--     cca.plates -->
<!--     }) -->
<!-- names(cca.plates.time.interval) <- levels( -->
<!--   read.cancer.se$time.points) -->
<!-- plate.time.interval.cca <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     as.data.frame(cca.plates.time.interval[[x]]) -->
<!--   }) %>% -->
<!--   do.call(rbind, .) %>% -->
<!--   as.data.frame() %>% -->
<!--     dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III) %>% -->
<!--   dplyr::mutate( -->
<!--     pcs = c(1:10, 1:10), -->
<!--     time = rep(c('2010', '2011_2014'), -->
<!--                each = 10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(pcs,time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'corr') %>% -->
<!--     dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III')) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- # Plot -->
<!-- ggplot() + geom_point( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time == '2010', ], -->
<!--     aes( -->
<!--       x = pcs, -->
<!--       y = corr , -->
<!--       group = datasets, -->
<!--       color = datasets -->
<!--     ), size = 3) + -->
<!--   geom_line(data = plate.time.interval.cca[plate.time.interval.cca$time == '2010',], -->
<!--             aes( -->
<!--               x = pcs, -->
<!--               y = corr , -->
<!--               group = datasets, -->
<!--               color = datasets -->
<!--             ), -->
<!--             linetype = "dashed", -->
<!--             size = 1) + -->
<!--   geom_point( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time != '2010',], -->
<!--     aes( -->
<!--     x = pcs, -->
<!--     y = corr , -->
<!--     group = datasets, -->
<!--     color = datasets -->
<!--   ), size = 3) + -->
<!--   scale_color_manual( -->
<!--     values = dataSets.colors, -->
<!--     name = 'Datasets', -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III')) + -->
<!--   geom_line( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time != '2010',], -->
<!--     aes( -->
<!--     x = pcs, -->
<!--     y = corr, -->
<!--     group = datasets, -->
<!--     color = datasets -->
<!--   ), size = 1) + -->
<!--    xlab('') + -->
<!--   ylab("Vector correlation") + -->
<!--   scale_color_manual( -->
<!--     values = dataSets.colors, -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III'), -->
<!--     name = 'Datasets') + -->
<!--   scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 14, angle = 25, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 12), -->
<!--     legend.title = element_text(size = 16)) -->
<!-- ``` -->
<!-- ### Association between gene expression and plates -->
<!-- Further, we use ANOVA to evaluate the plate effects on individual genes expression. -->
<!-- ```{r, warning=F, message=F, error=F, fig.cap='Boxplots of log2 F statistics obtained from ANOVA within each time points, for gene expression with plate as a factor'} -->
<!-- ftest.genePlates.time.interval <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     index.time <- read.cancer.se$time.points == x -->
<!--     ftest.plates.time <- lapply( -->
<!--       normalizations, -->
<!--       function(x){ -->
<!--         ftest.plates <- -->
<!--           .Ftest( -->
<!--             data = as.matrix(SummarizedExperiment::assay( -->
<!--               read.cancer.se[, index.time], -->
<!--               x) -->
<!--               ), -->
<!--             variable = read.cancer.se$plate_RNAseq[index.time], -->
<!--             is.log = TRUE, -->
<!--             n.cores = n.cores -->
<!--           ) -->
<!--       }) -->
<!--     names(ftest.plates.time) <- normalizations -->
<!--     ftest.plates.time -->
<!--   }) -->
<!-- names(ftest.genePlates.time.interval) <- levels( -->
<!--   read.cancer.se$time.points -->
<!--   ) -->
<!-- gene.plate.time.interval.ftest <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     sub <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         ftest.genePlates.time.interval[[x]][[y]]$FValue -->
<!--       }) -->
<!--    sub <-  do.call(cbind, sub) -->
<!--    colnames(sub) <- normalizations -->
<!--    sub -->
<!--   }) %>% -->
<!--   do.call(rbind, .) %>% -->
<!--     as.data.frame() %>% -->
<!--     dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III) %>% -->
<!--   dplyr::mutate( -->
<!--     time = rep(c('2010', '2011:2014'), -->
<!--                each = nrow(read.cancer.se)) ) %>% -->
<!--     tidyr::pivot_longer( -->
<!--     -c(time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'f.val') %>% -->
<!--     dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III')) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- ### Plot -->
<!-- ggplot( -->
<!--   data = gene.plate.time.interval.ftest, -->
<!--   aes(y = log2(f.val), x = datasets, fill = time)) + -->
<!--   geom_split_violin() + -->
<!--   scale_fill_manual(values = major.times.colors, name = 'Time (Years)') + -->
<!--   ylab(expression(Log[2]~'F statistics')) + -->
<!--   xlab('') + -->
<!--     theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 14, angle = 25, hjust = 1, vjust = 1), -->
<!--     axis.text.y = element_text(size = 14), -->
<!--     legend.text = element_text(size = 12), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- ``` -->
<!-- ## CMS clusters -->
<!-- Here, we evaluate the performance of different normalization methods in separating the CMS clusters. -->
<!-- ### PCA plots -->
<!-- PCA plots of the RUV-III normalized data show distinct clusters of the consensus molecular subtypes (CMS) for the READ RNA-seq samples, whereas these subtypes are not as clearly separated in the TCGA normalized datasets (Figure \@ref(fig:CmsAllPca)). -->
<!-- ```{r CmsAllPca, message=FALSE, warning=FALSE, error=FALSE, results=FALSE, fig.dim=c(12,12), fig.cap='PCA plots coloured by the CMS in the TCGA READ RNA-seq data normalized by different methods. From top to bottom: Raw counts, FPKM, FPKM.UQ and RUV-III'} -->
<!-- cms.cols <- c( -->
<!--   'cms.cancer.rawCounts', -->
<!--   'cms.cancer.fpkm', -->
<!--   'cms.cancer.fpkmUq', -->
<!--   'cms.cancer.ruv' -->
<!--   ) -->
<!-- pp <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     pcs <- pca.all[[x]] -->
<!--     p <- .scatter.density.pc( -->
<!--       pcs = pcs$sing.val$u[,1:3], -->
<!--       pc.var = pcs$var, -->
<!--       group.name = 'CMS subtypes', -->
<!--       group = read.sampleAnnot[ ,cms.cols[x]], -->
<!--       color = cms.colors, -->
<!--       strokeSize = .2, -->
<!--       pointSize = 3, -->
<!--       strokeColor = 'gray30', -->
<!--       alpha = .5) -->
<!--     p -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[1]], -->
<!--     pp[[2]], -->
<!--     pp[[3]], -->
<!--     pp[[4]], -->
<!--     ncol = 4)) -->
<!-- ``` -->
<!-- ### Association between PCs and CMS -->
<!-- Figure \@ref(fig:CmsAllCca) shows the vector correlation coefficient between CMS subtypes and the first 10 principal components for different normalization methods. Ideally, we should see high association between the PCs and the CMS. -->
<!-- ```{r CmsAllCca, message=FALSE, warning=FALSE, fig.cap='A plot showing the vector correlation coefficient between CMS subtypes and up to the first 10 principal components'} -->
<!-- cca.cms <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     cms.dummies <- fastDummies::dummy_cols(read.sampleAnnot[ , cms.cols[x]]) -->
<!--     cms.dummies <- cms.dummies[, c(2:ncol(cms.dummies))] -->
<!--     cms.caa.allPcs <- sapply( -->
<!--       1:10, -->
<!--       function(y) { -->
<!--         cms.caa <- stats::cancor( -->
<!--           x = pca.all[[x]]$sing.val$u[, 1:y, drop = FALSE], -->
<!--           y = cms.dummies) -->
<!--         1 - prod(1 - cms.caa$cor^2) -->
<!--     }) -->
<!--   }) -->
<!-- names(cca.cms) <- normalizations -->
<!-- pcs.cms.cca <- as.data.frame(cca.cms) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'vec.corr') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III')) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- # Plot -->
<!-- ggplot(pcs.cms.cca, aes(x = pcs, y = vec.corr, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = 1) + -->
<!--   geom_point(aes(color = datasets), size = 3) + -->
<!--   xlab('') + -->
<!--   ylab (expression("Vector correlation")) + -->
<!--   scale_color_manual( -->
<!--     values=c(dataSets.colors), -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III'), -->
<!--     name = 'Datastes') + -->
<!--   scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 14, angle = 25,vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 12), -->
<!--     legend.title = element_text(size = 16)) -->
<!-- ``` -->
<!-- ### Silhouette coefficient and ARI index analyses -->
<!-- Figure \@ref(fig:CmsAllSliho) shows that the RUV-III normalization results in a better separation of the CMS compared to the other normalizations.  -->
<!-- ```{r CmsAllSliho, message=FALSE, warning=FALSE, fig.dim=c(8,4), fig.cap='Silhouette coefficients and ARI index for mixing samples from two different key time intervals'} -->
<!-- silCoef.cms <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     .silhouette.coeff( -->
<!--       pcs = pca.all[[x]]$sing.val$u, -->
<!--       variable = read.sampleAnnot[, cms.cols[x]], -->
<!--       nPCs = 3) -->
<!--     }) -->
<!-- names(silCoef.cms) <- normalizations -->
<!-- pcs.cms.silCoef <- as.data.frame(silCoef.cms) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'silCoef') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) -->
<!--     ) -->
<!-- p1 <- ggplot(pcs.cms.silCoef, aes(x = datasets, y = silCoef)) + -->
<!--   geom_col() + -->
<!--   ylab('Silhouette coefficient') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     plot.title = element_text(size = 15), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- # ARI -->
<!-- nPCs <- 3 -->
<!-- set.seed(2011110837) -->
<!-- ari.cms <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     pcs <- pca.all[[x]]$sing.val$u[,1:nPCs] -->
<!--     BIC <- mclust::mclustBIC(data = pcs) -->
<!--     mod <- mclust::Mclust(data = pcs, x = BIC, G = 4) -->
<!--     mclust::adjustedRandIndex( -->
<!--       mod$classification, -->
<!--       read.sampleAnnot[, cms.cols[x]] -->
<!--       ) -->
<!--     }) -->
<!-- names(ari.cms) <- normalizations -->
<!-- pcs.cms.ari <- as.data.frame(ari.cms) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'ari') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) -->
<!--     ) -->
<!-- # Plot -->
<!-- p2 <- ggplot(pcs.cms.ari, aes(x = datasets, y = ari)) + -->
<!--   geom_col() + -->
<!--   ylab('ARI') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     plot.title = element_text(size = 15), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1,  -->
<!--   p2,  -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- ## Gene-level counts are not proportional to library size -->
<!-- The FPKM and FPKM.UQ normalizations rely on global scale factors computed based on library size or upper quartiles of samples in the raw count data to remove library size effects (Figure \@ref(fig:scaleFactors)). These methods assume that gene-level counts all are proportional to the global scale factors. -->
<!-- ```{r scaleFactors, message=F, fig.cap='Scale factors obtained by library size and upper-quartile of the raw counts.'} -->
<!-- uq <- apply( -->
<!--   raw.count.data,  -->
<!--   2,  -->
<!--   function(x) quantile(x, probs = c(.75))) -->
<!-- ls <- colSums(raw.count.data) -->
<!-- scale.factors <- data.frame( -->
<!--   LS = ls/mean(ls), -->
<!--   UQ = uq/mean(uq), -->
<!--   samples = c(1:166), -->
<!--   time = read.cancer.se$time.points) %>%  -->
<!--   tidyr::pivot_longer( -->
<!--     -c(samples, time),  -->
<!--     names_to = 'method',  -->
<!--     values_to = 'scaler') %>%  -->
<!--   data.frame(.) -->
<!-- ggplot(scale.factors, aes(x = samples, y = scaler, color = time)) + -->
<!--   geom_point(size = 2) + -->
<!--   facet_wrap(~method) + -->
<!--   scale_color_manual(values = major.times.colors, name = 'Time (Years)') + -->
<!--   ylab('Scale factor') + -->
<!--   xlab('Samples') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 24), -->
<!--     axis.title.y = element_text(size = 24), -->
<!--     axis.text.x = element_text(size = 14), -->
<!--     axis.text.y = element_text(size = 14), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 22)) -->
<!-- ``` -->
<!-- However, we show that in the READ raw count data different groups of genes exhibit different relationships to the global scale factors used in the FPKM and FPKM.UQ normalizations.\ -->
<!-- The first group consists of genes whose counts are proportional to the global scale factors. For these genes, the FPKM and FPKM.UQ normalizations are adequate to remove the association between library size variation and gene expression. The DEAD-Box Helicase 23 (DDX23) gene is an example from this group (Figure \@ref(fig:lsGeneLevels), first row).\ -->
<!-- The second group includes genes whose expression levels are greater than those expected using the global scaling factors, and so those factors are insufficient for adjusting their expression levels to be independent of library size. The La Ribonucleoprotein 7 (LARP7)_gene represents the behaviour of genes in this group (Figure \@ref(fig:lsGeneLevels), second row).\ -->
<!-- The third group contains genes such as AlkB Homolog 7 (ALKBH7), whose expression levels are not associated with library size in the raw count data. Then, the FPKM and FPKM.UQ normalizations introduce the library size variation to the expression levels of genes in this group (Figure \@ref(fig:lsGeneLevels), third row).\ -->
<!-- Finally, there are genes such as Transmembrane Protein 160 (TMEM160) whose expression levels relate to library size in a manner opposite to that motivating the use of global scaling factors. Applying scaling factors to such genes exacerbates rather than removes variation associated with library size (Figure \@ref(fig:lsGeneLevels), fourth row).\ -->
<!-- Note that we found the same issue in the TCGA RNA-seq datasets such as kidney chromophobe and uveal melanoma, where samples were profiled using a single plate. We refer to our vignette for on library size normalization [R.Molania, bioRxiv, 2021](https://www.biorxiv.org/content/10.1101/2021.11.01.466731v1). -->
<!-- ```{r lsGeneLevels, warning=F, message=F, fig.dim=c(12,12), fig.cap='Expression patterns of four genes (DDX23, LARP7, ALKBH7, TMEM160) whose counts have different relationships with the global scaling factors calculated from the READ raw count data.'} -->
<!-- row.names(read.cancer.se) <- SummarizedExperiment::rowData( -->
<!--   read.cancer.se -->
<!--   )$hgnc_symbol_BioMart -->
<!-- selected.genes <- c( -->
<!--   'DDX23', -->
<!--   'LARP7', -->
<!--   'ALKBH7', -->
<!--   'TMEM160') -->
<!-- pp <- lapply( -->
<!--   selected.genes, -->
<!--   function(i){ -->
<!--     if(i == 'DDX23') -->
<!--       { -->
<!--       df <- data.frame( -->
<!--     Raw.counts = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_counts'))[i, ], -->
<!--     FPKM = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM'))[i , ], -->
<!--     FPKM.UQ = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM.UQ'))[i , ], -->
<!--     RUV.III = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'RUV_III'))[i , ], -->
<!--     samples = c(1:166), -->
<!--     time = read.cancer.se$time.points -->
<!--   ) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(samples, time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'expr') %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, -->
<!--       grep('RUV.III', datasets), 'RUV-III')) %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, grep( -->
<!--         'Raw.counts', datasets), 'Raw counts')) %>% -->
<!--     dplyr::mutate(datasets = factor(datasets,  levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) ) %>% -->
<!--   as.data.frame(.) -->
<!--   p <- ggplot(df, aes(x = samples, y = expr, color = time)) + -->
<!--     geom_point(size = 3) + -->
<!--     scale_color_manual(values = major.times.colors) + -->
<!--     ylab(expression(Log[2] ~ 'gene expression'))  + -->
<!--     xlab('') + -->
<!--     facet_wrap( ~ datasets, scale = 'free', ncol = 4) + -->
<!--     ggtitle(i) + -->
<!--     theme( -->
<!--       panel.background = element_blank(), -->
<!--       axis.line = element_line(colour = 'black', size = 1), -->
<!--       plot.title = element_text(size = 14), -->
<!--       axis.title.x =  element_blank(), -->
<!--       axis.title.y =  element_blank(), -->
<!--       axis.text.x = element_text(size = 10), -->
<!--       axis.text.y = element_text(size = 12), -->
<!--       legend.text = element_text(size = 10), -->
<!--       legend.title = element_text(size = 14), -->
<!--       strip.text.x = element_text(size = 20), -->
<!--       legend.position = 'none' -->
<!--     ) +  -->
<!--     guides(color = guide_legend(title = "Time (years)")) -->
<!--   p -->
<!--     }else { -->
<!--       df <- data.frame( -->
<!--     Raw.counts = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_counts'))[i, ], -->
<!--     FPKM = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM'))[i , ], -->
<!--     FPKM.UQ = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM.UQ'))[i , ], -->
<!--     RUV.III = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'RUV_III'))[i , ], -->
<!--     samples = c(1:166), -->
<!--     time = read.cancer.se$time.points -->
<!--   ) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(samples, time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'expr') %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, -->
<!--       grep('RUV.III', datasets), 'RUV-III')) %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, grep( -->
<!--         'Raw.counts', datasets), 'Raw counts')) %>% -->
<!--     dplyr::mutate(datasets = factor(datasets,  levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III')) ) %>% -->
<!--   as.data.frame(.) -->
<!--   p <- ggplot(df, aes(x = samples, y = expr, color = time)) + -->
<!--     geom_point(size = 3) + -->
<!--     scale_color_manual(values = major.times.colors) + -->
<!--     ylab(expression(Log[2] ~ 'gene expression'))  + -->
<!--     xlab('Samples') + -->
<!--     facet_wrap( ~ datasets, scale = 'free', ncol = 4) + -->
<!--     ggtitle(i) + -->
<!--     theme( -->
<!--       panel.background = element_blank(), -->
<!--       axis.line = element_line(colour = 'black', size = 1), -->
<!--       plot.title = element_text(size = 14), -->
<!--       axis.title.x =  element_blank(), -->
<!--       axis.title.y =  element_blank(), -->
<!--       axis.text.x = element_text(size = 12), -->
<!--       axis.text.y = element_text(size = 12), -->
<!--       legend.text = element_text(size = 10), -->
<!--       legend.title = element_text(size = 14), -->
<!--       strip.text.x = element_text(size = 0)) +  -->
<!--     guides(color = guide_legend(title = "Time (years)")) -->
<!--   p -->
<!--     } -->
<!-- }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp, -->
<!--     ncol = 1)) -->
<!-- ``` -->
<!-- ## Gene co-expression analysis -->
<!-- The large sample library size differences in the data can compromise down-stream analyses such as gene co-expression. This variation can have two effects on gene co-expression analysis. It can lead to apparent correlations between genes that are most likely un-correlated. For example, the correlation between the TATA Element Modulatory Factor 1 (TMF1) and Bcl-2-associated transcription factor 1 (BCLAF1) genes are ρ=0.8 and ρ=0.7 in the TCGA FPKM and FPKM.UQ normalized data, respectively (Figure \@ref(fig:ArtiGeneCorr)). The role of the TMF1 gene has not been characterized in colon adenocarcinoma. Though, the BCLAF1 gene shows a pro-tumorigenic role in this cancer type [Xuexia Zhou et.al](https://www.nature.com/articles/ncomms5581). Then, one might suggest that the TMF1 gene expression may have a role in tumorigenesis in colon cancer due to its high correlation with the BCLAF1 gene expression. However, we see no such correlation in the RUV-III normalized data, which is consistent with the correlation obtain from an independent platform, namely the TCGA READ microarray data. -->
<!-- ```{r ArtiGeneCorr, warning=F, message=F, error=F, fig.dim=c(10,3), fig.cap='Gene co-expression analyses of TCGA READ RNA-seq data using different normalizations. Scatter plots of the gene expression levels of the MDH2 and EIF4H genes in the TCGA READ raw counts and differently normalized datasets. The red line shows overall association, and the green and yellow lines show associations between the gene expression within 2010 samples and within the rest of the samples, respectively.'} -->
<!-- genes.x <- c('MDH2', 'TMF1') -->
<!-- genes.y <- c('EIF4H', 'BCLAF1') -->
<!-- plot.title <- c( -->
<!--   'Raw counts', -->
<!--   'FPKM', -->
<!--   'FPKM.UQ', -->
<!--   'RUV-III' -->
<!--   ) -->
<!-- pp <- lapply( -->
<!--   c(1:2), -->
<!--   function(i){ -->
<!--   p.all <- lapply( -->
<!--     c(1:5), -->
<!--     function(x){ -->
<!--       if(x == 1){ -->
<!--         df <- SummarizedExperiment::assay(read.cancer.se, normalizations[x]) -->
<!--         df <- data.frame( -->
<!--           gene1 = df[genes.x[i] , ], -->
<!--           gene2 = df[genes.y[i] , ], -->
<!--           time = read.cancer.se$time.points -->
<!--         ) -->
<!--       p <- ggplot(df, aes(x = gene1, y = gene2, color = time)) + -->
<!--         geom_point(pch = 19, size = 1) + -->
<!--         scale_color_manual(values = major.times.colors) + -->
<!--         xlab(bquote(Log[2] ~ .(paste0('expression ', genes.x[i])))) + -->
<!--         ylab(bquote(Log[2] ~ .(paste0('expression ', genes.y[i])))) + -->
<!--         ggtitle(plot.title[x]) + -->
<!--         ggpubr::stat_cor( -->
<!--           aes(color = time, label = ..r.label..), -->
<!--           method = "spearman" , -->
<!--           label.x.npc = .87, -->
<!--           label.y.npc = .2, -->
<!--           hjust = 0, -->
<!--           r.accuracy = 0.1, -->
<!--           size = 3, -->
<!--           cor.coef.name = "rho" -->
<!--         ) + -->
<!--         geom_smooth( -->
<!--           method = 'lm', -->
<!--           formula = y ~ x, -->
<!--           col = 'red', -->
<!--           se = FALSE -->
<!--         ) + -->
<!--         geom_smooth( -->
<!--           aes(group = time), -->
<!--           method = 'lm', -->
<!--           formula = y ~ x, -->
<!--           se = FALSE) + -->
<!--         ggpubr::stat_cor( -->
<!--           aes(label = ..r.label..), -->
<!--           label.x.npc = .87, -->
<!--           label.y.npc = .33, -->
<!--           hjust = 0, -->
<!--           size = 3, -->
<!--           r.accuracy = 0.1, -->
<!--           col = 'red', -->
<!--           cor.coef.name = "rho" -->
<!--         )  + -->
<!--         gg.theme.2 + -->
<!--         guides(color = guide_legend(override.aes = list(size = 1))) -->
<!--       p -->
<!--       } else if(x > 1 & x< 5){ -->
<!--         df <- SummarizedExperiment::assay(read.cancer.se,  normalizations[x]) -->
<!--         df <- data.frame( -->
<!--           gene1 = df[genes.x[i] , ], -->
<!--           gene2 = df[genes.y[i] , ], -->
<!--           time = read.cancer.se$time.points -->
<!--         ) -->
<!--           p <- ggplot(df, aes(x = gene1, y = gene2, color = time)) + -->
<!--             geom_point(pch = 19, size = 1) + -->
<!--             ylab('') + -->
<!--             xlab('') + -->
<!--             scale_color_manual(values = major.times.colors) + -->
<!--             ggtitle(plot.title[x]) + -->
<!--             ggpubr::stat_cor( -->
<!--               aes(color = time, label = ..r.label..), -->
<!--               method = "spearman" , -->
<!--               label.x.npc = .87, -->
<!--               label.y.npc = .2, -->
<!--               hjust = 0, -->
<!--               size = 3, -->
<!--               r.accuracy = 0.1, -->
<!--               cor.coef.name = "rho" -->
<!--             ) + -->
<!--             geom_smooth( -->
<!--               method = 'lm', -->
<!--               formula = y ~ x, -->
<!--               col = 'red', -->
<!--               se = FALSE -->
<!--             ) + -->
<!--             geom_smooth(aes(group = time), -->
<!--                         method = 'lm', -->
<!--                         formula = y ~ x, -->
<!--                         se = FALSE) + -->
<!--             ggpubr::stat_cor( -->
<!--               aes(label = ..r.label..), -->
<!--               label.x.npc = .87, -->
<!--               label.y.npc = .33, -->
<!--               hjust = 0, -->
<!--               size = 3, -->
<!--               col = 'red', -->
<!--               r.accuracy = 0.1, -->
<!--               cor.coef.name = "rho" -->
<!--             )  + -->
<!--             gg.theme.2 + -->
<!--             guides(color = guide_legend(override.aes = list(size = 1))) -->
<!--           p -->
<!--           } else { -->
<!--           df <- data.frame( -->
<!--             gene1 = read.micro.array[genes.x[i] , ], -->
<!--             gene2 = read.micro.array[genes.y[i] , ] -->
<!--             ) -->
<!--           p <- ggplot(df, aes(x = gene1, y = gene2)) + -->
<!--             geom_point(pch = 19, size = 1) + -->
<!--             ggtitle('Microarray') + -->
<!--             ylab('') + -->
<!--             xlab('') + -->
<!--             geom_smooth( -->
<!--               method = 'lm', -->
<!--               formula = y ~ x, -->
<!--               col = 'red', -->
<!--               se = FALSE -->
<!--             ) + -->
<!--             ggpubr::stat_cor( -->
<!--               aes(label = ..r.label..), -->
<!--               label.x.npc = .87, -->
<!--               label.y.npc = .13, -->
<!--               hjust = 0, -->
<!--               size = 3, -->
<!--               col = 'red', -->
<!--               r.accuracy = 0.1, -->
<!--               cor.coef.name = "rho" -->
<!--             )  + -->
<!--             gg.theme.2 + -->
<!--             guides(color = guide_legend(override.aes = list(size = 1))) -->
<!--           p -->
<!--         } -->
<!--     }) -->
<!--   p.all -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[2]], ncol = 5)) -->
<!-- ``` -->
<!-- On the other hand, the unwanted variation can obscure correlations between gene-gene expression levels that likely to be truly correlated. For example, the overall correlation between the Malate Dehydrogenase 2 (MDH2) and Eukaryotic Translation Initiation Factor 4H (EIF4H) genes is ρ = -0.05, whereas they exhibit a high correlation within each key time interval in the TCGA normalized data (Figure \@ref(fig:ArtiGeneCorrExample)). The overall correlation of these genes was 0.7 in the RUV-III normalized data, consistent with what was seen in the TCGA READ microarray data. The MDH2 and EIF4H genes show important roles in cancer growth and metastasis, then, they are of clinically importance for cancer treatment [Zhan-Hong Chen et.al](https://pubmed.ncbi.nlm.nih.gov/31088567/), [Hyun Seung Ban et.al](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0162568)]. The high correlation between these two genes revealed by RUV-III may suggest that they are involved in a co-expression network, which would be a novel finding. -->
<!-- ```{r ArtiGeneCorrExample, warning=F, message=F, error=F, fig.dim=c(10,3), fig.cap='Gene co-expression analyses of TCGA READ RNA-seq data using different normalizations. Scatter plots of the gene expression levels of the BCLAF1 and EIF4H genes in the TCGA READ raw counts and differently normalized datasets. The red line shows overall association, and the green and yellow lines show associations between the gene expression within 2010 samples and within the rest of the samples, respectively.'} -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[1]], ncol = 5)) -->
<!-- ``` -->
<!-- ### Artifactual gene co-expression -->
<!-- We extend this analysis to all possible gene-gene correlations of the genes that have the highest correlation with library size in the FPKM.UQ normalized data (\@ref(fig:ArtiGeneCorrHeatMap). Strikingly, the results show numerous strong but likely spurious correlations between gene pairs in the FPKM.UQ normalized data, whereas using RUV-III significantly reduced these correlations. -->
<!-- ```{r ArtiGeneCorrHeatMap, warning=F, message=F, error=F, fig.cap=c('Gene co-expression analyses of TCGA READ RNA-seq data normalized by FPKM.UQ and RUV-III')} -->
<!-- # highly affected genes by major times -->
<!-- ftest.geneMajorTime <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     .Ftest( -->
<!--       data = as.matrix(SummarizedExperiment::assay(read.cancer.se, x)), -->
<!--       variable = read.cancer.se$time.points, -->
<!--       is.log = TRUE, -->
<!--       n.cores = n.cores -->
<!--       ) -->
<!--   }) -->
<!-- names(ftest.geneMajorTime) <- normalizations -->
<!-- selected.genes <- ftest.geneMajorTime$HTseq_FPKM.UQ$FValue %>% -->
<!--   tidyr::replace_na(., replace = 0) -->
<!-- index.genes <- which(log2(selected.genes + 1) >  6.43) -->
<!-- selected.genes <- ftest.geneMajorTime$HTseq_FPKM.UQ$Genes[index.genes] -->
<!-- ### correlation between genes and library size in TCGA FPKM.UQ -->
<!-- corr.gene.ls.fpkm.uq <- corr.geneLs$HTseq_FPKM.UQ$ls_rho[index.genes] -->
<!-- names(corr.gene.ls.fpkm.uq) <- selected.genes -->
<!-- cor.matrix.fpkm.uq <- cor(t( -->
<!--   SummarizedExperiment::assay( -->
<!--     read.cancer.se[selected.genes, ], -->
<!--     'HTseq_FPKM.UQ'))) -->
<!-- ### corrlation heat map -->
<!-- col_fun <- circlize::colorRamp2( -->
<!--   c(-0.8, 0, 0.8), -->
<!--   c("navy", "white", "yellow3") -->
<!--   ) -->
<!-- ha.annot.tcag = ComplexHeatmap::HeatmapAnnotation( -->
<!--   Corr_coef = corr.gene.ls.fpkm.uq, -->
<!--   col = list(Corr_coef = col_fun), -->
<!--   show_legend = TRUE, -->
<!--   show_annotation_name = FALSE) -->
<!-- corr.matrix.heatmap.fpkm.uq <- ComplexHeatmap::Heatmap( -->
<!--   cor.matrix.fpkm.uq, -->
<!--   col = rev(RColorBrewer::brewer.pal( -->
<!--     n = 11, -->
<!--     name = 'BrBG')), -->
<!--   show_row_names = FALSE, -->
<!--   name = 'Corr_coef', -->
<!--   show_column_names = FALSE, -->
<!--   top_annotation = ha.annot.tcag, -->
<!--   show_column_dend = FALSE, -->
<!--   show_row_dend = FALSE, -->
<!--   show_heatmap_legend = TRUE, -->
<!--   cluster_rows = TRUE, -->
<!--   cluster_columns = TRUE, -->
<!--   column_title = "TCGA FPKM.UQ" -->
<!--   ) -->
<!-- h <- ComplexHeatmap::draw(corr.matrix.heatmap.fpkm.uq) -->
<!-- order.rows <- ComplexHeatmap::row_order(h ) -->
<!-- ### RUV -->
<!-- corr.gene.ls.ruv <- corr.geneLs$RUV_III$ls_rho[index.genes] -->
<!-- names(corr.gene.ls.ruv) <- selected.genes -->
<!-- corr.gene.ls.ruv <- corr.gene.ls.ruv[order.rows] -->
<!-- cor.matrix.ruv <- cor(t( -->
<!--   SummarizedExperiment::assay( -->
<!--     read.cancer.se[selected.genes, ], -->
<!--     'RUV_III') -->
<!--   )) -->
<!-- cor.matrix.ruv <- cor.matrix.ruv[order.rows , order.rows] -->
<!-- col_fun = circlize::colorRamp2( -->
<!--   c(-0.8, 0, 0.8), -->
<!--   c("navy", "white", "yellow3") -->
<!--   ) -->
<!-- ha.annot.ruv = ComplexHeatmap::HeatmapAnnotation( -->
<!--   Corr_coef = corr.gene.ls.ruv, -->
<!--   col = list(Corr_coef = col_fun ), -->
<!--   show_legend = TRUE, -->
<!--   show_annotation_name = FALSE -->
<!--   ) -->
<!-- corr.matrix.heatmap.RUVIII <- ComplexHeatmap::Heatmap( -->
<!--   cor.matrix.ruv, -->
<!--   col = rev(RColorBrewer::brewer.pal( -->
<!--     n = 11, -->
<!--     name = 'BrBG')), -->
<!--    name = 'Corr_coef', -->
<!--   show_row_names = FALSE, -->
<!--   show_column_names = FALSE, -->
<!--   top_annotation = ha.annot.ruv, -->
<!--   show_column_dend = FALSE, -->
<!--   show_row_dend = FALSE, -->
<!--   show_heatmap_legend = TRUE, -->
<!--   cluster_rows = FALSE, -->
<!--   cluster_columns = FALSE, -->
<!--   column_title = "RUV-III" -->
<!--   ) -->
<!-- corr.matrix.heatmap.RUVIII -->
<!-- ``` -->
<!-- ## Association between gene expression and survival -->
<!-- Association between gene expression and survival outcomes of patients is another downstream analysis that is influenced by the library size variation in the TCGA READ RNA-seq data. For example, RUV-III, as opposed to the TCGA normalized data, revealed that the expression of the Ras-Related in Brain 18 (RAB18), and F-Box And Leucine Rich Repeat Protein 14 (FBXL14) genes are highly associated with overall survival outcome of patients in the data (figure \@ref(fig:GeneAndSurvival)). -->
<!-- ```{r GeneAndSurvival, warning=F, message=F, error=F, fig.dim=c(8,6), fig.cap='Association between gene expression and overall survival in the raw data and differently normalized datasets of the TCGA READ RNA-Seq data. Plots are for raw counts, FPKM, FPKM.UQ and RUV-III.'} -->
<!-- ## ggplot theme -->
<!-- selected.genes <- c( -->
<!--   'RAB18',  -->
<!--   'FBXL14', -->
<!--   'PTPN14', -->
<!--   'CSGALNACT2' -->
<!--   ) -->
<!-- data.sets.tcga <- c( -->
<!--   'read.rawCounts.cancer', -->
<!--   'read.fpkm.cancer', -->
<!--   'read.fpkmUq.cancer', -->
<!--   'read.ruv' -->
<!--   ) -->
<!-- for(i in selected.genes){ -->
<!--   sur.gene <- lapply( -->
<!--   normalizations[1:4], -->
<!--   function(x){ -->
<!--     p <- survival_plot( -->
<!--       data = as.data.frame(SummarizedExperiment::assay(read.cancer.se, x)), -->
<!--       stratify = 'expr', -->
<!--       annot = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--       scoreCol =  NULL, -->
<!--       gene = i, -->
<!--       covariate = NULL, -->
<!--       isCategoricalCov = FALSE, -->
<!--       timeCol = "OS.time_liu", -->
<!--       eventCol = "OS_liu", -->
<!--       nGroup = 2, -->
<!--       confInt = FALSE, -->
<!--       mainTitle1 = i, -->
<!--       ylabel = "Survival", -->
<!--       cols = c( -->
<!--         brewer.pal(9, "Set1")[c(2, 3, 4, 5, 7, 8)], -->
<!--         brewer.pal(8, "Dark2")[c(8, 1, 4, 6)]), -->
<!--       nColLegend = 1, -->
<!--       plotType = "autoplot") -->
<!--     return(p$plot) -->
<!--   }) -->
<!--   gridExtra::grid.arrange( -->
<!--     sur.gene[[1]] + ggplot.them, -->
<!--     sur.gene[[2]] + ggplot.them, -->
<!--     sur.gene[[3]] + ggplot.them, -->
<!--     sur.gene[[4]] + ggplot.them, -->
<!--     ncol = 2 -->
<!--   ) -->
<!-- } -->
<!-- ``` -->
<!-- The reason is clear from the expression patterns across time: dividing samples based on median expression mainly resulted in two groups with low and high library size, which was not biologically meaningful for the TCGA normalization (figure \@ref(fig:GeneAndSurvivalExamples)). The RAB18 gene expression plays pivotal roles in cell proliferation and metastasis, and high expression is associated with poor survival in different cancer types [Keng Zhong et.al](https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-14-703). The FBXL14 gene expression mediates the epithelial-mesenchymal transition in cancer, which indicates that the FBXL14 could function as an EMT inhibitor to suppress metastasis in human cancers [Yizuo Song et.al](https://stemcellres.biomedcentral.com/articles/10.1186/s13287-019-1222-0).\ -->
<!-- There are many more of these examples in the TCGA READ RNA-Seq data. -->
<!-- ```{r GeneAndSurvivalExamples, warning=F, message=F, error=F,fig.dim=c(8,6),fig.cap='Gene expression patterns of the RAB18 and FBXL14 genes in the raw data and differently normalized  datasets of the TCGA READ RNA-Seq data.'} -->
<!-- selected.genes <- c('RAB18', 'FBXL14') -->
<!-- p <- lapply(selected.genes, function(i){ -->
<!--   df <- data.frame( -->
<!--     Raw.counts = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_counts'))[i, ], -->
<!--     FPKM = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM'))[i , ], -->
<!--     FPKM.UQ = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'HTseq_FPKM.UQ'))[i , ], -->
<!--     RUV.III = unlist(SummarizedExperiment::assay( -->
<!--       read.cancer.se, 'RUV_III'))[i , ], -->
<!--     samples = c(1:166), -->
<!--     time = read.cancer.se$time.points -->
<!--   ) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(samples, time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'expr') %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, -->
<!--       grep('RUV.III', datasets), 'RUV-III')) %>% -->
<!--     dplyr::mutate(datasets = replace( -->
<!--       datasets, grep( -->
<!--         'Raw.counts', datasets), 'Raw counts')) %>% -->
<!--     dplyr::mutate(datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III'))) %>% -->
<!--     data.frame(.) -->
<!--   p <- ggplot(df, aes(x = samples, y = expr, color = time)) + -->
<!--     geom_point(size = 2) + -->
<!--     scale_color_manual(values = major.times.colors) + -->
<!--     ylab(expression(Log[2] ~ 'gene expression'))  + -->
<!--     xlab('Samples') + -->
<!--     ggtitle(i) + -->
<!--     facet_wrap( ~ datasets, scale = 'free', ncol = 4) + -->
<!--     theme( -->
<!--       panel.background = element_blank(), -->
<!--       axis.line = element_line(colour = 'black', size = 1), -->
<!--       plot.title = element_text(size = 14), -->
<!--       axis.title.x = element_text(size = 12), -->
<!--       axis.title.y = element_text(size = 12), -->
<!--       axis.text.x = element_text(size = 12), -->
<!--       axis.text.y = element_text(size = 12), -->
<!--       legend.text = element_text(size = 10), -->
<!--       legend.title = element_text(size = 12), -->
<!--       strip.text.x = element_text(size = 12), -->
<!--       legend.position = 'none' -->
<!--     ) + -->
<!--     guides(color = guide_legend(title = "Time (years)")) -->
<!--   p -->
<!-- }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange,  -->
<!--   p) -->
<!-- ``` -->
<!-- ## Tumour purity variation -->
<!-- Note that here we have not attempted to remove variation caused by tumor purity in the data. Consequently, the tumor purity estimates obtained from the RUV-III and FPKM.UQ normalized data were highly correlated (figure \@ref(fig:PurityEstimate)). This illustrates the ability of RUV-III to only remove the variation the user wishes to remove and no more, i.e. to retain other variation that is of biological origin. -->
<!-- ```{r PurityEstimate, message=F, warning=F, fig.dim=c(7,3), fig.cap='Scatter plot (left-hand side) shows tumour purity scores obtained from the TCGA FPKM.UQ and RUV-III normalized data. Scatter plot (right-hand side) shows tumour purity scores obtained from the TCGA FPKM.UQ (using the ESTIMATE method) and RUV-III normalized data (using the singscore method). Note, no attempt was made to remove tumour purity variation from the data by RUV-III normalization.'} -->
<!-- read.geneAnnot <- SummarizedExperiment::rowData(read.cancer.se) -->
<!-- purity.gene.sig <- read.geneAnnot$stromal == 'yes' | -->
<!--   read.geneAnnot$immnue == 'yes' -->
<!-- purity.gene.sig <- row.names(read.cancer.se)[purity.gene.sig] -->
<!-- row.names(ruviii.prps.norm) <- row.names(read.cancer.se) -->
<!-- rankDatas.ruv <- singscore::rankGenes(ruviii.prps.norm) -->
<!-- read.cancer.se$purity_RUV <- singscore::simpleScore( -->
<!--   rankData = rankDatas.ruv, -->
<!--   upSet = purity.gene.sig, -->
<!--   centerScore = F -->
<!--   )$TotalScore -->
<!-- p1 <- ggplot( -->
<!--   data = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--   aes(x = purity_HTseq_FPKM, y = 1-purity_RUV)) + -->
<!--   geom_point(size = 2) + -->
<!--   xlab('Tumour purity scores (FPKM.UQ)') + -->
<!--   ylab('Tumour purity scores (RUV-III)') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 12), -->
<!--     axis.title.y = element_text(size = 12), -->
<!--     axis.text.x = element_text(size = 10), -->
<!--     axis.text.y = element_text(size = 10) -->
<!--   ) -->
<!-- p2 <- ggplot( -->
<!--   data = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--   aes(x = ESTIMATE_Aran, y = 1-purity_RUV)) + -->
<!--   geom_point(size = 2) + -->
<!--   xlab('Tumour purity (ESTIMATE)') + -->
<!--   ylab('Tumour purity scores (RUV-III)') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 12), -->
<!--     axis.title.y = element_text(size = 12), -->
<!--     axis.text.x = element_text(size = 10), -->
<!--     axis.text.y = element_text(size = 10)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1, -->
<!--   p2, -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- # How robust RUV-III is to poorly chosen PRPS -->
<!-- Here, we assess the performance of RUV-III with poorly chosen PRPS on the TCGA READ data. To simulate poorly chosen PRPS, we randomly shuffle (20%, 40%, 60% and 80%) of the CMS subtypes that were originally used to create PRPS for RUV-III. The shuffling steps were repeated 10 times for each proportion and the results were averaged for normalization performance assessments. -->
<!-- ## PRPS with shuffling of CMS -->
<!-- ```{r PrpsShuffling, message=F, error=F, warning=F, results=F} -->
<!-- samples.to.use <- -->
<!--   read.cancer.se$cms.cancer.fpkmUq != 'Not classified' & -->
<!--   read.cancer.se$msi.status != 'Indeterminate' & -->
<!--   read.cancer.se$cms.use == 'yes' -->
<!-- sample.info <- droplevels( -->
<!--   as.data.frame( -->
<!--     SummarizedExperiment::colData(read.cancer.se[ , samples.to.use])) -->
<!--   ) -->
<!-- raw.counts <- raw.count.data[ , samples.to.use] -->
<!-- set.seed(703220347) -->
<!-- prps.random <- lapply( -->
<!--   c(seq(.2, .8, .2)),  -->
<!--   function(x){ -->
<!--     prps <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         a <- sample( -->
<!--           x = 1:nrow(sample.info),  -->
<!--           size = round(x = x*nrow(sample.info), digits = 0)) -->
<!--           b <- sample( -->
<!--             x = 1:nrow(sample.info), -->
<!--             size = round(x = x * nrow(sample.info), digits = 0)) -->
<!--           sample.info$cms.cancer.fpkmUq.2 <- sample.info$cms.cancer.fpkmUq -->
<!--           sample.info$msi.status.2 <- sample.info$msi.status -->
<!--           sample.info$cms.cancer.fpkmUq.2[a] <- sample.info$cms.cancer.fpkmUq[b] -->
<!--           sample.info$msi.status.2[a] <- sample.info$msi.status[b] -->
<!--           read.prps <- .CreatePseudoSamplesForLsPurityBatch( -->
<!--             expr.data = raw.counts, -->
<!--             sample.info = sample.info, -->
<!--             batch = 'PlateId_mda', -->
<!--             biology = c('cms.cancer.fpkmUq.2', 'msi.status.2'), -->
<!--             purity = FALSE, -->
<!--             include.ls = FALSE, -->
<!--             include.purity = FALSE, -->
<!--             minSamplesPerBatchPS = 2) -->
<!--       }) -->
<!--     names(prps) <- paste0('rep', 1:10) -->
<!--     prps -->
<!-- }) -->
<!-- names(prps.random) <- paste0( -->
<!--   "shuffle",  -->
<!--   seq(.2, .8, .2) -->
<!--   ) -->
<!-- prps.random.data <- lapply( -->
<!--   c(1:4),  -->
<!--   function(y){ -->
<!--     prps <- lapply( -->
<!--       c(1:10), -->
<!--       function(x){ -->
<!--     prps.batch <- prps.random[[y]][[x]]$ps.batch -->
<!--     dim(prps.batch) -->
<!--     colnames(prps.batch) <- paste( -->
<!--       lapply( -->
<!--         colnames(prps.batch), -->
<!--         function(x){ -->
<!--           unlist(strsplit(x, '[_]'))[1]}), -->
<!--       lapply( -->
<!--         colnames(prps.batch), -->
<!--         function(x){ -->
<!--           unlist(strsplit(x, '[_]'))[2]}), -->
<!--       sep = '_' ) -->
<!--     prps.batch -->
<!--   }) -->
<!--     names(prps) <- paste0('rep', 1:10) -->
<!--     prps -->
<!-- }) -->
<!-- names(prps.random.data) <- names(prps.random) -->
<!-- ``` -->
<!-- ## RUV-III normalization -->
<!-- We apply RUV-III with different sets of poorly chosen PRPS. -->
<!-- ```{r} -->
<!-- row.names(raw.count.data) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$hgnc_symbol_BioMart -->
<!-- ruviii.prps.random <- lapply( -->
<!--   c(1:4),  -->
<!--   function(x){ -->
<!--     ruv <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         prps <- prps.random.data[[x]][[y]] -->
<!--         ruv.data.input <- cbind( -->
<!--           raw.count.data, -->
<!--           prps) -->
<!--     ruv.rep.matrix.2 <- ruv::replicate.matrix(colnames(ruv.data.input)) -->
<!--     ruviii.norm.random <- RUV_III_PRPS( -->
<!--       Y = t(log2(ruv.data.input + 1)), -->
<!--       M = ruv.rep.matrix.2, -->
<!--       ctl = negative.control.genes, -->
<!--       k = 20, -->
<!--       eta = NULL, -->
<!--       return.info = TRUE) -->
<!--     ruviii.prps.random <- t(ruviii.norm.random$newY[1:ncol(read.cancer.se), ]) -->
<!--     ruviii.prps.random -->
<!--     }) -->
<!--     names(ruv) <- paste0('rep', 1:10) -->
<!--     ruv -->
<!--   }) -->
<!-- names(ruviii.prps.random) <- names(prps.random.data) -->
<!-- ``` -->
<!-- ## Performance assessments for normalizations -->
<!-- ### Library size effects -->
<!-- #### Association between PCs and library size -->
<!-- As we have mentioned above, the first 5-10 PCs should have weak association with library size in an well- normalized dataset. The linear regression between the first ten PC, taken cumulatively, and library size clearly shows that the RUV-III with poorly chosen PRPS outperforms the FPKM and FPKM.UQ normalizations in removing library size effects from the data (Figure \@ref(fig:LsPcsRandom)). -->
<!-- ```{r LsPcsRandom, message=FALSE, warning=FALSE, fig.cap='A plot showing the R-squared of linear regression between library size and up to the first 10 principal components (taken cumulatively) for different normalization methods.'} -->
<!-- ## PCA -->
<!-- pca.ruv.random <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     pcs <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         .pca( -->
<!--           data = ruviii.prps.random[[x]][[y]], -->
<!--           is.log = TRUE) -->
<!--         }) -->
<!--     names(pcs) <- paste0('rep', 1:10) -->
<!--     pcs -->
<!--   }) -->
<!-- names(pca.ruv.random) <- names(prps.random.data) -->
<!-- ### linear regression -->
<!-- lreg.pcs.ls.ruv.random <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     rSquared <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         pcs <- pca.ruv.random[[x]][[y]]$sing.val$u -->
<!--         ls.rSquared <- sapply( -->
<!--           1:10, -->
<!--           function(z) { -->
<!--             lm.ls <- summary(lm(read.sampleAnnot$libSize ~ pcs[, 1:z]))$r.squared -->
<!--           }) -->
<!--       }) -->
<!--     names(rSquared) <-  paste0('rep', 1:10) -->
<!--     rSquared -->
<!--   }) -->
<!-- names(lreg.pcs.ls.ruv.random) <- names(prps.random.data) -->
<!-- lreg.pcs.ls.ruv.random <- lapply( -->
<!--   c(1:4),  -->
<!--   function(x){ -->
<!--     rowMeans(do.call(cbind, lreg.pcs.ls.ruv.random[[x]])) -->
<!--   }) -->
<!-- names(lreg.pcs.ls.ruv.random) <- paste0( -->
<!--   'RUV_III_',  -->
<!--   seq(.2,.8, .2) -->
<!--   ) -->
<!-- lreg.pcs.ls.ruv.random <- do.call( -->
<!--   c,  -->
<!--   list( -->
<!--     lreg.pcs.ls,  -->
<!--     lreg.pcs.ls.ruv.random) -->
<!--   ) -->
<!-- pcs.ls.lnreg <- as.data.frame(lreg.pcs.ls.ruv.random) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-0.2' = RUV_III_0.2, -->
<!--     'RUV-III-0.4' = RUV_III_0.4, -->
<!--     'RUV-III-0.6' = RUV_III_0.6, -->
<!--     'RUV-III-0.8' = RUV_III_0.8, -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'r.sq') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-0.2', -->
<!--         'RUV-III-0.4', -->
<!--         'RUV-III-0.6', -->
<!--         'RUV-III-0.8' -->
<!--         )) -->
<!--     ) -->
<!-- ### plot -->
<!-- ggplot(pcs.ls.lnreg, aes(x = pcs, y = r.sq, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = .5) + -->
<!--   geom_point(aes(color = datasets), size = 2) + -->
<!--   xlab('PCs') + ylab (expression("R"^"2")) + -->
<!--   scale_color_manual( -->
<!--     values = c(dataSets.colors.2), -->
<!--     name = 'Datasets', -->
<!--     labels = names(dataSets.colors.2)) + -->
<!--   scale_x_continuous( -->
<!--     breaks = (1:10), -->
<!--     labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous( -->
<!--     breaks = scales::pretty_breaks(n = 5), -->
<!--     limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 35, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- ``` -->
<!-- ### Separation of CMS subtypes -->
<!-- Here, we evaluate the performance of different normalization methods in separating the CMS clusters. -->
<!-- #### Association between PCs and CMS -->
<!-- Figure \@ref(fig:CmsCcaPrpsRandom) shows the vector correlation coefficient between CMS subtypes and  the first 10 principal components for different normalization methods. The results show that the RUV-III with poorly chosen PRPS shows better performance compared to the TCGA normalizations in separating the CMS.  -->
<!-- ```{r CmsCcaPrpsRandom, message=FALSE, warning=FALSE, fig.cap='A plot showing the vector correlation coefficient between CMS subtypes and up to the first 10 principal components'} -->
<!-- cms.dummies <- fastDummies::dummy_cols(read.sampleAnnot$cms.cancer.ruv) -->
<!-- cms.dummies <- cms.dummies[, c(2:ncol(cms.dummies))] -->
<!-- cca.cms.ruv.random <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     cca.all <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         pcs <- pca.ruv.random[[x]][[y]] -->
<!--         cms.caa <- sapply( -->
<!--           1:10, -->
<!--           function(z) { -->
<!--             cms.caa <- -->
<!--               stats::cancor( -->
<!--                 x = pcs$sing.val$u[, 1:z, drop = FALSE], -->
<!--                 y = cms.dummies) -->
<!--             1 - prod(1 - cms.caa$cor ^ 2) -->
<!--             }) -->
<!--       }) -->
<!--     names(cca.all) <- paste0('rep', 1:10) -->
<!--     cca.all -->
<!--   }) -->
<!-- cca.cms.ruv.random <- lapply( -->
<!--   c(1:4),  -->
<!--   function(x){ -->
<!--     rowMeans(do.call(cbind, cca.cms.ruv.random[[x]])) -->
<!--   }) -->
<!-- names(cca.cms.ruv.random) <- paste0( -->
<!--   'RUV_III_',  -->
<!--   seq(.2,.8, .2) -->
<!--   ) -->
<!-- cca.cms.ruv.random  <- do.call( -->
<!--   c,  -->
<!--   list( -->
<!--     cca.cms,  -->
<!--     cca.cms.ruv.random)) -->
<!-- pcs.cms.cca <- as.data.frame(cca.cms.ruv.random) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-0.2' = RUV_III_0.2, -->
<!--     'RUV-III-0.4' = RUV_III_0.4, -->
<!--     'RUV-III-0.6' = RUV_III_0.6, -->
<!--     'RUV-III-0.8' = RUV_III_0.8 -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'vec.corr') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-0.2', -->
<!--         'RUV-III-0.4', -->
<!--         'RUV-III-0.6', -->
<!--         'RUV-III-0.8' -->
<!--         )) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- # Plot -->
<!-- ggplot(pcs.cms.cca, aes(x = pcs, y = vec.corr, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = .5) + -->
<!--   geom_point(aes(color = datasets), size = 2) + -->
<!--   xlab('PCs') + -->
<!--   ylab (expression("Vector correlation")) + -->
<!--   scale_color_manual( -->
<!--     values=c(dataSets.colors.2), -->
<!--     labels = names(dataSets.colors.2), -->
<!--     name = 'Datasets') + -->
<!--   scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     plot.title = element_text(size = 15), -->
<!--     axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 10) -->
<!--   ) -->
<!-- ``` -->
<!-- #### Silhouette coefficient  ann ARI index analyses -->
<!-- Silhouette coefficient  ann ARI index analyses (Figure \@ref(fig:SilAriCmsPrpsRandom)) confirm that the RUV-III with poorly chosen PRPS show better performance compared to the TCGA normalizations in separating the CMS.  -->
<!-- ```{r SilAriCmsPrpsRandom, message=FALSE, warning=FALSE, fig.dim=c(8,4), fig.cap='Silhouette coefficients and ARI index for separating the CMS'} -->
<!-- silCoef.cms.ruv.random <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     sil.coef <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         pcs <- pca.ruv.random[[x]][[y]] -->
<!--         .silhouette.coeff( -->
<!--           pcs = pcs$sing.val$u, -->
<!--           variable = read.sampleAnnot$cms.cancer.ruv, -->
<!--           nPCs = 3) -->
<!--       }) -->
<!--     names(sil.coef) <- paste0('rep', 1:10) -->
<!--     sil.coef -->
<!-- }) -->
<!-- silCoef.cms.ruv.random <- lapply( -->
<!--   c(1:4),  -->
<!--   function(x){ -->
<!--     mean(unlist(silCoef.cms.ruv.random[[x]])) -->
<!--   }) -->
<!-- names(silCoef.cms.ruv.random) <- paste0( -->
<!--   'RUV_III_',  -->
<!--   seq(.2,.8, .2) -->
<!--   ) -->
<!-- silCoef.cms.ruv.random <- do.call( -->
<!--   c,  -->
<!--   c(silCoef.cms,  -->
<!--     silCoef.cms.ruv.random)) -->
<!-- pcs.cms.silCoef <- as.data.frame(silCoef.cms.ruv.random) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(),  -->
<!--     names_to = 'silCoef.cms',  -->
<!--     values_to = 'silCoef') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III',  -->
<!--     'RUV-III-0.2', -->
<!--     'RUV-III-0.4', -->
<!--     'RUV-III-0.6', -->
<!--     'RUV-III-0.8' -->
<!--     )) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-0.2', -->
<!--       'RUV-III-0.4', -->
<!--       'RUV-III-0.6', -->
<!--       'RUV-III-0.8'))) -->
<!-- p1 <- ggplot(pcs.cms.silCoef, aes(x = datasets, y = silCoef)) + -->
<!--   geom_col() + -->
<!--   ylab('Silhouette coefficient') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- # ARI -->
<!-- set.seed(2011110837) -->
<!-- ari.cms.ruv.random <- lapply( -->
<!--   c(1:4), -->
<!--   function(x){ -->
<!--     ari <- lapply( -->
<!--       c(1:10),  -->
<!--       function(y){ -->
<!--         pcs <- pca.ruv.random[[x]][[y]]$sing.val$u[, 1:3] -->
<!--         BIC <- mclust::mclustBIC(data = pcs) -->
<!--         mod <- mclust::Mclust(data = pcs, x = BIC, G = 4) -->
<!--         mclust::adjustedRandIndex( -->
<!--           mod$classification, -->
<!--           read.sampleAnnot$cms.cancer.ruv) -->
<!--       }) -->
<!--     names(ari) <-  paste0('rep', 1:10) -->
<!--     ari -->
<!--   }) -->
<!-- names(ari.cms.ruv.random) <- names(prps.random.data) -->
<!-- ari.cms.ruv.random <- lapply( -->
<!--   c(1:4),  -->
<!--   function(x){ -->
<!--     mean(unlist(ari.cms.ruv.random[[x]])) -->
<!--   }) -->
<!-- names(ari.cms.ruv.random) <- paste0( -->
<!--   'RUV_III_',  -->
<!--   seq(.2,.8, .2)) -->
<!-- ari.cms.ruv.random <- do.call( -->
<!--   c,  -->
<!--   c(ari.cms,  -->
<!--     ari.cms.ruv.random)) -->
<!-- pcs.cms.ari <- as.data.frame(ari.cms.ruv.random) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(),  -->
<!--     names_to = 'silCoef.cms',  -->
<!--     values_to = 'ari') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III',  -->
<!--     'RUV-III-0.2', -->
<!--     'RUV-III-0.4', -->
<!--     'RUV-III-0.6', -->
<!--     'RUV-III-0.8' -->
<!--     )) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-0.2', -->
<!--       'RUV-III-0.4', -->
<!--       'RUV-III-0.6', -->
<!--       'RUV-III-0.8'))) -->
<!-- # Plot -->
<!-- p2 <- ggplot(pcs.cms.ari, aes(x = datasets, y = ari)) + -->
<!--   geom_col() + -->
<!--   ylab('ARI') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1,  -->
<!--   p2,  -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- ### Association between gene expression and survival -->
<!-- We assess the performance of RUV-III with poorly chosen PRPS in revealing the association between gene expression and survival .Figure \@ref(fig:GeneSurPrpsRandom) shows the p-values obtained from Kaplan Meier survival analyses for several genes. -->
<!-- ```{r GeneSurPrpsRandom, warning=F, message=F, error=F, fig.dim=c(8,6), fig.cap='Association between gene expression and overall survival in the raw data and differently normalized  datasets of the TCGA READ RNA-Seq data.'} -->
<!-- selected.genes <- c( -->
<!--   'RAB18',  -->
<!--   'FBXL14', -->
<!--   'PTPN14', -->
<!--   'CSGALNACT2') -->
<!-- data.sets.tcga <- c( -->
<!--   'read.rawCounts.cancer', -->
<!--   'read.fpkm.cancer', -->
<!--   'read.fpkmUq.cancer', -->
<!--   'read.ruv') -->
<!-- pval_1 <- lapply( -->
<!--   selected.genes,  -->
<!--   function(x){ -->
<!--     pval <- lapply( -->
<!--       normalizations,  -->
<!--       function(y){ -->
<!--         p <- survival_plot( -->
<!--           data = as.data.frame( -->
<!--             SummarizedExperiment::assay(read.cancer.se,y)), -->
<!--           stratify = 'expr', -->
<!--           annot = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--           scoreCol =  NULL, -->
<!--           gene = x, -->
<!--           covariate = NULL, -->
<!--           isCategoricalCov = FALSE, -->
<!--           timeCol = "OS.time_liu", -->
<!--           eventCol = "OS_liu", -->
<!--           nGroup = 2, -->
<!--           confInt = FALSE, -->
<!--           mainTitle1 = x, -->
<!--           ylabel = "Survival", -->
<!--           cols = brewer.pal(9, "Set1")[c(2, 3)], -->
<!--           nColLegend = 1, -->
<!--           plotType = "autoplot") -->
<!--         p$pval -->
<!--       }) -->
<!--     names(pval) <- normalizations.names -->
<!--     pval -->
<!--   }) -->
<!-- names(pval_1) <- selected.genes -->
<!-- ##  -->
<!-- pval_2 <- lapply( -->
<!--   selected.genes,  -->
<!--   function(g){ -->
<!--     lapply( -->
<!--       c(1:4), -->
<!--       function(x){ -->
<!--         all.pval <- lapply( -->
<!--           c(1:10),  -->
<!--           function(y){ -->
<!--             p <- survival_plot( -->
<!--               data = ruviii.prps.random[[x]][[y]], -->
<!--               stratify = 'expr', -->
<!--               annot = read.sampleAnnot, -->
<!--               scoreCol =  NULL, -->
<!--               gene = g, -->
<!--               covariate = NULL, -->
<!--               isCategoricalCov = FALSE, -->
<!--               timeCol = "OS.time_liu", -->
<!--               eventCol = "OS_liu", -->
<!--               nGroup = 2, -->
<!--               confInt = FALSE, -->
<!--               mainTitle1 = g, -->
<!--               ylabel = "Survival", -->
<!--               nColLegend = 1, -->
<!--               plotType = "autoplot" -->
<!--             ) -->
<!--         p$pval -->
<!--       }) -->
<!--     names(all.pval) <- paste0('rep', 1:10) -->
<!--     all.pval -->
<!--   }) -->
<!-- }) -->
<!-- names(pval_2) <- selected.genes -->
<!-- mean.pval_2 <- lapply( -->
<!--   selected.genes,  -->
<!--   function(g){ -->
<!--     pval <- unlist(lapply( -->
<!--       c(1:4),  -->
<!--       function(x) mean(unlist(pval_2[[g]][[x]])) -->
<!--       )) -->
<!--     names(pval) <- c(paste0( -->
<!--       'RUV-III.',  -->
<!--       seq(.2, .8, .2))) -->
<!--     pval -->
<!--   }) -->
<!-- names(mean.pval_2) <- selected.genes -->
<!-- pval.genes <- lapply( -->
<!--   selected.genes,  -->
<!--   function(x){ -->
<!--     c(unlist(pval_1[[x]]), mean.pval_2[[x]]) -->
<!--   }) -->
<!-- names(pval.genes) <- selected.genes -->
<!-- pval.genes <- data.frame(pval.genes) %>%  -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       row.names(.),  -->
<!--       levels = row.names(.))) %>%  -->
<!--   tidyr::pivot_longer( -->
<!--     -datasets,  -->
<!--     names_to = 'genes',  -->
<!--     values_to = 'pval') -->
<!-- ggplot(pval.genes, aes(x = datasets, y = pval)) + -->
<!--   geom_point() + -->
<!--   xlab('') + -->
<!--   ylab('P-value') + -->
<!--   facet_wrap(~genes) + -->
<!--   theme_bw() + -->
<!--   geom_hline(yintercept = 0.05, col = 'black', linetype = 'dotted') + -->
<!--   theme( -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text( -->
<!--       size = 12, -->
<!--       angle = 45, -->
<!--       vjust = 1, -->
<!--       hjust = 1 -->
<!--     ), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 16)) -->
<!-- ``` -->
<!-- ### Co-gene expression -->
<!-- Here, we explore the correlation between two pairs of genes including MDH2_EIF4H and TMF1_BCLAF1. The results show that the RUV-III with poorly chosen PRPS results in satisfactory normalization (Figure \@ref(fig:GeneSurPrpsRandom)). -->
<!-- ```{r GeneCoExpr, message=F, warning=F, error=F, fig.cap='Gene co-expression analyses of TCGA READ RNA-seq data using different normalizations. Plots show Spearman correlation analyses between two pairs of genes including MDH2_EIF4H and TMF1_BCLAF1 for different datasets.'}  -->
<!-- pair.genes <- list( -->
<!--   pair1 = c('MDH2', 'EIF4H'),  -->
<!--   pair2 = c('TMF1', 'BCLAF1')) -->
<!-- g.g.corr_1 <- lapply( -->
<!--   names(pair.genes),  -->
<!--   function(x){ -->
<!--     corr.coef <- lapply( -->
<!--       normalizations,  -->
<!--       function(y){ -->
<!--         data <- SummarizedExperiment::assay(read.cancer.se, y) -->
<!--         cor.test( -->
<!--           data[pair.genes[[x]][1] , ],  -->
<!--           data[pair.genes[[x]][2] , ],  -->
<!--           method = 'spearman')[[4]] -->
<!--         }) -->
<!--     names(corr.coef) <- normalizations.names -->
<!--     corr.coef -->
<!--     }) -->
<!-- names(g.g.corr_1) <- names(pair.genes) -->
<!-- g.g.corr_2 <- lapply( -->
<!--   names(pair.genes),  -->
<!--   function(g){ -->
<!--     corr <- lapply( -->
<!--       c(1:4),  -->
<!--       function(x){ -->
<!--         corr.all <- lapply( -->
<!--           c(1:10),  -->
<!--           function(y){ -->
<!--             data <- ruviii.prps.random[[x]][[y]] -->
<!--             cor.test( -->
<!--               data[pair.genes[[g]][1] ,], -->
<!--               data[pair.genes[[g]][2] ,], -->
<!--               method = 'spearman')[[4]]}) -->
<!--         mean(unlist(corr.all)) -->
<!--         }) -->
<!--     names(corr) <- paste0( -->
<!--       'RUV-III.',  -->
<!--       seq(.2, .8, .2)) -->
<!--     corr -->
<!--     }) -->
<!-- names(g.g.corr_2) <- names(pair.genes) -->
<!-- g.g.corr <- lapply( -->
<!--   names(pair.genes),  -->
<!--   function(x){ -->
<!--     c(unlist(g.g.corr_1[[x]]), unlist(g.g.corr_2[[x]])) -->
<!--   }) -->
<!-- names(g.g.corr) <- c('MDH2_EIF4H', 'TMF1_BCLAF1') -->
<!-- g.g.corr <- as.data.frame(g.g.corr) %>%  -->
<!--   dplyr::mutate(datasets = factor(row.names(.), row.names(.))) %>%  -->
<!--   dplyr::mutate(datasets = gsub('.rho', '', datasets)) %>%  -->
<!--   tidyr::pivot_longer(-datasets, names_to = 'genes', values_to = 'corr') -->
<!-- ggplot(g.g.corr, aes(x = datasets, y = corr)) + -->
<!--   geom_point() + -->
<!--   xlab('') + -->
<!--   ylab('Spearman correlation') + -->
<!--   facet_wrap(~genes) + -->
<!--   theme_bw() + -->
<!--   theme( -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text( -->
<!--       size = 12, -->
<!--       angle = 45, -->
<!--       vjust = 1, -->
<!--       hjust = 1 -->
<!--     ), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 10) -->
<!--   ) -->
<!-- ``` -->
<!-- # RUV-III with partially known biological labels for PRPS -->
<!-- We assess the performance of RUV-III with PRPS in situations where the biological labels are partially known (hereafter called the RUV-III-P). To simulate such situations, we only used the CMS4 (one of the CMS subtypes) to create RRPS for RUV-III normalization of the TCGA READ RNA-seq data. Note that, this subtype is not present across all the plates.  -->
<!-- ```{r PrPsGenerationP, message=FALSE, warning=FALSE, results=F} -->
<!-- samples.to.use <- -->
<!--   read.cancer.se$cms.cancer.fpkmUq == 'CMS4' & -->
<!--   read.cancer.se$msi.status != 'Indeterminate' & -->
<!--   read.cancer.se$cms.use == 'yes' -->
<!-- sample.info <- droplevels(read.sampleAnnot[samples.to.use , ]) -->
<!-- raw.counts <- raw.count.data[ , samples.to.use] -->
<!-- read.prps.par <- -->
<!--   .CreatePseudoSamplesForLsPurityBatch( -->
<!--     expr.data = raw.counts, -->
<!--     sample.info = sample.info, -->
<!--     batch = 'PlateId_mda', -->
<!--     biology = c('cms.cancer.fpkmUq', 'msi.status'), -->
<!--     purity = FALSE, -->
<!--     librarySize = 'libSize', -->
<!--     include.ls = TRUE, -->
<!--     include.purity = FALSE, -->
<!--     minSamplesPerBatchPS = 2, -->
<!--     minSamplesForLibrarySizePerBatch = 5, -->
<!--     minSamplesForLibrarySizePS = 2) -->
<!-- ``` -->
<!-- ## RUV-III-PRPS normalization -->
<!-- We apply RUV-III normalization on the TCGA READ RNA-seq data with the PRPS that were generated using the CMS4 only. -->
<!-- ```{r RuviiiNormP, message=FALSE, warning=FALSE} -->
<!-- ### prps -->
<!-- prps.batch <- read.prps.par$ps.batch -->
<!-- colnames(prps.batch) <- paste( -->
<!--   lapply( -->
<!--   colnames(prps.batch), -->
<!--   function(x){ -->
<!--     unlist(strsplit(x, '[_]'))[1] -->
<!--   }), -->
<!--   lapply( -->
<!--   colnames(prps.batch), -->
<!--   function(x){ -->
<!--     unlist(strsplit(x, '[_]'))[2] -->
<!--   }), -->
<!--   sep = '_' -->
<!--   ) -->
<!-- ### ruv input data -->
<!-- ruv.data.input <- cbind( -->
<!--   raw.count.data, -->
<!--   prps.batch, -->
<!--   read.prps.par$ps.ls -->
<!--   ) -->
<!-- ### replicate matrix -->
<!-- ruv.rep.matrix <- ruv::replicate.matrix( -->
<!--   colnames(ruv.data.input) -->
<!--   ) -->
<!-- ruviii.prps.par.norm <- RUV_III_PRPS( -->
<!--   Y = t(log2(ruv.data.input + 1)), -->
<!--   M = ruv.rep.matrix, -->
<!--   ctl = negative.control.genes, -->
<!--   k = 20, -->
<!--   eta = 1, -->
<!--   return.info = TRUE -->
<!--   ) -->
<!-- ruviii.prps.par <- t(ruviii.prps.par.norm$newY[1:ncol(read.cancer.se), ]) -->
<!-- ``` -->
<!-- ## Performance assessments for normalizations -->
<!-- Here, we compare the performance of the RUV-III-P normalized data with the original RUV-III normalized data and TCGA FPKM and FPKM.UQ datasets.\ -->
<!-- We create a new SummarizedExperiment object that contains the RUV-III normalized data as well as the TCGA normalized datasets. -->
<!-- ```{r SeOnAllDataP, message=FALSE, warning=FALSE} -->
<!-- row.names(read.cancer.se) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$gene_id.v -->
<!-- row.names(ruviii.prps.norm) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$gene_id.v -->
<!-- row.names(ruviii.prps.par) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$gene_id.v -->
<!-- read.cancer.se <- SummarizedExperiment::SummarizedExperiment( -->
<!--   assays = list( -->
<!--     HTseq_counts = SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_counts'), -->
<!--     HTseq_FPKM = SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_FPKM'), -->
<!--     HTseq_FPKM.UQ = SummarizedExperiment::assay( -->
<!--       read.cancer.se, -->
<!--       'HTseq_FPKM.UQ'), -->
<!--     RUV_III = ruviii.prps.norm, -->
<!--     RUV_III_p = ruviii.prps.par -->
<!--     ), -->
<!--   colData = S4Vectors::DataFrame( -->
<!--     SummarizedExperiment::colData(read.cancer.se)), -->
<!--   rowData = as.data.frame( -->
<!--     SummarizedExperiment::rowData(read.cancer.se))) -->
<!-- ``` -->
<!-- ### Library size effects -->
<!-- Figure \@ref(fig:LsPcsPlotPrpsPar) shows that  RUV-III-P removes the library size variation across samples. -->
<!-- ```{r LsPcsPlotPrpsPar, message=FALSE, warning=FALSE, fig.align='center', fig.dim=c(12,12), fig.cap='The scatter plots of first three principal components for raw counts, FPKM, FPKM.UQ and RUV-III normalized data coloured by key time points (2010 vs. 2011-2014).'} -->
<!-- normalizations <- names( -->
<!--   SummarizedExperiment::assays(read.cancer.se) -->
<!--   ) -->
<!-- pca.ruv.partial <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     .pca( -->
<!--       data = as.matrix( -->
<!--         SummarizedExperiment::assay(read.cancer.se, x) -->
<!--         ), -->
<!--       is.log = TRUE) -->
<!--   }) -->
<!-- names(pca.ruv.partial) <- normalizations -->
<!-- pp <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     pcs <- pca.ruv.partial[[x]] -->
<!--     p <- .scatter.density.pc( -->
<!--       pcs = pcs$sing.val$u[,1:3], -->
<!--       pc.var = pcs$var, -->
<!--       group.name = 'Time (years)', -->
<!--       group = read.cancer.se$time.points, -->
<!--       color = major.times.colors, -->
<!--       strokeSize = .2, -->
<!--       pointSize = 3, -->
<!--       strokeColor = 'gray30', -->
<!--       alpha = .5) -->
<!--     p -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[1]], -->
<!--     pp[[2]], -->
<!--     pp[[3]], -->
<!--     pp[[4]], -->
<!--     pp[[5]], -->
<!--     ncol = 4)) -->
<!-- ``` -->
<!-- #### Association between PCs and library size -->
<!-- The linear regression between the first ten PC, taken cumulatively, and library size clearly shows that the RUV-III-P outperforms TCGA normalizations in removing library size effects from the data (Figure \@ref(fig:LsPcsLregPrpsPar)). Note, the performance of the RUV-III-P is slightly lower than RUV-III. This might be due to presence of CMS4 in 8 out 14 plates.  -->
<!-- ```{r LsPcsLregPrpsPar, message=FALSE, warning=FALSE, fig.cap='A plot showing the R-squared (R2) of linear regression between library size and up to the first 10 principal components (taken cumulatively) for different normalization methods.'} -->
<!-- lreg.pcs.ls <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     pcs <- pca.ruv.partial[[x]]$sing.val$u -->
<!--     tcga.ls.rSquared <- sapply( -->
<!--       1:10, -->
<!--       function(y) { -->
<!--         lm.ls <- summary(lm( -->
<!--           read.sampleAnnot$libSize ~ pcs[, 1:y]) -->
<!--           )$r.squared -->
<!--     }) -->
<!--   }) -->
<!-- names(lreg.pcs.ls) <- normalizations -->
<!-- pcs.ls.lnreg <- as.data.frame(lreg.pcs.ls) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--      'RUV-III-P'= RUV_III_p -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'r.sq') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-P' -->
<!--         ))) -->
<!-- ggplot(pcs.ls.lnreg, aes(x = pcs, y = r.sq, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = .5) + -->
<!--   geom_point(aes(color = datasets), size = 2) + -->
<!--   xlab('PCs') + ylab (expression("R"^"2")) + -->
<!--   scale_color_manual( -->
<!--     values = c(dataSets.colors3), -->
<!--     name = 'Datasets', -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III', 'RUV-III-P')) + -->
<!--   scale_x_continuous( -->
<!--     breaks = (1:10), -->
<!--     labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous( -->
<!--     breaks = scales::pretty_breaks(n = 5), -->
<!--     limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 35, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- ``` -->
<!-- #### DE analysis between sample with low and high library size -->
<!-- Figure \@ref(fig:DeHighLowLsPrpsPar) shows that the RUV-III-P largely removes that association between genes and library size in the TCGA READ RNA-seq data. -->
<!-- ```{r DeHighLowLsPrpsPar, message=FALSE, warning=FALSE, fig.cap='P-value histograms obtained from differential expression analysis between samples with low and high library size.'} -->
<!-- de.ls.high.low <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     data <- SummarizedExperiment::assay(read.cancer.se, x) -->
<!--     de <- .wilcoxon.test( -->
<!--       expr.data = data, -->
<!--       is.log = TRUE, -->
<!--       variable = read.sampleAnnot$time.points, -->
<!--       n.cores = n.cores) -->
<!--     de -->
<!--   }) -->
<!-- names(de.ls.high.low) <- normalizations -->
<!-- pval.de.time.interval <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     de.ls.high.low[[x]]$pvalue -->
<!--   }) -->
<!-- names(pval.de.time.interval) <- normalizations -->
<!-- pval.de.time.interval <- pval.de.time.interval %>% -->
<!--   as.data.frame() %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-P'= RUV_III_p -->
<!--     ) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'p.val') %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-P' -->
<!--       ))) -->
<!-- ### Plot -->
<!-- ggplot(pval.de.time.interval, aes( p.val)) + -->
<!--   geom_histogram(binwidth = .1) + -->
<!--   scale_x_continuous(breaks = c(seq(0, 1, .5))) + -->
<!--   xlab('p_values') +  -->
<!--   ylab('Frequency') + -->
<!--   facet_wrap( ~ datasets, ncol = 5) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 16), -->
<!--     axis.title.y = element_text(size = 16), -->
<!--     plot.title = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 16)) -->
<!-- ``` -->
<!-- #### Association between gene expression and library size -->
<!-- Further, figure \@ref(fig:LsGeneCorrPrpsPar) shows that the RUV-III-P largely removes that association between genes and library size in the TCGA READ RNA-seq data. -->
<!-- ```{r LsGeneCorrPrpsPar, message=FALSE, warning=FALSE, fig.dim=c(12, 6), fig.cap='histograms of Spearman correlation coefficients between the gene expression levels and library size.'} -->
<!-- corr.geneLs <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     .corr.gene.variable( -->
<!--       expr.data = as.matrix(SummarizedExperiment::assay(read.cancer.se, x)), -->
<!--       is.log = TRUE, -->
<!--       variable = read.sampleAnnot$libSize, -->
<!--       method = 'spearman', -->
<!--       n.cores = n.cores, -->
<!--       group = 'ls') -->
<!--     }) -->
<!-- names(corr.geneLs) <- normalizations -->
<!-- gene.ls.corr.coeff <- lapply( -->
<!--   normalizations, -->
<!--   function(x) corr.geneLs[[x]]$ls_rho -->
<!--   ) -->
<!-- names(gene.ls.corr.coeff) <- normalizations -->
<!-- gene.ls.corr.coeff <- gene.ls.corr.coeff %>% -->
<!--   as.data.frame() %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-P' = RUV_III_p) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     everything(), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'corr.coeff') %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-P')) -->
<!--     ) -->
<!-- # plot -->
<!-- ggplot(gene.ls.corr.coeff, aes(corr.coeff, y = datasets, fill = datasets)) + -->
<!--   geom_violin(alpha = 0.7 ) + -->
<!--   xlab("Spearman correlation") + -->
<!--   ylab('') + -->
<!--   scale_fill_manual(values = dataSets.colors3, name = 'Datasets') + -->
<!--     theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 10)) -->
<!-- ``` -->
<!-- ### Plates effects -->
<!-- To examine plate effects and separate this variation from the large library size variation in the data, we perform our evaluation within each key time interval. -->
<!-- #### Association between PCs and plates -->
<!-- Figure \@ref(fig:PlateEffectsAllCaaPrpsPar) shows that the RUV-III-P performs better compared to the TCGA normalizations in removing plate effects. -->
<!-- ```{r PlateEffectsAllCaaPrpsPar, message=FALSE, warning=FALSE, fig.cap='A plot showing the vector correlation coefficient between plates and the first 10 principal components within each time interval'} -->
<!-- pca.main.time.all <- lapply( -->
<!--   levels(read.sampleAnnot$time.points), -->
<!--   function(x){ -->
<!--     index <- read.sampleAnnot$time.points == x -->
<!--     pca.times  <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         pcs <- .pca( -->
<!--           data = as.matrix(SummarizedExperiment::assay(read.cancer.se[ , index], y)), -->
<!--           is.log = TRUE) -->
<!--     }) -->
<!--     names(pca.times) <- normalizations -->
<!--     return(pca.times) -->
<!--   }) -->
<!-- names(pca.main.time.all) <- paste0( -->
<!--   'Time_', -->
<!--   levels(read.sampleAnnot$time.points)) -->
<!-- cca.plates.time.interval <- lapply( -->
<!--   c(1:2), -->
<!--   function(x){ -->
<!--     pca.times <- pca.main.time.all[[x]] -->
<!--     index.time <- read.cancer.se$time.points == levels(read.cancer.se$time.points)[x] -->
<!--     plates.dummies <- fastDummies::dummy_cols(read.cancer.se$plate_RNAseq[index.time]) -->
<!--     plates.dummies <- plates.dummies[, c(2:ncol(plates.dummies))] -->
<!--     cca.plates <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         pcs <- pca.times[[y]]$sing.val$u -->
<!--         sapply( -->
<!--           1:10, -->
<!--           function(z) { -->
<!--             cca.plates <- stats::cancor( -->
<!--               x = pcs[, 1:z, drop = FALSE], -->
<!--               y = plates.dummies) -->
<!--             1 - prod(1 - cca.plates$cor ^ 2) -->
<!--           }) -->
<!--         }) -->
<!--     names(cca.plates) <- normalizations -->
<!--     cca.plates -->
<!--     }) -->
<!-- names(cca.plates.time.interval) <- levels(read.cancer.se$time.points) -->
<!-- plate.time.interval.cca <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     as.data.frame(cca.plates.time.interval[[x]]) -->
<!--   }) %>% -->
<!--   do.call(rbind, .) %>% -->
<!--   as.data.frame() %>% -->
<!--     dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-P' = RUV_III_p) %>% -->
<!--   dplyr::mutate( -->
<!--     pcs = c(1:10, 1:10), -->
<!--     time = rep(c('2010', '2011_2014'), -->
<!--                each = 10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -c(pcs,time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'corr') %>% -->
<!--     dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-P' -->
<!--         )) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- # Plot -->
<!-- ggplot() + geom_point( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time == '2010', ], -->
<!--     aes( -->
<!--       x = pcs, -->
<!--       y = corr , -->
<!--       group = datasets, -->
<!--       color = datasets -->
<!--     ), size = 2) + -->
<!--   geom_line(data = plate.time.interval.cca[plate.time.interval.cca$time == '2010',], -->
<!--             aes( -->
<!--               x = pcs, -->
<!--               y = corr , -->
<!--               group = datasets, -->
<!--               color = datasets -->
<!--             ), -->
<!--             linetype = "dashed", -->
<!--             size = .5) + -->
<!--   geom_point( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time != '2010',], -->
<!--     aes( -->
<!--     x = pcs, -->
<!--     y = corr , -->
<!--     group = datasets, -->
<!--     color = datasets -->
<!--   ), size = 2) + -->
<!--   scale_color_manual( -->
<!--     values = dataSets.colors3, -->
<!--     name = 'Datasets', -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III', 'RUV-III-P')) + -->
<!--   geom_line( -->
<!--     data = plate.time.interval.cca[plate.time.interval.cca$time != '2010',], -->
<!--     aes( -->
<!--     x = pcs, -->
<!--     y = corr, -->
<!--     group = datasets, -->
<!--     color = datasets -->
<!--   ), size = .5) + -->
<!--    xlab('PCs') + -->
<!--   ylab("Vector correlation") + -->
<!--   scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 35, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- ``` -->
<!-- #### Association between gene expression and plates -->
<!-- Further, we use ANOVA to evaluate the plates effects on individual genes expression. Figure \@ref(fig:AnovaPlatesPrpsPar) shows that the RUV-III-P largely removes the plate effects on individual genes. -->
<!-- ```{r AnovaPlatesPrpsPar, warning=F, message=F, error=F, fig.cap='Boxplots of log2 F statistics obtained from ANOVA within each time points, for gene expression with plate as a factor'} -->
<!-- ftest.genePlates.time.interval <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     index.time <- read.cancer.se$time.points == x -->
<!--     ftest.plates.time <- lapply( -->
<!--       normalizations, -->
<!--       function(x){ -->
<!--         ftest.plates <- -->
<!--           .Ftest( -->
<!--             data = as.matrix(SummarizedExperiment::assay( -->
<!--               read.cancer.se[, index.time], -->
<!--               x) -->
<!--               ), -->
<!--             variable = read.cancer.se$plate_RNAseq[index.time], -->
<!--             is.log = TRUE, -->
<!--             n.cores = 5 -->
<!--           ) -->
<!--       }) -->
<!--     names(ftest.plates.time) <- normalizations -->
<!--     ftest.plates.time -->
<!--   }) -->
<!-- names(ftest.genePlates.time.interval) <- levels( -->
<!--   read.cancer.se$time.points -->
<!--   ) -->
<!-- gene.plate.time.interval.ftest <- lapply( -->
<!--   levels(read.cancer.se$time.points), -->
<!--   function(x){ -->
<!--     sub <- lapply( -->
<!--       normalizations, -->
<!--       function(y){ -->
<!--         ftest.genePlates.time.interval[[x]][[y]]$FValue -->
<!--       }) -->
<!--    sub <- do.call(cbind, sub) -->
<!--    colnames(sub) <- normalizations -->
<!--    sub}) %>% -->
<!--   do.call(rbind, .) %>% -->
<!--     as.data.frame() %>% -->
<!--     dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-P' = RUV_III_p) %>% -->
<!--   dplyr::mutate( -->
<!--     time = rep(c('2010', '2011:2014'), -->
<!--                each = nrow(read.cancer.se)) ) %>% -->
<!--     tidyr::pivot_longer( -->
<!--     -c(time), -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'f.val') %>% -->
<!--     dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-P')) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- ### Plot -->
<!-- ggplot( -->
<!--   data = gene.plate.time.interval.ftest, -->
<!--   aes(x = datasets, y = log2(f.val), fill = time)) + -->
<!--   geom_boxplot( -->
<!--     position = position_dodge(width = .7), -->
<!--     width = 1, -->
<!--     alpha = 0.8, -->
<!--     lwd = .7 -->
<!--   ) + -->
<!--   scale_fill_manual(values = major.times.colors, name = ('Time (years)')) + -->
<!--   ylab(expression(Log[2]~'F statistics')) + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 16), -->
<!--     axis.title.y = element_text(size = 16), -->
<!--     plot.title = element_text(size = 15), -->
<!--     axis.text.x = element_text(size = 16, angle = 25, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 14), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14), -->
<!--     strip.text.x = element_text(size = 10)) -->
<!-- ``` -->
<!-- ### CMS clusters -->
<!-- Here, we evaluate the performance of different normalization methods in separating the CMS clusters. -->
<!-- #### PCA plots -->
<!-- Figure \@ref(fig:CmsAllPcaPrpsPar) shows the RUV-III-P normalization led to distinct clusters of the consensus molecular subtypes (CMS) for the READ RNA-seq samples, whereas these subtypes are not as clearly separated in the TCGA normalized datasets. -->
<!-- ```{r CmsAllPcaPrpsPar, message=FALSE, warning=FALSE, error=FALSE, results=FALSE, fig.dim=c(12,12), fig.cap='PCA plots of different datasets (from top to bottom: raw counts, FPKM, FPKM.UQ, RUV-III and RUV-III-P normalized data coloured by the CMS.'} -->
<!-- cms.cols <- c( -->
<!--   'cms.cancer.rawCounts', -->
<!--   'cms.cancer.fpkm', -->
<!--   'cms.cancer.fpkmUq', -->
<!--   'cms.cancer.ruv', -->
<!--   'cms.cancer.ruv' -->
<!--   ) -->
<!-- pp <- lapply( -->
<!--   c(1:5), -->
<!--   function(x){ -->
<!--     pcs <- pca.ruv.partial[[x]] -->
<!--     p <- .scatter.density.pc( -->
<!--       pcs = pcs$sing.val$u[,1:3], -->
<!--       pc.var = pcs$var, -->
<!--       group.name = 'CMS', -->
<!--       group = read.sampleAnnot[ ,cms.cols[x]], -->
<!--       color = cms.colors, -->
<!--       strokeSize = .2, -->
<!--       pointSize = 3, -->
<!--       strokeColor = 'gray30', -->
<!--       alpha = .5) -->
<!--     p -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   c(pp[[1]], -->
<!--     pp[[2]], -->
<!--     pp[[3]], -->
<!--     pp[[4]], -->
<!--     pp[[5]], -->
<!--     ncol = 4)) -->
<!-- ``` -->
<!-- #### Association between PCs and CMS -->
<!-- Figure (\@ref(fig:CmsAllCcaPrpsPar)) shows the vector correlation coefficient between CMS subtypes and  the first 10 principal components for different normalization methods. Ideally, we should see high association between the PCs and the CMS. The RUV-III-P normalization outperforms the TCGA normalizations in separating the CMS subtypes. -->
<!-- ```{r CmsAllCcaPrpsPar, message=FALSE, warning=FALSE, fig.cap='A plot showing the vector correlation coefficient between CMS subtypes and up to the first 10 principal components'} -->
<!-- cms.cols <- c( -->
<!--   'cms.cancer.rawCounts', -->
<!--   'cms.cancer.fpkm', -->
<!--   'cms.cancer.fpkmUq', -->
<!--   'cms.cancer.ruv', -->
<!--   'cms.cancer.ruv' -->
<!--   ) -->
<!-- nPCs <- 10 -->
<!-- cca.cms <- lapply( -->
<!--   c(1:5), -->
<!--   function(x){ -->
<!--     cms.dummies <- fastDummies::dummy_cols(read.sampleAnnot[ , cms.cols[x]]) -->
<!--     cms.dummies <- cms.dummies[, c(2:ncol(cms.dummies))] -->
<!--     cms.caa.allPcs <- sapply( -->
<!--       1:nPCs, -->
<!--       function(y) { -->
<!--         cms.caa <- stats::cancor( -->
<!--           x = pca.ruv.partial[[x]]$sing.val$u[, 1:y, drop = FALSE], -->
<!--           y = cms.dummies) -->
<!--         1 - prod(1 - cms.caa$cor^2) -->
<!--     }) -->
<!--   }) -->
<!-- names(cca.cms) <- normalizations -->
<!-- pcs.cms.cca <- as.data.frame(cca.cms) %>% -->
<!--   dplyr::rename( -->
<!--     'Raw counts' = HTseq_counts, -->
<!--     FPKM = HTseq_FPKM, -->
<!--     FPKM.UQ = HTseq_FPKM.UQ, -->
<!--     'RUV-III' = RUV_III, -->
<!--     'RUV-III-P'= RUV_III_p -->
<!--   ) %>% -->
<!--   dplyr::mutate(pcs = c(1:10)) %>% -->
<!--   tidyr::pivot_longer( -->
<!--     -pcs, -->
<!--     names_to = 'datasets', -->
<!--     values_to = 'vec.corr') %>% -->
<!--   dplyr::mutate( -->
<!--     datasets = factor( -->
<!--       datasets, -->
<!--       levels = c( -->
<!--         'Raw counts', -->
<!--         'FPKM', -->
<!--         'FPKM.UQ', -->
<!--         'RUV-III', -->
<!--         'RUV-III-P')) -->
<!--     ) %>% -->
<!--   data.frame(.) -->
<!-- # Plot -->
<!-- ggplot(pcs.cms.cca, aes(x = pcs, y = vec.corr, group = datasets)) + -->
<!--   geom_line(aes(color = datasets), size = .5) + -->
<!--   geom_point(aes(color = datasets), size = 2) + -->
<!--   xlab('PCs') + -->
<!--   ylab (expression("Vector correlation")) + -->
<!--   scale_color_manual( -->
<!--     values=c(dataSets.colors3), -->
<!--     labels = c('Raw counts', 'FPKM','FPKM.UQ', 'RUV-III', 'RUV-III-P'), -->
<!--     name = 'Datasets') + -->
<!--   scale_x_continuous(breaks = (1:10), labels = c('PC1', paste0('PC1:', 2:10)) ) + -->
<!--   scale_y_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(0,1)) + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = 1), -->
<!--     axis.title.x = element_text(size = 18), -->
<!--     axis.title.y = element_text(size = 18), -->
<!--     axis.text.x = element_text(size = 12, angle = 35, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12), -->
<!--     legend.text = element_text(size = 10), -->
<!--     legend.title = element_text(size = 14)) -->
<!-- ``` -->
<!-- #### Silhouette coefficient  ann ARI index analyses -->
<!-- Silhouette coefficients analysis and Adjusted Rand Index show that the RUV-III-P performs better than the TCGA normalization in separating the CMS subtypes. -->
<!-- ```{r CmsSilAriPrpsPar, message=FALSE, warning=FALSE, fig.dim=c(8,4), fig.cap='Silhouette coefficients and ARI index for mixing samples from two different key time intervals'} -->
<!-- silCoef.cms <- lapply( -->
<!--   c(1:5), -->
<!--   function(x){ -->
<!--     .silhouette.coeff( -->
<!--       pcs = pca.ruv.partial[[x]]$sing.val$u, -->
<!--       variable = read.sampleAnnot[, cms.cols[x]], -->
<!--       nPCs = 3) -->
<!--     }) -->
<!-- names(silCoef.cms) <- normalizations -->
<!-- pcs.cms.silCoef <- as.data.frame(silCoef.cms) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'silCoef') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III', -->
<!--     'RUV-III-P')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-P')) -->
<!--     ) -->
<!-- p1 <- ggplot(pcs.cms.silCoef, aes(x = datasets, y = silCoef)) + -->
<!--   geom_col() + -->
<!--   ylab('Silhouette coefficient') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 12, angle = 25, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- # ARI -->
<!-- nPCs <- 3 -->
<!-- set.seed(2011110837) -->
<!-- ari.cms <- lapply( -->
<!--   c(1:5), -->
<!--   function(x){ -->
<!--     pcs <- pca.ruv.partial[[x]]$sing.val$u[,1:nPCs] -->
<!--     BIC <- mclust::mclustBIC(data = pcs) -->
<!--     mod <- mclust::Mclust(data = pcs, x = BIC, G = 4) -->
<!--     mclust::adjustedRandIndex( -->
<!--       mod$classification, -->
<!--       read.sampleAnnot[, cms.cols[x]] -->
<!--       ) -->
<!--     }) -->
<!-- names(ari.cms) <- normalizations -->
<!-- pcs.cms.ari <- as.data.frame(ari.cms) %>% -->
<!--   tidyr::pivot_longer(everything(), names_to = 'silCoef.cms', values_to = 'ari') %>% -->
<!--   dplyr::mutate(datasets = c( -->
<!--     'Raw counts', -->
<!--     'FPKM', -->
<!--     'FPKM.UQ', -->
<!--     'RUV-III', -->
<!--     'RUV-III-P')) %>% -->
<!--   dplyr::mutate(datasets = factor( -->
<!--     datasets, -->
<!--     levels = c( -->
<!--       'Raw counts', -->
<!--       'FPKM', -->
<!--       'FPKM.UQ', -->
<!--       'RUV-III', -->
<!--       'RUV-III-P')) -->
<!--     ) -->
<!-- # Plot -->
<!-- p2 <- ggplot(pcs.cms.ari, aes(x = datasets, y = ari)) + -->
<!--   geom_col() + -->
<!--   ylab('ARI') + -->
<!--   xlab('') + -->
<!--   theme( -->
<!--     panel.background = element_blank(), -->
<!--     axis.line = element_line(colour = 'black', size = .85), -->
<!--     axis.title.x = element_text(size = 14), -->
<!--     axis.title.y = element_text(size = 14), -->
<!--     axis.text.x = element_text(size = 12, angle = 25, vjust = 1, hjust = 1), -->
<!--     axis.text.y = element_text(size = 12)) -->
<!-- gridExtra::grid.arrange( -->
<!--   p1,  -->
<!--   p2,  -->
<!--   ncol = 2) -->
<!-- ``` -->
<!-- ## Association between gene expression and survival -->
<!-- Figure \@ref(fig:GeneSurPrpsPar) shows that the RUV-III-P reveals the association between gene expression and overall survival in the TCGA READ RNA-seq data. However, we did not see such association for FBXK14 in the RUV-III-P normalized data. This might be due to presences of our PRPS in 8 out of 14 plates. -->
<!-- ```{r GeneSurPrpsPar, warning=F, message=F, error=F, fig.dim=c(8,6), fig.cap='Association between gene expression and overall survival in the raw data and differently normalized  datasets of the TCGA READ RNA-Seq data. from left to right: raw counts, FPKM, FPKM.UQ, RUV-III and RUV-III-P'} -->
<!-- names(read.cancer.se) <- as.data.frame( -->
<!--   SummarizedExperiment::rowData(read.cancer.se) -->
<!--   )$hgnc_symbol_BioMart -->
<!-- selected.genes <- c( -->
<!--   'RAB18',  -->
<!--   'FBXL14', -->
<!--   'PTPN14', -->
<!--   'CSGALNACT2' -->
<!--   ) -->
<!-- for(i in selected.genes){ -->
<!--   sur.gene <- lapply( -->
<!--   normalizations, -->
<!--   function(x){ -->
<!--     p <- survival_plot( -->
<!--       data = as.data.frame(SummarizedExperiment::assay(read.cancer.se, x)), -->
<!--       stratify = 'expr', -->
<!--       annot = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--       scoreCol =  NULL, -->
<!--       gene = i, -->
<!--       covariate = NULL, -->
<!--       isCategoricalCov = FALSE, -->
<!--       timeCol = "OS.time_liu", -->
<!--       eventCol = "OS_liu", -->
<!--       nGroup = 2, -->
<!--       confInt = FALSE, -->
<!--       mainTitle1 = i, -->
<!--       ylabel = "Survival", -->
<!--       cols = c( -->
<!--         brewer.pal(9, "Set1")[c(2, 3, 4, 5, 7, 8)], -->
<!--         brewer.pal(8, "Dark2")[c(8, 1, 4, 6)]), -->
<!--       nColLegend = 1, -->
<!--       plotType = "autoplot") -->
<!--     return(p$plot) -->
<!--   }) -->
<!--   gridExtra::grid.arrange( -->
<!--     sur.gene[[1]] + ggplot.them, -->
<!--     sur.gene[[2]] + ggplot.them, -->
<!--     sur.gene[[3]] + ggplot.them, -->
<!--     sur.gene[[4]] + ggplot.them, -->
<!--     sur.gene[[5]] + ggplot.them, -->
<!--     ncol = 3) -->
<!-- } -->
<!-- ``` -->
<!-- # Other normalizations -->
<!-- There are other RNA-seq normalizations including SVAseq, ComBat-seq and RUVg methods that are not specifically designed for normalization, although they can be helpful for that task when the unwanted variation is orthogonal to the biology, something that is rarely known in advance. The same applies to the RUVs method provided in the [RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html) package. Although if there are true replicates (missing from TCGA and most large cancer RNAseq studies), it can be used to normalize RNA-seq datasets.\ -->
<!-- Here, we apply ComBat-seq, RUVg and RUVs on the TCGA READ RNA-seq data and assess their performance. -->
<!-- ## ComBat-seq -->
<!-- We apply the ComBat-seq method on the TCGA READ RNA-seq data. We specify either plates or major times intervals as known sources of batch effects. -->
<!-- ```{r ComBatSeq, message=F, warning=F, results=F} -->
<!-- combat_seq.plates <- sva::ComBat_seq( -->
<!--   counts = raw.count.data,  -->
<!--   batch = read.cancer.se$plate_RNAseq -->
<!--   ) -->
<!-- combat_seq.time <- sva::ComBat_seq( -->
<!--   counts = raw.count.data,  -->
<!--   batch = read.cancer.se$time.points -->
<!--   ) -->
<!-- combat.norm <- list( -->
<!--   combat_seq.plates = combat_seq.plates,  -->
<!--   combat_seq.time = combat_seq.time) -->
<!-- ``` -->
<!-- ### Association between PCs and library size -->
<!-- Figure \@ref(fig:ComBatSeqPca) shows that the Combat-seq method did not remove the library size variation in the TCGA READ RNA-seq data. This results may suggest that the Combat-seq method is not specifically designed for normalization. -->
<!-- ```{r ComBatSeqPca,  message=F, warning=F, results=F, fig.cap='PCA plots of the Combat-seq normalized data coloured by the major time intervals. Top: known source of batch effects is plates. Bottom: known source of batch effects is the major time intervals.'} -->
<!-- pca.combat <- lapply( -->
<!--   c(1:2),  -->
<!--   function(x){ -->
<!--     pca <- .pca(data = combat.norm[[x]], is.log = FALSE) -->
<!--   }) -->
<!-- pp <- lapply( -->
<!--   c(1:2), -->
<!--   function(x){ -->
<!--     pcs <- pca.combat[[x]] -->
<!--     p <- .scatter.density.pc( -->
<!--       pcs = pcs$sing.val$u[,1:3], -->
<!--       pc.var = pcs$var, -->
<!--       group.name = 'Time (years)', -->
<!--       group = read.sampleAnnot$time.points, -->
<!--       color = major.times.colors, -->
<!--       strokeSize = .2, -->
<!--       pointSize = 3, -->
<!--       strokeColor = 'gray30', -->
<!--       alpha = .5) -->
<!--     p -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange,  -->
<!--   c(pp[[1]], pp[[2]], ncol = 4)) -->
<!-- ``` -->
<!-- ## RUVg -->
<!-- In the [RUVSeq](https://bioconductor.org/packages/release/bioc/html/RUVSeq.html) R Bioconductor package there are two methods including RUVg and RUVs that produce both the estimated factors of unwanted variation and normalized counts data. The RUVg method uses a set of negative control genes to estimates the factors of unwanted variation which can be used as covariates in differential gene expression analysis. The authors recommend: “the normalized counts can be used just for exploration, as removing the unwanted factors from the counts can also remove part of a factor of interest.”[Ref](https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf) -->
<!-- However, we apply RUVg on the TCGA READ RNA-seq data. -->
<!-- ```{r ruvg, warning=F,message=F} -->
<!-- ruv.g <- RUVSeq::RUVg( -->
<!--   x = raw.count.data,  -->
<!--   cIdx = negative.control.genes, -->
<!--   k = 12) -->
<!-- ruv.g <- ruv.g$normalizedCounts -->
<!-- ``` -->
<!-- ### Association between gene expression and survival -->
<!-- The survival analyses (Figure \@ref(fig:geneSurRuvg)) show that the RUVg method removes the association between gene expression and survival in the TCGA READ RNA-seq data. -->
<!-- ```{r geneSurRuvg, message=F, error=FALSE, fig.dim=c(8,4),fig.cap='Association between gene expression and overall survival in the RUVg normalized data of the TCGA READ RNA-seq data.'} -->
<!-- selected.genes <- c( -->
<!--   'FBXL14', -->
<!--   'CSGALNACT2') -->
<!-- pp <- lapply( -->
<!--   selected.genes,  -->
<!--   function(x) { -->
<!--     p <- survival_plot( -->
<!--       data = ruv.g, -->
<!--       stratify = 'expr', -->
<!--       annot = as.data.frame(SummarizedExperiment::colData(read.cancer.se)), -->
<!--       scoreCol =  NULL, -->
<!--       gene = x, -->
<!--       covariate = NULL, -->
<!--       isCategoricalCov = FALSE, -->
<!--       timeCol = "OS.time_liu", -->
<!--       eventCol = "OS_liu", -->
<!--       nGroup = 2, -->
<!--       confInt = FALSE, -->
<!--       mainTitle1 = x, -->
<!--       ylabel = "Survival", -->
<!--       cols = brewer.pal(9, "Set1")[c(2, 3)], -->
<!--       nColLegend = 1, -->
<!--       plotType = "autoplot" -->
<!--     ) -->
<!--     p$plot -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange,  -->
<!--   c(pp, ncol = 2)) -->
<!-- ``` -->
<!-- ## RUVs -->
<!-- Another method in the RUVSeq package is RUVs. This method uses replicate/negative control samples for which the covariates of interest are constant and negative control genes to estimate the factors of unwanted variation. RUVs produces both the factors of unwanted variation that can be used for differential expression analysis and normalized counts data. Generally, this method is used in situations where the genuine replicate/negative control samples are available.  For example, an RNA-Seq study that involves control and treatment samples.\ -->
<!-- However, we apply RUVs on the TCGA READ RNA-Seq data and use the consensus molecular subtypes as replicate/negative control samples for RUVs. We  also use a set of negative control genes that was used for the RUV-III-PRPS normalization.  -->
<!-- ```{r RuvsNorm, message=F, warning=F, error=F} -->
<!-- rep.samples <- RUVSeq::makeGroups(read.cancer.se$cms.cancer.ruv) -->
<!-- ruv.s <- RUVSeq::RUVs( -->
<!--   x = raw.count.data,  -->
<!--   cIdx = negative.control.genes, -->
<!--   k = 12,  -->
<!--   scIdx = rep.samples) -->
<!-- ruv.s <- ruv.s$normalizedCounts -->
<!-- ``` -->
<!-- ### Expression levels of several genes -->
<!-- Figure \@ref(fig:ruvsExamples) shows that RUVs converts the expression levels of some genes with reasonable expression to zero across all samples. For examples, the log2 expression of the LARP7 gene is between 10 to 12.5 in the TCGA raw counts data, while the expression of this gene is almost zero in the RUVs normalized data. These results may suggest that  RUVs  is not designed for RNA-seq normalization in situations where true technical replicates are not available. -->
<!-- ```{r ruvsExamples, message=F, warning=F, error=F, fig.cap='Expression levels of several genes in the TCGA READ raw counts and the RUVs normalized data.'} -->
<!-- selected.genes <- c( -->
<!--   'RAB18',  -->
<!--   'LARP7', -->
<!--   'PTPN14', -->
<!--   'CSGALNACT2') -->
<!-- pp <- lapply( -->
<!--   selected.genes,  -->
<!--   function(g){ -->
<!--     df <- data.frame( -->
<!--       Raw.count = raw.count.data[g ,], -->
<!--       RUVs = ruv.s[g ,], -->
<!--       sample = c(1:166)) %>% -->
<!--       tidyr::pivot_longer( -->
<!--         -sample, -->
<!--         values_to = 'expr', -->
<!--         names_to = 'genes') -->
<!--     ggplot(df, aes(x = sample, y = log2(expr + 1))) + -->
<!--       geom_point() + -->
<!--       ylab('Gene expression') + -->
<!--       xlab('Samples') + -->
<!--       ggtitle(g) + -->
<!--       facet_wrap(~ genes) + -->
<!--       theme( -->
<!--         panel.background = element_blank(), -->
<!--         axis.line = element_line(colour = 'black', size = .85), -->
<!--         axis.title.x = element_text(size = 14), -->
<!--         axis.title.y = element_text(size = 14), -->
<!--         axis.text.x = element_text( -->
<!--           size = 12, -->
<!--           angle = 25, -->
<!--           vjust = 1, -->
<!--           hjust = 1), -->
<!--         axis.text.y = element_text(size = 12), -->
<!--         legend.text = element_text(size = 10), -->
<!--         legend.title = element_text(size = 14), -->
<!--         strip.text.x = element_text(size = 10)) -->
<!--   }) -->
<!-- do.call( -->
<!--   gridExtra::grid.arrange, -->
<!--   pp) -->
<!-- ``` -->
<!-- # R Session information -->
<!-- ```{r} -->
<!-- options(max.print = 10^4) -->
<!-- sessionInfo() -->
<!-- ``` -->