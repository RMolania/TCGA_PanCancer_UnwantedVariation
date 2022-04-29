This repository contains R scripts, RMarkdown reports, R package and R shiny app associated with the manuscript: [Molania R et al., Removing Unwanted Variation from Larg-scale Cancer RNA-seq Data. bioRxiv, 2021.]()

Available folders contain the below files:

- **Scripts**: contains all Rmarkdown (.Rmd) scripts containing the main analyses codes. 
- **Vigettes**: Two comprehensive vignettes showing all the steps to identify and remove unwanted variation from the TCGA READ abd BRCA RNA-seq data 
- **TCGA_RShiny_**: The shiny app for exploration sources of unwanted variation in all TCGA RNA-seq datasets.

Data availability
We have created a Summarized experiment objects contacting expression data (raw counts, FPKM and FPKM.UQ), clinical and batch information, and gene annotations for all the TCGA RNA-Seq data. These files are deposited here (https://zenodo.org/record/6326542#.YlN56y8Rquo). The datasets that are required for the vignettes can be found here (https://zenodo.org/record/6392171#.YlN6Yi8Rquo). The RUV-III normalized data of the TCGA READ, COAD and BRCA RNA-Seq datasets are deposited (https://zenodo.org/record/6459560#.YldJ4S8Rquo). The TCGA gene expression microarray data was obtained from (https://gdac.broadinstitute.org), data version 2016/01/28. The breast cancer laser captured microdissected microarray and RNA-seq datasets were downloaded from the NCBI Gene Expression Omnibus with accession number GSE78958, GSE96058 and GSE81538. 

[All the TCGA RNA-seq data](https://zenodo.org/record/6326542#.Ymu9wS8Rquo)

Vigettes:
[The TCGA READ RNA-seq data]()


Please email molania.r@wehi.edu.au if you have any questions.


*Unwanted variation in TCGA RNA-seq datasets.*

<img src="./TCGA_Main1.png">

Unwanted variation in individual TCGA RNA-seq datasets. A) Distribution of (log2) library size coloured by years for each cancer type. The year information was not available for the LAML RNA-seq data. The library sizes are calculated after removing lowly expressed genes for each cancer type. B) R2 obtained from linear regression between the first, first and second, etc., cumulatively to the fifth principal component (PC) and library size (first panel), tumor purity (second panel), and Relative Log Expression (RLE) median (third panel) in the raw count, FPKM and FPKM.UQ normalized datasets. The fourth panel shows the vector correlation between the first 5 PC cumulatively and sequencing-plates in the datasets. Grey colour indicates that samples were profiled across a single plate. C) Spearman correlation coefficients between individual gene expression levels and library size (first panel), tumor purity (second panel) and RLE median (third panel) in the datasets. The fourth panel shows log2 F statistics obtained from ANOVA of gene expression levels by the factor: plate variable. Plates with less than three samples were excluded from the analysis. ANOVA was not possible for cancer types whose samples were profiled using a single plate.