This repository contains R scripts, RMarkdown reports, R package and R shiny app associated with the manuscript [Molania R et al, Removing Unwanted Variation from Larg-scale Cancer RNA-seq Data. bioRxiv, 2021.]()

Available folders contain the below files:

- **report_script**: several Rmarkdown (.Rmd) reports and R scripts containing the main analyses codes. 
- **package**: R package regenrating the results of this data
- **shinyApp**: The shiny app for data exploration and normalisation

Please email molania.r@wehi.edu.au if you have any questions.


*Unwanted variation in individual TCGA RNA-seq datasets.*

<img src="./TCGA_Main1.png">

Unwanted variation in individual TCGA RNA-seq datasets. A) Distribution of (log2) library size coloured by years for each cancer type. The year information was not available for the LAML RNA-seq data. The library sizes are calculated after removing lowly expressed genes for each cancer type. B) R2 obtained from linear regression between the first, first and second, etc., cumulatively to the fifth principal component (PC) and library size (first panel), tumor purity (second panel), and Relative Log Expression (RLE) median (third panel) in the raw count, FPKM and FPKM.UQ normalized datasets. The fourth panel shows the vector correlation between the first 5 PC cumulatively and sequencing-plates in the datasets. Grey colour indicates that samples were profiled across a single plate. C) Spearman correlation coefficients between individual gene expression levels and library size (first panel), tumor purity (second panel) and RLE median (third panel) in the datasets. The fourth panel shows log2 F statistics obtained from ANOVA of gene expression levels by the factor: plate variable. Plates with less than three samples were excluded from the analysis. ANOVA was not possible for cancer types whose samples were profiled using a single plate.