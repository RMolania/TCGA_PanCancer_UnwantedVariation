# Load Data
# Due to data size issues the actual raw data is not included as part of the package. A condensed
# version of the original dataset containing 100 observations is used as part of the package functionality.
brca.se <- base::readRDS(here::here("data_raw","ForShiny_TCGA_SummarizedExperiment_HTseq_BRCA.rds"))

# Data Subset
brca.data <- brca.se[1:100, ]

#
#library(SummarizedExperiment)
# Gene level data
gene.annot <-  as.data.frame(SummarizedExperiment::rowData(brca.data))
# Sample level data
sample.info <-  as.data.frame(SummarizedExperiment::colData(brca.data))
# Raw and normalized values
raw.count <- as.data.frame(SummarizedExperiment::assay(brca.data, 'HTseq_counts'))
fpkm <- as.data.frame(SummarizedExperiment::assay(brca.data, 'HTseq_FPKM'))
fpkm.uq <- as.data.frame(SummarizedExperiment::assay(brca.data, 'HTseq_FPKM.UQ'))

# Save Data in Package
usethis::use_data(brca.data, compress = "xz", overwrite = TRUE)

