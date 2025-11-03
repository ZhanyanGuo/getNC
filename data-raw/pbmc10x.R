# data-raw/pbmc10x.R

library(Seurat)

pbmc <- Read10X(data.dir = "data-raw//filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc,project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc_small <- subset(pbmc, cells = colnames(pbmc)[1:1000])

usethis::use_data(pbmc_small, overwrite = TRUE)
