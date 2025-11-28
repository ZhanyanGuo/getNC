#' pbmc_small
#'
#' A dataset of Peripheral Blood Mononuclear Cells (PBMC)
#' freely available from 10X Genomics.
#' There are 2,700 single cells that were sequenced on the Illumina NextSeq 500.
#' This is coming from the cellranger pipeline from 10X, returning a unique
#' molecular identified (UMI) count matrix.
#' For example dataset we take the first 1000 cells.
#' For detailed generating step, see data-raw/pbmc10x.R
#'
#' @format gene-by-cell UMI count matrix (dgCMatrix)RDS-compressed sparse matrix
#' @source Seurat https://satijalab.org/seurat/articles/pbmc3k_tutorial
#' @references
#' Hao, Y. et al. (2021). Integrated analysis of multimodal single-cell data.
#'   \emph{Cell} 184(13):3573–3587.
#' 
"pbmc_small"