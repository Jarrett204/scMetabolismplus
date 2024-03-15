#' Single-Cell Metabolism Analysis
#'
#' This function performs metabolism analysis on single-cell gene expression data.
#' It supports multiple methods for quantifying metabolism activity, including VISION, AUCell, and GSVA.
#' Users can choose between KEGG and REACTOME metabolism pathways and optionally apply imputation to handle missing data.
#'
#' @param countexp A matrix of count expression data where rows are genes and columns are single cells.
#' @param method The method to use for quantifying metabolism activity. One of "VISION", "AUCell", "ssGSEA", or "GSVA". Default is "VISION".
#' @param imputation Logical indicating whether to apply imputation to the count expression data. Default is FALSE.
#' @param ncores The number of cores to use for parallel processing. Default is 2.
#' @param metabolism.type The type of metabolism pathway to analyze. One of "KEGG" or "REACTOME". Default is "KEGG".
#'
#' @return A data frame with metabolism activity scores for each cell.
#'
#' @examples
#' # Example dataset (not real data)
#'
#' # Perform metabolism analysis using default settings
#' results <- sc.metabolism(countexp)
#'
#' @export
#'
#' @references
#' Yingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. <https://pubmed.ncbi.nlm.nih.gov/34417225/>
#' George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. <https://doi.org/10.1101/397588>
sc.metabolism <- function(countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG") {
  # Function body
}



sc.metabolism <- function(countexp, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG") {

  #signatures_KEGG_metab <- "./data/KEGG_metabolism_nc.gmt"
  #signatures_REACTOME_metab <- "./data/REACTOME_metabolism.gmt"

  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolismplus")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "scMetabolismplus")

  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}

  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")

    #Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588
    #Github: https://github.com/KlugerLab/ALRA
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")


    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }

  #signature method
  cat("Start quantify the metabolism activity...\n")

  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)

    options(mc.cores = ncores)

    vis <- analyze(vis)

    signature_exp<-data.frame(t(vis@SigScores))
  }

  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }

  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }

  #GSVA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }

  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  signature_exp
}
