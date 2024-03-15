#' scMetabolismplus
#'
#' scMetabolism
#' @param obj seruat 对象
#' @param Cancer 数据属于的癌种 对象
#' @keywords scMetabolismplus
#' @examples
#' sc.metabolism.Seurat()
#' @export sc.metabolism.Seurat


sc.metabolism.Seurat <- function(obj, method = "AUCell", imputation = F,Cancer="BRCA", metabolism.type = "KEGG") {

  countexp<-obj@assays$RNA@counts

  countexp<-data.frame(as.matrix(countexp))

  #signatures_KEGG_metab <- "./data/KEGG_metabolism_nc.gmt"
  #signatures_REACTOME_metab <- "./data/REACTOME_metabolism.gmt"

  signatures_KEGG_metab <- system.file("extdata/KEGG/", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_REACTOME_metab <- system.file("extdata/Reactome/", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_GO_metab <- system.file("extdata/GO/", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_HALLMARK_metab <- system.file("extdata/Hallmark/", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_HMDB_metab <- system.file("extdata/HMDB/", paste0(Cancer,".gmt"), package = "scMetabolismplus")


  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  if (metabolism.type == "GO")  {gmtFile<-signatures_GO_metab; cat("Your choice is: GO\n")}
  if (metabolism.type == "Hallmark")  {gmtFile<-signatures_HALLMARK_metab; cat("Your choice is: Hallmark\n")}
  if (metabolism.type == "HMDB")  {gmtFile<-signatures_HMDB_metab; cat("Your choice is: HMDB\n")}

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

  #VISION 暂时看不好使
  #if (method == "VISION") {
  #  library(VISION)
  #  n.umi <- colSums(countexp2)
  #  scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
  #  vis <- Vision(scaled_counts, signatures = gmtFile)

  #  options(mc.cores = ncores)

  #  vis <- analyze(vis)

   # signature_exp<-data.frame(t(vis@SigScores))
  #}

  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), plotStats=F) #rank
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

  obj@assays$METABOLISM$score<-signature_exp
  obj
}

