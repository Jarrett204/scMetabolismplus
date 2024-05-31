#' scMetabolismplus
#'
#' scMetabolism
#' @param obj seruat 对象
#' @param Cancer 数据属于的癌种 对象
#' @keywords scMetabolismplus
#' @examples
#' sc.metabolism.Seurat.pathway()
#' @export sc.metabolism.Seurat.pathway


sc.metabolism.Seurat.pathway <- function(obj, method = "AUCell", imputation = F,Cancer="BRCA", metabolism.type = "KEGG",ncore=20) {
  library(GSEABase)
  countexp<-obj@assays$RNA@counts

  countexp<-data.frame(as.matrix(countexp))

  #signatures_KEGG_metab <- "./data/KEGG_metabolism_nc.gmt"
  #signatures_REACTOME_metab <- "./data/REACTOME_metabolism.gmt"

  signatures_KEGG_metab <- system.file("extdata/KEGG_path", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_REACTOME_metab <- system.file("extdata/Reactome_path", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_GO_metab <- system.file("extdata/GO_path", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures_HMDB_metab <- system.file("extdata/HMDB_path", paste0(Cancer,".gmt"), package = "scMetabolismplus")
  signatures__metab <- system.file("extdata/GO_path", paste0(Cancer,".gmt"), package = "scMetabolismplus")


  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "Reactome")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  if (metabolism.type == "GO")  {gmtFile<-signatures_GO_metab; cat("Your choice is: GO\n")}
  if (metabolism.type == "HMDB")  {gmtFile<-signatures_HMDB_metab; cat("Your choice is: HMDB\n")}
  if (metabolism.type == "GO")  {gmtFile<-signatures_GO_metab; cat("Your choice is: GO\n")}

  file.exists(gmtFile)
  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }

  #signature method
  cat("Start quantify the path activity...\n")

  #VISION 暂时看不好使
  if (method == "VISION") {
   library(VISION)
  n.umi <- colSums(countexp2)
  scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
  vis <- Vision(scaled_counts, signatures = gmtFile)
  # 检查数据中NA和零值的数量
  options(mc.cores = 1)

  vis <- analyze(vis)

  signature_exp<-data.frame(t(vis@SigScores))
  }

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

  cat("\ Thanks to:XUDONG,MENGDI,YILIN,WANGQI,RENHE,JIANYU3.26.1\
      You guys are brilliant!
      \n\n")
  obj@meta.data$Cancer <- Cancer
  obj@meta.data$dataset <- metabolism.type
  obj@assays$METABOLISM$score<-signature_exp
  obj
}

