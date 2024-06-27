#' scMetabolismplus
#'
#' scMetabolism
#' @param obj seruat 对象
#' @param ncores
#' @param input_pathway #输入的通路和所属类别
#'  使用的核心数量，默认为 20
#' @keywords scMetabolismplus
#' @examples
#' sc.metabolism.customized()
#' @export sc.metabolism.customized

sc.metabolism.customized <- function(obj, method = "AUCell", imputation = F,ncores=20) {
  library(GSEABase)
  library(dplyr)
  countexp<-obj@assays$RNA@counts
  countexp<-data.frame(as.matrix(countexp))
  signatures_metab=system.file("extdata/", "all_metabolism.csv", package = "scMetabolismplus")%>%read.csv()
  select_sig <- merge(input_pathway, signatures_metab, by = c("Pathway", "Category"))
  select_sig <- select_sig %>%
    # 在第一列添加新列，组合 Category 和 Pathway
    mutate(Combined = paste(Category, Pathway, sep = "_")) %>%
    # 移动 Combined 列到第一列位置
    dplyr::select(Combined, everything()) %>%
    # 删除原始的 Pathway 和 Category 列
    dplyr::select(-Pathway, -Category)
  colnames(select_sig)[1]="Pathway"
  # 转换为GMT格式
  gmt_lines <- apply(select_sig, 1, function(x) {
    paste(x["Pathway"], gsub(",", "\t", x["Genes"]), sep = "\t")
  })
  # 输出到控制台，或写入到文件
  #cat(gmt_lines, sep = "\n")
  # 可选：将结果写入文件
  writeLines(gmt_lines, "./pathways_select.gmt")
  gmtFile="./pathways_select.gmt"
if (!file.exists("./pathways_select.gmt")) {
  # 如果文件不存在，停止执行并返回错误信息
  stop("Error: The file './pathways_select.gmt' does not exist. Please check the file path.")
}
  #imputation
  if (imputation) {
    cat("Start imputation...\n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  } else {
    countexp2 <- countexp
  }

  #signature method
  cat("Start quantify the path activity...\n")

  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    # 检查数据中NA和零值的数量
    options(mc.cores =ncores)
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
    # 设置并行计算参数
    bpparam <- MulticoreParam(workers = ncores)
    geneSets <- getGmt(gmtFile) #signature read
    # 创建 ssgseaParam 参数对象
    ssgsea_param <- ssgseaParam(expr = as.matrix(countexp2), geneSets = geneSets)

    # 计算 GSVA 富集分数
    gsva_es <- gsva(ssgsea_param, BPPARAM = bpparam)
    signature_exp<-data.frame(gsva_es)
  }


  cat("\ Prceeding!
      \n\n")
  obj@meta.data$Cancer <- Cancer
  obj@assays$METABOLISM$score<-signature_exp
  obj
}

