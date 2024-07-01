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

sc.metabolism.customized <- function(obj, method = "AUCell",input_pathway, imputation = F,ncores=20,input_dir) {
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
    paste(x["Pathway"],NA, gsub(",", "\t", x["Genes"]), sep = "\t")
  })
  # 输出到控制台，或写入到文件
  #cat(gmt_lines, sep = "\n")
  # 可选：将结果写入文件
  writeLines(gmt_lines, paste0(input_dir,"pathways_select.gmt"))
  gmtFile=paste0(input_dir,"pathways_select.gmt")
  geneSets <- getGmt(gmtFile) #signature read
  geneSets_list <- lapply(geneSets, geneIds)
  names(geneSets_list)=names(geneSets)
if (!file.exists(paste0(input_dir,"pathways_select.gmt"))) {
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
    # 示例基因表达矩阵
    # scaled_counts <- your_scaled_counts_matrix  # 请替换为实际的标准化基因表达矩阵
    # 1. 计算每个基因在多少个细胞中被检测到（非零表达）
    gene_detection_counts <- rowSums(scaled_counts > 0)
    # 2. 选择在至少一定比例的细胞中检测到的基因
    # 例如，在至少0.10%的细胞中检测到的基因
    min_cells <- 0.001 * ncol(scaled_counts)
    aproved_genes <- gene_detection_counts >= min_cells
    # 3. 输出结果
    genes_used <-names(aproved_genes)[which(aproved_genes)]

    # 检查每个基因集的基因是否存在于 genes_used 向量中
    gene_existence_check <- lapply(geneSets_list, function(gene_list) {
      all(gene_list %in% genes_used)
    })
    # 打印所有基因都不在 genes_used 中的基因集的名称
    all_false_gene_set_names <- names(geneSets_list)[unlist(gene_existence_check)]

    if (length(all_false_gene_set_names) == 0) {
      message <- "All pathways passed testing."
      writeLines(message, "./pathways_testing.txt")
      }
    if (length(all_false_gene_set_names) > 0) {
      if (length(all_false_gene_set_names) == length(names(geneSets_list))) {
        # 如果 all_false_gene_set_names 的长度等于 geneSets_list 的名称数量，则停止运行并输出特定消息
        message <- "All pathway(s)' genes do not meet the optimal range of the VISION algorithm, we have to stop. Please change the algorithm or pathway for calculation."
        writeLines(message, "./pathways_testing.txt")
        print(message)
        return(F)
      } else {
        # 如果 all_false_gene_set_names 的长度大于 0 但不等于 geneSets_list 的名称数量，则输出另一条消息
        message <- "These/This pathway(s)' genes do not meet the optimal range of the VISION algorithm:\n"
        message <- paste0(message, paste(all_false_gene_set_names, collapse = "\n"))
        writeLines(message, "pathways_testing.txt")
      }
    }

    vis <- Vision(scaled_counts, signatures = gmtFile,min_signature_genes=0,sig_gene_threshold = 0.001)
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

    # 获取排名数据中的基因ID
    rankings_ids <- rownames(cells_rankings)

    available_ratios <- sapply(geneSets_list, function(genes) {
      available_genes <- genes %in% rankings_ids
      ratio <- sum(available_genes) / length(genes)
      return(ratio * 100)  # 返回百分比))
    })
    names(available_ratios)=names(geneSets)

    gene_existence_check <- available_ratios<20
    # 打印所有基因都不在 genes_used 中的基因集的名称
    all_false_gene_set_names <- names(geneSets)[gene_existence_check]

    if (length(all_false_gene_set_names) == 0) {
      message <- "All pathways passed testing."
      writeLines(message, "./pathways_testing.txt")
    }
    if (length(all_false_gene_set_names) > 0) {
      if (length(all_false_gene_set_names) == length(names(geneSets_list))) {
        # 如果 all_false_gene_set_names 的长度等于 geneSets_list 的名称数量，则停止运行并输出特定消息
        message <- "All pathways' genes do not meet the optimal range of the AUCell algorithm, we have to stop. Please change the algorithm or pathway for calculation."
        writeLines(message, "./pathways_testing.txt")
        print(message)
        return(F)
      } else {
        # 如果 all_false_gene_set_names 的长度大于 0 但不等于 geneSets_list 的名称数量，则输出另一条消息
        message <- "These pathways' genes do not meet the optimal range of the AUCell algorithm:\n"
        message <- paste0(message, paste(all_false_gene_set_names, collapse = "\n"))
        writeLines(message, "pathways_testing.txt")
      }
    }





    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc


    signature_exp <- data.frame(getAUC(cells_AUC))


  }

  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    # 设置并行计算参数

    genes_used <-rownames(countexp2)
    # 检查每个基因集的基因是否存在于 genes_used 向量中
    gene_existence_check <- lapply(geneSets_list, function(gene_list) {
      all(gene_list %in% genes_used)
    })
    # 打印所有基因都不在 genes_used 中的基因集的名称
    all_false_gene_set_names <- names(geneSets_list)[unlist(gene_existence_check)]

    if (length(all_false_gene_set_names) == 0) {
      message <- "All pathways passed testing."
      writeLines(message, "./pathways_testing.txt")
    }
    if (length(all_false_gene_set_names) > 0) {
      if (length(all_false_gene_set_names) == length(names(geneSets_list))) {
        # 如果 all_false_gene_set_names 的长度等于 geneSets_list 的名称数量，则停止运行并输出特定消息
        message <- "All pathways' genes do not present in the ssGSEA, we have to stop. Please change the algorithm or pathway for calculation."
        writeLines(message, "./pathways_testing.txt")
        print(message)
        return(F)
      } else {
        # 如果 all_false_gene_set_names 的长度大于 0 但不等于 geneSets_list 的名称数量，则输出另一条消息
        message <- "These pathways' genes do not meet the optimal range of the AUCell algorithm:\n"
        message <- paste0(message, paste(all_false_gene_set_names, collapse = "\n"))
        writeLines(message, "pathways_testing.txt")
      }
    }




    bpparam <- MulticoreParam(workers = ncores)


    # 创建 ssgseaParam 参数对象
    ssgsea_param <- ssgseaParam(expr = as.matrix(countexp2), geneSets = geneSets, minSize = 1)

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

