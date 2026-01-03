#' scMetabolism
#' @param obj objective
#' @param pathway objective pathway
#' @param phenotype phenotype
#' @keywords scMetabolism
#' @examples
#' DimPlot.metabolism()
#' DotPlot.metabolism()
#' BoxPlot.metabolism()
#' PathUmp.metabolism
#' @export DimPlot.metabolism
#' @export DotPlot.metabolism
#' @export BoxPlot.metabolism
#' @export PathUmp.metabolism

library(ggplot2)
library(wesanderson)
library(data.table)
library(rsvd)
library(dplyr)

DimPlot.metabolism <- function(obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, size= 0.01,Width=6,Height=5,dynamic=F){

  cat("Establishing connection\n\n")
  library(wesanderson)
  library(ggplot2)
  library(progress)

  base_colors <- c("#4D457E", "#CD506B", "#E9DF45")
  # 创建渐变色调色板函数
  pal <- colorRampPalette(base_colors)(100)

  # Initialize progress bar with total steps
  total_steps <- length(pathway) * 2  # One for creating data, one for plotting
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )

  # umap
  if (dimention.reduction.type == "umap"){

    if (dimention.reduction.run == T) obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:40)
    umap.loc <- obj@reductions$umap@cell.embeddings

    row.names(umap.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score

    signature_ggplot2 <- data.frame()
    for (input.pathway in pathway) {
      signature_ggplot <- data.frame(umap.loc, Pathway = input.pathway, Score = unlist(signature_exp[input.pathway, , drop = T]))
      if(input.pathway == pathway[1]){
        signature_ggplot2 <- signature_ggplot
      } else {
        signature_ggplot2 <- rbind(signature_ggplot2, signature_ggplot)
      }
      pb$tick()  # Update progress bar after data processing
    }

    if(!dynamic){
    output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dimplot")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    }
    Dimplot_all <- list()
    for (input.pathway in pathway) {
      pathway_data <- subset(signature_ggplot2, Pathway == input.pathway)
      plot <- ggplot(data = pathway_data, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
        geom_point(size = size) +
        scale_color_gradientn(colours = pal) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        ggtitle(input.pathway) +
        theme_bw()
      if(!dynamic){
      ggsave(filename = paste0(output_dir,"/",input.pathway, ".png"), plot = plot, width = Width, height = Height)}
      pb$tick()  # Update progress bar after plotting
      Dimplot_all[[input.pathway]]=plot
    }
    #print(plot)  # 打印或保存plot
  }

  # tsne
  if (dimention.reduction.type == "tsne"){
    if (dimention.reduction.run == T) obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:40)
    tsne.loc <- obj@reductions$tsne@cell.embeddings

    row.names(tsne.loc) <- colnames(obj)
    signature_exp <- obj@assays$METABOLISM$score

    signature_ggplot2 <- data.frame()
    for (input.pathway in pathway) {
      signature_ggplot <- data.frame(tsne.loc, Pathway = input.pathway, Score = unlist(signature_exp[input.pathway, , drop = T]))
      if(input.pathway == pathway[1]){
        signature_ggplot2 <- signature_ggplot
      } else {
        signature_ggplot2 <- rbind(signature_ggplot2, signature_ggplot)
      }
      pb$tick()  # Update progress bar after data processing
    }
    if(!dynamic){
    output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dimplot")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    }
    for (input.pathway in pathway) {
      pathway_data <- subset(signature_ggplot2, Pathway == input.pathway)
      plot <- ggplot(data = pathway_data, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
        geom_point(size = size) +
        scale_color_gradientn(colours = pal) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        ggtitle(paste("Pathway:", input.pathway)) +
        theme_bw()
      ggsave(filename = paste0(output_dir, "/", "plot_", input.pathway, ".png"), plot = plot, width = Width, height = Height)
      pb$tick()  # Update progress bar after plotting
    }
  }
  plot
  return(Dimplot_all)
}

DotPlot.metabolism <- function(obj, pathway, phenotype, norm = "y",Width=6,Height=4,dynamic=F){
  library(viridis)
  library(progress)
  library(dplyr)
  library(ggplot2)

  input.norm = norm
  input.pathway <- pathway
  input.parameter <- phenotype
  pal <- viridis::viridis(100)

  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score

  cat("Starting Dotplots\n\n")

  metadata[,input.parameter] <- as.character(metadata[,input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway,])

  gg_table <- c()

  # Calculate total steps for progress bar
  total_steps_pathway <- length(input.pathway)
  input.group.x <- unique(as.character(metadata[,input.parameter]))
  input.group.y <- unique(as.character(input.pathway))
  total_steps_median <- length(input.group.x) * length(input.group.y)
  if (input.norm == "y") {
    total_steps_norm <- length(input.group.y)
  } else if (input.norm == "x") {
    total_steps_norm <- length(input.group.x)
  } else {
    total_steps_norm <- 0
  }
  total_steps <- total_steps_pathway + total_steps_median + total_steps_norm

  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )

  # Process pathways
  for (i in 1:length(input.pathway)){
    gg_table <- rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
    pb$tick()
  }
  gg_table <- data.frame(gg_table)

  # Get median values
  gg_table_median <- c()
  for (x in 1:length(input.group.x)){
    for (y in 1:length(input.group.y)){
      gg_table_sub <- subset(gg_table, gg_table[,1] == input.group.x[x] & gg_table[,2] == input.group.y[y])
      gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], median(as.numeric(as.character(gg_table_sub[,3])))))
      pb$tick()
    }
  }
  gg_table_median <- data.frame(gg_table_median)
  gg_table_median[,3] <- as.numeric(as.character(gg_table_median[,3]))

  # Normalize data
  gg_table_median_norm <- c()
  range01 <- function(x){(x - min(x)) / (max(x) - min(x))}

  if (input.norm == "y"){
    for (y in 1:length(input.group.y)){
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[,2] == input.group.y[y])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, gg_table_median_sub)
      pb$tick()
    }
  }

  if (input.norm == "x"){
    for (x in 1:length(input.group.x)){
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[,1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, gg_table_median_sub)
      pb$tick()
    }
  }

  if (input.norm == "na") gg_table_median_norm <- gg_table_median

  gg_table_median_norm <- data.frame(gg_table_median_norm)
  gg_table_median_norm[,3] <- as.numeric(as.character(gg_table_median_norm[,3]))
  gg_table_median_norm <- dplyr::filter(gg_table_median_norm, !is.na(X3))
  library(pheatmap)
  if(!dynamic){
  output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dotplot")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  }
  # 使用 pheatmap 进行行聚类
  # 将长格式数据转换为宽格式
  wide_df <- gg_table_median_norm %>%
    ungroup() %>%
    select(X1, X2, X3) %>%
    tidyr::pivot_wider(names_from = X1, values_from = X3) %>%
    as.data.frame()

  # --- 2. 提取通路名称并创建纯数值矩阵 ---
  # wide_df 的第一列是 X2 (通路名)，其余列是样本/分组
  row_names_vec <- as.character(wide_df$X2)
  wide_matirx <- wide_df[, -1, drop = FALSE] # 👈 关键：drop = FALSE 保证单行不坍塌

  # 转为纯数值并重新赋行名
  wide_matirx <- as.data.frame(lapply(wide_matirx, as.numeric))
  rownames(wide_matirx) <- row_names_vec

  # --- 3. 健壮的聚类判断逻辑 ---
  if (nrow(wide_matirx) > 1) {
    # 只有多个通路时才聚类
    clustering <- pheatmap(as.matrix(wide_matirx), silent = TRUE)
    row_order <- rownames(wide_matirx)[clustering$tree_row$order]
  } else {
    # 只有一个通路时，顺序就是它自己
    row_order <- row_names_vec
    cat("Only one pathway detected, skip clustering.\n")
  }

  # --- 4. 这里的变量赋值确保后面 write.csv 不会报错 ---
  # 后面的代码会用到 row_order 来设置 factor levels
  gg_table_median_norm$X2 <- factor(gg_table_median_norm$X2, levels = row_order)
  gg_table_median_norm$X1 <- as.factor(gg_table_median_norm$X1)



   if(length(gg_table_median_norm$X1) == 0){
    cat("No pathway qualified\n\n")
  } else {
    plot_dot <- ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[,1], y = gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
      geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[,3])) +
      ylab("Metabolic Pathway") + xlab(input.parameter) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_line(color = "gray", size = 0.2),
            panel.grid.major = element_line(color = "gray", size = 0.2)) +
      scale_color_gradientn(colours = pal) +
      labs(color = "Value", size = "Value") +
      NULL

    cat("You can paste pathways here for boxplot and so on\n\n")
    print(gg_table_median_norm$X2 %>% unique())
    print(plot_dot)

    # Create output directory and save plot
    if(!dynamic){
    ggsave(filename = paste0(output_dir, "/", "plot_Dot", ".png"), plot = plot_dot, width =Width, height = Height,limitsize = FALSE)
    write.csv(wide_matirx,paste0(output_dir, "/", "wide_matrix", ".csv"),row.names = T)}
    result <- list(
      plot = plot_dot,
      pathway = row_order,
      level=levels(gg_table_median_norm[,1])
    )
    return(result)
  }
}

BoxPlot.metabolism <- function(obj, pathway, phenotype,levels,Width=6,Height=4,dynamic=F){
  library(wesanderson)
  library(RColorBrewer)
  library(ggsci)
  library(progress)

  input.pathway <- pathway
  input.parameter <- phenotype

  cat("Start BoxPlot\n\n")

  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  metadata[,input.parameter] <- as.character(metadata[,input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway,])

  # Initialize progress bar with total steps
  total_steps <- length(input.pathway) * 2  # One for creating data, one for plotting
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )

  # Arrange large table
  gg_table <- c()
  input.pathway <- as.character(input.pathway)
  for (i in 1:length(input.pathway)){
    gg_table <- rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
    pb$tick()  # Update progress bar after data processing
  }
  gg_table <- data.frame(gg_table)
  colnames(gg_table) <- c("cluster", "Pathway", "Score")
  gg_table$Score <- as.numeric(as.character(gg_table$Score))
  print(head(gg_table))

  # Combine multiple color palettes manually
  colors <- c(pal_jama()(4),pal_npg()(5),pal_jco()(5),pal_aaas()(5),pal_lancet()(5))
  if(!dynamic){
  output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Boxplot")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  }
  result <- list()
  for (select.pathway in input.pathway) {
    pathway_data <- subset(gg_table, Pathway == select.pathway)
    print(head(pathway_data))

    stats <- pathway_data %>%
      group_by(cluster) %>%
      summarise(
        ymin = min(Score),
        lower = quantile(Score, 0.25),
        middle = median(Score),
        upper = quantile(Score, 0.75),
        ymax = max(Score)
      )
    print(head(stats))
    pathway_data$cluster=factor(pathway_data$cluster,levels = levels)
    plot_box <- ggplot(data = pathway_data, aes(x = cluster, y = Score, fill = cluster)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      ylab("Metabolic Pathway") +
      xlab("Cell type") +
      theme_bw() +
      ggtitle(select.pathway) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_line(),
            plot.title = element_text(size = 3, face = "bold") ,
            panel.grid.major = element_line()) +
      labs(fill = "Cell type") +
      scale_fill_manual(values = colors)
      result[[select.pathway]] = plot_box
    # Calculate dynamic height based on the number of pathways
    if(!dynamic){
    ggsave(filename = paste0(output_dir, "/",select.pathway, ".png"), plot = plot_box, width = Width, height = Height)
    }
    pb$tick()  # Update p rogress bar after plotting

  }
  plot_box
  return(result)
}

PathUmp.metabolism <- function(obj, phenotype,n.neighbors=3,threshold = 3, top_n = 5,size=5,Width=5,Height=5,dynamic=F) {
  library(progress)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(ggrepel)
  library(ggforce)
  if (n.neighbors < 2) {
    stop("Error: 'n.neighbors' must be at least 2.")
  }
  t_test_t_score <- function(data, cluster_id, pathway_id) {
    current_cluster_scores <- data %>%
      filter(cluster == cluster_id & Pathway == pathway_id) %>%
      pull(Score)
    other_cluster_scores <- data %>%
      filter(cluster != cluster_id & Pathway == pathway_id) %>%
      pull(Score)

    # 检查每组的观察值数量是否足够
    if (length(current_cluster_scores) < 2) {
      message(paste("Cluster", cluster_id, "and Pathway", pathway_id, "has less than 2 current cluster scores."))
      return(NA)  # 如果任一组的观察值少于2，则返回NA
    }

    if (length(other_cluster_scores) < 2) {
      message(paste("Cluster", cluster_id, "and Pathway", pathway_id, "has less than 2 other cluster scores."))
      return(NA)  # 如果任一组的观察值少于2，则返回NA
    }

    # 进行t检验并返回t检验统计量
    t_test_result <- t.test(current_cluster_scores, other_cluster_scores)
    return(t_test_result$statistic)
  }

  # 检查并创建输出目录
  if(!dynamic){
  output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Path_Umap")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  }
  lens=length(rownames(obj@assays$METABOLISM$score))
  if(lens<=3){
  print("few pathways in the dataset detected")

  umap_df=data.frame(pathway=rownames(obj@assays$METABOLISM$score))
  umap_df$UMAP_1=c(0,-1,1)[1:nrow(umap_df)]
  umap_df$UMAP_2=c(0,0,0)[1:nrow(umap_df)]

  # 手动扩展边界
  x_limits <- range(umap_df$UMAP_1)
  y_limits <- range(umap_df$UMAP_2)
  x_range <- diff(x_limits)
  y_range <- diff(y_limits)
  expand_factor <- 0.2
  total_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(alpha = 0.4, size = size, color = "4DBBD5FF",stroke = 1.5) +
    geom_text_repel(aes(label = pathway), size = size/2, color = "black", alpha = 0.6,
                    box.padding = 0.5, point.padding = 0.5,
                    min.segment.length = 0, max.overlaps = Inf) +
    scale_color_viridis(discrete = TRUE) +  # 使用 viridis 调色板
    theme_bw() +
    theme(
      legend.position = "none",  # 移除图例
      panel.grid.major = element_blank(),  # 移除主网格线
      panel.grid.minor = element_blank(),  # 移除次网格线
      panel.border = element_blank(),  # 移除面板边框
      axis.line = element_blank(),  # 移除坐标轴线
      axis.ticks = element_blank(),  # 移除坐标轴刻度
      axis.text = element_blank(),  # 移除坐标轴文字
      axis.title = element_blank()  # 移除坐标轴标题
    )+
    expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                  y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

  print(total_plot)
  # 保存总的 UMAP 图
  result <- list()
  if(!dynamic){
  ggsave(filename = file.path(output_dir, "total_umap_plot.png"), plot = total_plot, width = Width, height = Height)
  write.csv(umap_df, file.path(output_dir, "Pathloc.csv"))
    }
  result[["total_umap_plot"]] <- total_plot
  result[["total_umap_df"]] <- umap_df
  cluster_present=levels(obj@meta.data[,phenotype])

  # 保存分的UMAP图
  if(lens==1){
    sep_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(alpha = 0.4, size = size, color = "grey",stroke = 1.5) +
      scale_color_viridis(discrete = TRUE) +  # 使用 viridis 调色板
      theme_bw() +
      theme(
        legend.position = "none",  # 移除图例
        panel.grid.major = element_blank(),  # 移除主网格线
        panel.grid.minor = element_blank(),  # 移除次网格线
        panel.border = element_blank(),  # 移除面板边框
        axis.line = element_blank(),  # 移除坐标轴线
        axis.ticks = element_blank(),  # 移除坐标轴刻度
        axis.text = element_blank(),  # 移除坐标轴文字
        axis.title = element_blank()  # 移除坐标轴标题
      )+expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                    y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

  for(selelct_cluster in cluster_present){
  if(!dynamic){
  ggsave(filename = file.path(output_dir, paste0(selelct_cluster,".png")), plot = sep_plot, width = Width, height = Height)
  }
  result[[selelct_cluster]] <- sep_plot
    }}else
    { # 获取 cluster 信息
    phenotype = phenotype
    metadata <- obj@meta.data
    metabolism.matrix <- obj@assays$METABOLISM$score
    metadata[, phenotype] <- as.character(metadata[, phenotype])
    metabolism.matrix_sub <- t(metabolism.matrix)

    # 构建数据表
    gg_table <- data.frame(cluster = rep(metadata[, phenotype], each = ncol(metabolism.matrix_sub)),
                           Pathway = rep(colnames(metabolism.matrix_sub), times = nrow(metadata)),
                           Score = as.numeric(metabolism.matrix_sub))
    cluster_pathway_means <- gg_table %>%
      group_by(cluster, Pathway) %>%
      summarise(mean_score = mean(Score, na.rm = TRUE))
    t_test_results_rows <- nrow(cluster_pathway_means)
    cluster_present <- unique(cluster_pathway_means$cluster)
    t_test_results <- cluster_pathway_means %>%
      rowwise() %>%
      mutate(t_score = t_test_t_score(gg_table, cluster, Pathway))
    # 设置阈值和 top_n
    top_pathways <- t_test_results %>%
      filter(t_score > threshold) %>%
      group_by(cluster) %>%
      arrange(desc(t_score)) %>%
      slice_head(n = top_n)
    if(!dynamic){
    write.csv(cluster_pathway_means, file.path(output_dir, "Pathscore.csv"))
    }

    for (selelct_cluster in cluster_present) {
      Pathway_vari <- filter(top_pathways, cluster == selelct_cluster) %>% .$Pathway

      umap_df <- umap_df %>%
        mutate(is_selected = ifelse(pathway %in% Pathway_vari, "selected", "not_selected"))

      center_x <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_1))
      center_y <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_2))

      p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = is_selected), alpha = 0.4, size = size,stroke = 1.5) +
        geom_mark_ellipse(data = umap_df %>% filter(is_selected == "selected"),
                          aes(x = UMAP_1, y = UMAP_2), expand = unit(0.5, "cm"), label.fill = NA) +  # 画圈
        geom_text_repel(data = umap_df %>% filter(is_selected == "selected"),
                        aes(label = pathway), size = 2*size/3, color = "black", alpha = 0.7) +
        annotate("text", x = center_x, y = center_y, label = paste(selelct_cluster),
                 size = size, color = scales::alpha("#4DBBD5FF", 0.7), fontface = "bold") +  # 中心标签
        scale_color_manual(values = c("selected" = "#E64B35FF", "not_selected" = "grey")) +
        theme_bw() +
        theme(
          legend.position = "none",  # 移除图例
          panel.grid.major = element_blank(),  # 移除主网格线
          panel.grid.minor = element_blank(),  # 移除次网格线
          panel.border = element_blank(),  # 移除面板边框
          axis.line = element_blank(),  # 移除坐标轴线
          axis.ticks = element_blank(),  # 移除坐标轴刻度
          axis.text = element_blank(),  # 移除坐标轴文字
          axis.title = element_blank()  # 移除坐标轴标题
        )+expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                      y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

      print(p)
      # 保存每个 cluster 的 UMAP 图
      if(!dynamic){
      ggsave(filename = file.path(output_dir, paste0(selelct_cluster, ".png")), plot = p, width = Width, height = Height)
      }
      result[[selelct_cluster]] <- p
    }
      }
  }


  else{
  # 提取通路得分矩阵
  pathway_scores <- obj@assays$METABOLISM$score %>% t()
  seurat_path <- CreateSeuratObject(counts = pathway_scores)
  #seurat_path <- ScaleData(seurat_path, features = rownames(seurat_path))
  #seurat_path <- RunPCA(seurat_path, features = rownames(seurat_path), npcs =2)  # 尽可能保证npcs小于维度数
  #seurat_path <- RunUMAP(seurat_path, dims = 1:2,n.neighbors=n.neighbors,n.components = 2)


  seurat_path <- RunUMAP(seurat_path, assay = "RNA", features = rownames(seurat_path), n.neighbors = n.neighbors)

  seurat_path@meta.data$Pathway <- colnames(seurat_path)
  seurat_path <- SetIdent(seurat_path, value = 'Pathway')
  umap_coords <- Embeddings(seurat_path, "umap")
  umap_df <- as.data.frame(umap_coords)
  umap_df$pathway <- colnames(seurat_path)

  # 获取 cluster 信息
  phenotype = phenotype
  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  metadata[, phenotype] <- as.character(metadata[, phenotype])
  metabolism.matrix_sub <- t(metabolism.matrix)

  # 构建数据表
  gg_table <- data.frame(cluster = rep(metadata[, phenotype], each = ncol(metabolism.matrix_sub)),
                         Pathway = rep(colnames(metabolism.matrix_sub), times = nrow(metadata)),
                         Score = as.numeric(metabolism.matrix_sub))

  # Initialize progress bar with total steps
  cluster_pathway_means <- gg_table %>%
    group_by(cluster, Pathway) %>%
    summarise(mean_score = median(Score, na.rm = TRUE)) #########重点改变
  t_test_results_rows <- nrow(cluster_pathway_means)
  cluster_present <- unique(cluster_pathway_means$cluster)
  total_steps <- 1 + 1 + t_test_results_rows + length(cluster_present)  # 数据表 + 平均分数 + t-score + UMAP图
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )
  pb$tick()  # Update progress bar after creating data table



  # 计算每个 cluster 和 pathway 的 t-score
  t_test_results <- cluster_pathway_means %>%
    rowwise() %>%
    mutate(t_score = t_test_t_score(gg_table, cluster, Pathway))
  pb$tick(t_test_results_rows)  # Update progress bar for t-score calculations

  # 设置阈值和 top_n
  top_pathways <- t_test_results %>%
    filter(t_score > threshold) %>%
    group_by(cluster) %>%
    arrange(desc(t_score)) %>%
    slice_head(n = top_n)

  # 生成总的 UMAP 图
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  # 手动扩展边界
  x_limits <- range(umap_df$UMAP_1)
  y_limits <- range(umap_df$UMAP_2)
  x_range <- diff(x_limits)
  y_range <- diff(y_limits)
  expand_factor <- 0.2
  # 生成总的 UMAP 图
  total_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(alpha = 0.4, size = size, color = "4DBBD5FF",stroke = 1.5) +
    geom_text_repel(aes(label = pathway), size = size/2, color = "black", alpha = 0.6,
                    box.padding = 0.5, point.padding = 0.5,
                    min.segment.length = 0, max.overlaps = Inf) +
    scale_color_viridis(discrete = TRUE) +  # 使用 viridis 调色板
    theme_bw() +
    theme(
      legend.position = "none",  # 移除图例
      panel.grid.major = element_blank(),  # 移除主网格线
      panel.grid.minor = element_blank(),  # 移除次网格线
      panel.border = element_blank(),  # 移除面板边框
      axis.line = element_blank(),  # 移除坐标轴线
      axis.ticks = element_blank(),  # 移除坐标轴刻度
      axis.text = element_blank(),  # 移除坐标轴文字
      axis.title = element_blank()  # 移除坐标轴标题
    ) +
    expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                  y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

  print(total_plot)
  result=list()
  # 保存总的 UMAP 图
  if(!dynamic){
  ggsave(filename = file.path(output_dir, "total_umap_plot.png"), plot = total_plot, width = Width, height = Height)
  write.csv(cluster_pathway_means, file.path(output_dir, "Pathscore.csv"))
  write.csv(umap_df, file.path(output_dir, "Pathloc.csv"))
  }
  result[["total_umap_plot"]] <- total_plot
  result[["total_umap_df"]] <- umap_df

  # 生成每个 cluster 的 UMAP 图
  for (selelct_cluster in cluster_present) {
    Pathway_vari <- filter(top_pathways, cluster == selelct_cluster) %>% .$Pathway

    umap_df <- umap_df %>%
      mutate(is_selected = ifelse(pathway %in% Pathway_vari, "selected", "not_selected"))

    center_x <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_1))
    center_y <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_2))

    p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = is_selected), alpha = 0.4, size = size,stroke = 1.5) +
      geom_mark_ellipse(data = umap_df %>% filter(is_selected == "selected"),
                        aes(x = UMAP_1, y = UMAP_2), expand = unit(0.5, "cm"), label.fill = NA) +  # 画圈
      geom_text_repel(data = umap_df %>% filter(is_selected == "selected"),
                      aes(label = pathway), size = 2*size/3, color = "black", alpha = 0.7) +
      annotate("text", x = center_x, y = center_y, label = paste(selelct_cluster),
               size = size, color = scales::alpha("#4DBBD5FF", 0.7), fontface = "bold") +  # 中心标签
      scale_color_manual(values = c("selected" = "#E64B35FF", "not_selected" = "grey")) +
      theme_bw() +
      theme(
        legend.position = "none",  # 移除图例
        panel.grid.major = element_blank(),  # 移除主网格线
        panel.grid.minor = element_blank(),  # 移除次网格线
        panel.border = element_blank(),  # 移除面板边框
        axis.line = element_blank(),  # 移除坐标轴线
        axis.ticks = element_blank(),  # 移除坐标轴刻度
        axis.text = element_blank(),  # 移除坐标轴文字
        axis.title = element_blank()  # 移除坐标轴标题
      ) +
      expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                    y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

    print(p)
    # 保存每个 cluster 的 UMAP 图
    if(!dynamic){
    ggsave(filename = file.path(output_dir, paste0(selelct_cluster, ".png")), plot = p, width = Width, height = Height)
    }
    pb$tick()  # Update progress bar after plotting
    result[[selelct_cluster]] <- p
  }
  }
  return(result)
}

