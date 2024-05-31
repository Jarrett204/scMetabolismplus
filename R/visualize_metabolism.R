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

DimPlot.metabolism <- function(obj, pathway, dimention.reduction.type = "umap", dimention.reduction.run = T, size= 1.5){

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
    output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dimplot")
    dir.create(output_dir, showWarnings = FALSE)

    for (input.pathway in pathway) {
      pathway_data <- subset(signature_ggplot2, Pathway == input.pathway)
      plot <- ggplot(data = pathway_data, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
        geom_point(size = size) +
        scale_color_gradientn(colours = pal) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        ggtitle(paste("Pathway:", input.pathway)) +
        theme_bw()
      ggsave(filename = paste0(output_dir, "/", "plot_", input.pathway, ".png"), plot = plot, width = 6, height = 5)
      pb$tick()  # Update progress bar after plotting
    }
    print(plot)  # 打印或保存plot
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
    output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dimplot")
    dir.create(output_dir, showWarnings = FALSE)

    for (input.pathway in pathway) {
      pathway_data <- subset(signature_ggplot2, Pathway == input.pathway)
      plot <- ggplot(data = pathway_data, aes(x = UMAP_1, y = UMAP_2, color = Score)) +
        geom_point(size = size) +
        scale_color_gradientn(colours = pal) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        ggtitle(paste("Pathway:", input.pathway)) +
        theme_bw()
      ggsave(filename = paste0(output_dir, "/", "plot_", input.pathway, ".png"), plot = plot, width = 6, height = 5)
      pb$tick()  # Update progress bar after plotting
    }
  }
  plot
}

DotPlot.metabolism <- function(obj, pathway, phenotype, norm = "y"){
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

  if(length(gg_table_median_norm$X1) == 0){
    cat("No pathway qualified\n\n")
  } else {
    plot_dot <- ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[,1], y = gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
      geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[,3])) +
      ylab("Metabolic Pathway") + xlab(input.parameter) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_line(color = "gray", size = 0.5),
            panel.grid.major = element_line(color = "gray", size = 0.8)) +
      scale_color_gradientn(colours = pal) +
      labs(color = "Value", size = "Value") +
      NULL

    cat("You can paste pathways here for boxplot and so on\n\n")
    print(gg_table_median_norm$X2 %>% unique())
    print(plot_dot)

    # Create output directory and save plot
    output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Dotplot")
    dir.create(output_dir, showWarnings = FALSE)
    ggsave(filename = paste0(output_dir, "/", "plot_Dot", ".pdf"), plot = plot_dot, width = 6, height = 5)

    result <- list(
      plot = plot_dot,
      pathway = gg_table_median_norm$X2 %>% unique()
    )
    return(result)
  }
}

BoxPlot.metabolism <- function(obj, pathway, phenotype, ncol = 1){
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
  for (i in 1:length(input.pathway)){
    gg_table <- rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i]))
    pb$tick()  # Update progress bar after data processing
  }
  gg_table <- data.frame(gg_table)
  colnames(gg_table) <- c("cluster", "Pathway", "Score")
  gg_table$Score <- as.numeric(as.character(gg_table$Score))
  print(head(gg_table))

  # Combine multiple color palettes manually
  colors <- c(pal_jama()(10),pal_aaas()(10),pal_frontiers()(10))
  output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Boxplot")
  dir.create(output_dir, showWarnings = FALSE)

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

    plot_box <- ggplot(data = pathway_data, aes(x = cluster, y = Score, fill = cluster)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.4) +
      ylab("Metabolic Pathway") +
      xlab("Input Parameter") +
      theme_bw() +
      ggtitle(paste("Pathway:", select.pathway)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid.minor = element_line(),
            panel.grid.major = element_line()) +
      labs(fill = "Input Parameter") +
      scale_fill_manual(values = colors)

    ggsave(filename = paste0(output_dir, "/", "plot_", select.pathway, ".pdf"), plot = plot_box, width = 6, height = 5)
    pb$tick()  # Update progress bar after plotting
  }
  plot_box
}

PathUmp.metabolism <- function(obj, threshold = 3, top_n = 5) {
  library(progress)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(ggrepel)
  library(ggforce)

  # 检查并创建输出目录
  output_dir <- paste0("./", unique(obj@meta.data$Cancer), "_", unique(obj@meta.data$dataset), "Path_Umap")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # 提取通路得分矩阵
  pathway_scores <- obj@assays$METABOLISM$score %>% t()
  seurat_path <- CreateSeuratObject(counts = pathway_scores)
  seurat_path <- RunUMAP(seurat_path, assay = "RNA", features = rownames(seurat_path), n.neighbors = 5)
  seurat_path@meta.data$Pathway <- colnames(seurat_path)
  seurat_path <- SetIdent(seurat_path, value = 'Pathway')
  umap_coords <- Embeddings(seurat_path, "umap")
  umap_df <- as.data.frame(umap_coords)
  umap_df$pathway <- colnames(seurat_path)

  # 获取 cluster 信息
  phenotype = "seurat_clusters"
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
    summarise(mean_score = mean(Score, na.rm = TRUE))
  t_test_results_rows <- nrow(cluster_pathway_means)
  cluster_present <- unique(cluster_pathway_means$cluster)
  total_steps <- 1 + 1 + t_test_results_rows + length(cluster_present)  # 数据表 + 平均分数 + t-score + UMAP图
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )
  pb$tick()  # Update progress bar after creating data table

  # 计算每个 cluster 中每个 pathway 的平均分数
  pb$tick()  # Update progress bar after calculating means

  # 定义一个函数来计算 t 检验的 t-score
  t_test_t_score <- function(data, cluster_id, pathway_id) {
    current_cluster_scores <- data %>%
      filter(cluster == cluster_id & Pathway == pathway_id) %>%
      pull(Score)
    other_cluster_scores <- data %>%
      filter(cluster != cluster_id & Pathway == pathway_id) %>%
      pull(Score)
    t_test_result <- t.test(current_cluster_scores, other_cluster_scores)
    return(t_test_result$statistic)
  }

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

  # 手动扩展边界
  x_limits <- range(umap_df$UMAP_1)
  y_limits <- range(umap_df$UMAP_2)
  x_range <- diff(x_limits)
  y_range <- diff(y_limits)
  expand_factor <- 0.2

  # 生成总的 UMAP 图
  total_plot <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = pathway), alpha = 0.6, size = 4) +
    geom_text_repel(aes(label = pathway), size = 3, color = "black", alpha = 0.6) +
    scale_color_lancet() +
    ggtitle("UMAP") +
    theme_bw() +
    theme(legend.position = "none") +
    expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                  y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

  # 保存总的 UMAP 图
  ggsave(filename = file.path(output_dir, "total_umap_plot.pdf"), plot = total_plot, width = 5, height = 5)
  write.csv(cluster_pathway_means, file.path(output_dir, "Pathscore.csv"))

  # 生成每个 cluster 的 UMAP 图
  for (selelct_cluster in cluster_present) {
    Pathway_vari <- filter(top_pathways, cluster == selelct_cluster) %>% .$Pathway

    umap_df <- umap_df %>%
      mutate(is_selected = ifelse(pathway %in% Pathway_vari, "selected", "not_selected"))

    center_x <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_1))
    center_y <- mean(umap_df %>% filter(is_selected == "selected") %>% pull(UMAP_2))

    p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(aes(color = is_selected), alpha = 0.6, size = 4) +
      geom_mark_ellipse(data = umap_df %>% filter(is_selected == "selected"),
                        aes(x = UMAP_1, y = UMAP_2), expand = unit(0.5, "cm"), label.fill = NA) +  # 画圈
      geom_text_repel(data = umap_df %>% filter(is_selected == "selected"),
                      aes(label = pathway), size = 3, color = "black", alpha = 0.7) +
      annotate("text", x = center_x, y = center_y, label = paste("UMAP - Cluster", selelct_cluster),
               size = 6, color = scales::alpha("#4DBBD5FF", 0.35), fontface = "bold") +  # 中心标签
      scale_color_manual(values = c("selected" = "#E64B35FF", "not_selected" = "grey")) +
      ggtitle(paste("UMAP - Cluster", selelct_cluster)) +
      theme_bw() +
      theme(legend.position = "none") +
      expand_limits(x = c(x_limits[1] - x_range * expand_factor, x_limits[2] + x_range * expand_factor),
                    y = c(y_limits[1] - y_range * expand_factor, y_limits[2] + y_range * expand_factor))  # 扩展网格边界

    # 保存每个 cluster 的 UMAP 图
    ggsave(filename = file.path(output_dir, paste0("umap_plot_cluster_", selelct_cluster, ".pdf")), plot = p, width = 5, height = 5)
    pb$tick()  # Update progress bar after plotting
  }
}


