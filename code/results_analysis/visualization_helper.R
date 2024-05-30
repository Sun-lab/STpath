# Load libraries

library(ggplot2)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(gtable)
library(grid)
library(patchwork)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

visualize_spots <- function(mapping_df, 
                            label, 
                            x_loc, y_loc, 
                            title, 
                            strata = NULL, 
                            colors = NULL, 
                            pointSize){
  
  param_list <- list(label, x_loc, y_loc, title, strata)
  chars <- which(sapply(param_list, is.character))
  
  for(i in chars){
    indx <- which(colnames(mapping_df) == param_list[[i]])
    if(is_empty(indx)){
      param_list[[i]] = ncol(mapping_df)+1
    }else{
      param_list[[i]] = indx
    }
  }
  
  if(sum(sapply(param_list, function(indx) indx > ncol(mapping_df))) > 0){
    stop("At least one column specified is not in the mapping dataframe.")
  }
  
  if(is.null(colors)){
    library(scales)
    # colors = RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(mapping_df[,label]))]
    # colors = grDevices::rainbow(length(unique(mapping_df[,label])), alpha = 0.7)
    colors <- scales::hue_pal()(length(unique(mapping_df[,label])))
    colors <- setNames(colors, sort(unique(mapping_df[,label])))
  }else{
    colors = colors
  }
  
  label_plot <- mapping_df[, as.numeric(param_list)]
  colnames(label_plot) <- c("label", "x", "y", "title", "strata")
  
  label_plot$x <- as.numeric(label_plot$x)
  label_plot$x <- label_plot$x - min(label_plot$x)
  label_plot$y <- as.numeric(label_plot$y)
  label_plot$y <- max(label_plot$y) - label_plot$y
  
  p = suppressMessages(ggplot(label_plot, aes(x = x, y = y)) +
                         geom_point(aes(colour = label),size = pointSize) +
                         scale_color_manual(values = colors) +
                         xlim(range_x) +
                         ylim(range_y) +
                         # scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
                         # scale_y_reverse() +
                         # scale_x_reverse() +
                         coord_fixed()+
                         facet_wrap(~strata) +
                         ggtitle(unique(label_plot[,"title"])) +
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}

visualize_spots_continous <- function(mapping_df, value, x_loc, y_loc, title, colors = NULL, NumCols = 4, pointSize = 1, rescale = F){
  
  param_list <- c(value, list(x_loc, y_loc, title))
  chars <- which(sapply(param_list, is.character))
  
  for(i in chars){
    indx <- which(colnames(mapping_df) == param_list[[i]])
    if(is_empty(indx)){
      param_list[[i]] = ncol(mapping_df)+1
    }else{
      param_list[[i]] = indx
    }
  }
  
  # if(sum(sapply(param_list, function(indx) indx > ncol(mapping_df))) > 0){
  #   stop("At least one column specified is not in the mapping dataframe.")
  # }
  # 
  if(is.null(colors)){
    colors = c("lightblue","lightyellow","red")
  }else{
    colors = colors
  }
  
  plot_data <- mapping_df[, as.numeric(param_list)]
  
  plot_value = plot_data[, value]
  
  if(rescale){
    if(!is.null(ncol(plot_value))){
      plot_value_scale = as.data.frame(apply(plot_value,2,function(x){
        (x - min(x)) / (max(x) - min(x))
      } ))
    }else{
      plot_value_scale = as.data.frame((plot_value - min(plot_value)) / (max(plot_value) - min(plot_value)))
      colnames(plot_value_scale) = value
    }
  }else{
    plot_value_scale = as.data.frame(plot_value) 
    colnames(plot_value_scale) = value
  }
  
  plot_value_scale$x <- as.numeric(plot_data[,x_loc])
  plot_value_scale$x <- plot_value_scale$x - min(plot_value_scale$x)
  plot_value_scale$y <- as.numeric(plot_data[,y_loc])
  plot_value_scale$y <- max(plot_value_scale$y) - plot_value_scale$y
  
  mData = melt(plot_value_scale,id.vars = c("x","y"))
  colnames(mData)[3] <- "feature"
  
  
  p = suppressMessages(ggplot(mData, aes(x, y)) + 
                         geom_point(aes(colour = value),size = pointSize) +
                         scale_color_gradientn(colours = colors) + 
                         xlim(range_x) +
                         ylim(range_y) +
                         # scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
                         facet_wrap(~feature, ncol = NumCols)+ 
                         ggtitle(unique(plot_data[,title])) +
                         coord_fixed()+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}

visualize_umap <- function(cluster_methods,
                           seurat_obj,
                           fn_prefix = "umap_integrated",
                           fn_suffix,
                           pdf_height = 6,
                           pdf_width = 7, 
                           split.by.1 = NULL, 
                           split.by.2 = NULL,
                           n_col.1 = NULL,
                           n_col.2 = NULL,
                           alpha = 0.4,
                           strip.text.size = 12,
                           axis.text.size = 8,
                           legend.text.size = 9,
                           axis.title.size = 12,
                           legend.point.alpha = 1,
                           legend.point.size = 2,
                           tag.size = 13,
                           plot_titles = NULL,
                           includes_number = F,
                           includes_tag = F){
  
  pdf(file.path(cluster_dir, paste0(fn_prefix, "_", fn_suffix, ".pdf")), 
      width = pdf_width, height = pdf_height, onefile=FALSE)
  
  num = 1
  for (cm in cluster_methods) {
    m <- strsplit(cm, "_")[[1]][1]
    norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
    
    if(norm == "lognorm"){
      DefaultAssay(seurat_obj) <- "RNA"
    }else if(norm == "sct"){
      DefaultAssay(seurat_obj) <- "SCT"
    }else{
      break
    }
    
    Idents(seurat_obj) <- cm
    colors_cluster <- scales::hue_pal()(length(unique(seurat_obj@meta.data[, cm])))
    colors_cluster <- setNames(colors_cluster, levels(seurat_obj@meta.data[,cm]))
            
    p1 <- DimPlot(seurat_obj,
                  reduction = paste0("umap.", m),
                  label = F,
                  cols = colors_cluster,
                  alpha = alpha,
                  split.by = split.by.1,
                  ncol = n_col.1) +
      ggtitle(cm) +
      theme(strip.text = element_text(size = strip.text.size),
            axis.text = element_text(size = axis.text.size),
            legend.text = element_text(size = legend.text.size),
            axis.title = element_text(size = axis.title.size),
            legend.position = 'right',
            plot.title = element_text(hjust = 0.5)) +
      guides(color=guide_legend(nrow=ceiling(length(colors_cluster)/12),
                                byrow=TRUE,
                                override.aes = list(size = legend.point.size,
                                                    alpha = legend.point.alpha))) +
      labs(x = "umap_1", y = "umap_2", color = "Cluster")
    
    if(!is.null(plot_titles)){
      p1 <- p1 + ggtitle(plot_titles[num])
    }
    
    if(includes_number){
      counts_1 <- seurat_obj@meta.data %>% 
        group_by(get(cm)) %>% 
        summarise(count = n()) %>%
        as.data.frame()
      colnames(counts_1)[1] <- cm
      
      p1 <- p1 + geom_text(data = counts_1, 
                           aes(x = c(-2, 3, 8), y = Inf, 
                               label = paste0("Count: ", count),
                               color = get(cm)),
                           vjust = 2,
                           alpha = 1,
                           size = 5, 
                           inherit.aes = FALSE,
                           show.legend = FALSE)
        
    }
    
    if(!is.null(split.by.2)){
      p2 <- DimPlot(seurat_obj,
                    reduction = paste0("umap.", m),
                    label = F,
                    split.by = split.by.2,
                    cols = colors_cluster,
                    alpha = alpha,
                    ncol = n_col.2) +
        theme(strip.text = element_text(size = strip.text.size),
              axis.text = element_text(size = axis.text.size),
              legend.text = element_text(size = legend.text.size),
              axis.title = element_text(size = axis.title.size),
              legend.position = 'right') +
        guides(color=guide_legend(nrow=ceiling(length(colors_cluster)/12),
                                  byrow=TRUE,
                                  override.aes = list(size = legend.point.size,
                                                      alpha = legend.point.alpha))) +
        labs(x = "umap_1", y = "umap_2", color = "Cluster")
 
      if(includes_number){
        counts_2 <- seurat_obj@meta.data %>% 
          group_by(get(cm), label) %>% 
          summarise(count = n()) %>%
          as.data.frame()
        colnames(counts_2)[1] <- cm
        counts_2$x <- ifelse(counts_2[, cm] == 0, -2,
                             ifelse(counts_2[, cm] == 1, 3, 8))
        counts_2$y <- Inf 
          
        p2 <- p2 + geom_text(data = counts_2, 
                             aes(x = x, y = y, 
                                 label = paste0("Count: ", count),
                                 color = get(cm)),
                             alpha = 1,
                             vjust = 2,
                             size = 4, 
                             inherit.aes = FALSE,
                             show.legend = FALSE)
        
      }
      p <- (p1/p2) + plot_layout(ncol = 1) +
        plot_layout(guides = 'collect') & 
        theme(legend.position = 'bottom',
              legend.justification = "center")
      
      if(includes_tag){
        p <- p + plot_annotation(tag_levels = 'A',
                               tag_prefix = "(",
                               tag_suffix = ")") & 
          theme(plot.tag = element_text(size = tag.size))
      }
      print(p)
    
    }else{
      print(p1)
    }
    num = num + 1
  }
  dev.off()
}

visualize_markers_heatmap <- function(cluster_methods, 
                                      seurat_obj, 
                                      markers_list, 
                                      fn_prefix = "heatmap_top10_DE_markers", 
                                      fn_suffix = "filt_zeros_all",
                                      pdf_width = 20, 
                                      pdf_height = 15){
  
  pdf(file.path(cluster_dir, paste0(fn_prefix, "_", fn_suffix, ".pdf")), 
      width = pdf_width, height = pdf_height)
  for (cm in cluster_methods) {
    print(cm)
    m <- strsplit(cm, "_")[[1]][1]
    norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
    
    Idents(object = seurat_obj) <- cm
    
    top10_cm <- markers_list[[cm]]
    if(nrow(top10_cm)==0) next
    
    if(norm == "lognorm"){
      DefaultAssay(seurat_obj) <- "RNA"
      scale_genes <- rownames(seurat_obj@assays$RNA$scale.data)
    }else if(norm == "sct"){
      DefaultAssay(seurat_obj) <- "SCT"
      scale_genes <- rownames(seurat_obj@assays$SCT$scale.data)
    }else{
      break
    }
    
    if(all(!top10_cm$gene %in% scale_genes)){
      plot.new()
      text(x = 0.35, y = 0.5, adj = c(0,1), 
           labels = paste0(cm, ": No genes in scale.data"))
      next
    }
    
    top10_cm <- top10_cm %>%
      arrange(cluster, rev(avg_log2FC))
    top10_cm$label <- paste0("(Cluster ", top10_cm$cluster, ", ", "avg LFC = ", 
                             round(top10_cm$avg_log2FC,2), ")")
    top10_cm$duplicated <- 0
    
    for(i in 1:nrow(top10_cm)){
      gene <- top10_cm[i,"gene"]
      if(gene %in% top10_cm$gene[c(i+1: nrow(top10_cm))]){
        loc = which(top10_cm$gene %in% gene)[-1]
        duplicated <- top10_cm[top10_cm$gene %in% gene,]
        top10_cm[i, "label"] <- paste0(duplicated$label, collapse = "\n")
        top10_cm[i, "duplicated"] <- 1
        top10_cm[loc, "duplicated"] <- 2
      }
    }
    
    top10_cm$label <- paste(top10_cm$gene_symbol,  top10_cm$label)
    top10_cm[top10_cm$duplicated == 0,]$label <- paste0(top10_cm[top10_cm$duplicated == 0,]$label, "\n") 
    top10_cm <- top10_cm[top10_cm$duplicated != 2,]
    
    label <- top10_cm$label
    names(label) <- top10_cm$gene
    
    top10_cm$caption <- paste0(top10_cm$gene_symbol, " (Cluster ", top10_cm$cluster, ")")
    caption_cm <- paste("Not in scale data: ", 
                        paste(top10_cm[!top10_cm$gene %in% scale_genes, "caption"], 
                              collapse = ", "))
    
    p <- DoHeatmap(seurat_obj, 
                   features = top10_cm$gene, 
                   combine = TRUE) &
      ggtitle(cm) &
      scale_y_discrete(labels = label) &
      plot_annotation(caption = caption_cm)
    
    print(p)
    rm(p)
  }
  dev.off()
}

visualize_markers_vlnplot <- function(cluster_methods, 
                                      seurat_obj, 
                                      markers_list, 
                                      annotations_df, 
                                      n_row = 3, 
                                      n_col = 10, 
                                      fn_prefix = "vlnplot_top10_DE_markers", 
                                      fn_suffix = "filt_zeros_all",
                                      pdf_width = 40, 
                                      pdf_height = 30, 
                                      exclude_zeros = T){
  
  for (cm in cluster_methods) {
    print(cm)
    colors_cluster <- scales::hue_pal()(length(unique(seurat_obj@meta.data[,cm])))
    # colors_cluster <- adjustcolor(colors_cluster, alpha.f = 0.5)
    colors_cluster <- setNames(colors_cluster, levels(seurat_obj@meta.data[,cm]))
    
    m <- strsplit(cm, "_")[[1]][1]
    norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
    
    Idents(object = seurat_obj) <- cm
    
    top10_cm <- markers_list[[cm]]
    if(nrow(top10_cm)==0) next
    
    if(norm == "lognorm"){
      DefaultAssay(seurat_obj) <- "RNA"
      scale_genes <- rownames(seurat_obj@assays$RNA$scale.data)
    }else if(norm == "sct"){
      DefaultAssay(seurat_obj) <- "SCT"
      scale_genes <- rownames(seurat_obj@assays$SCT$scale.data)
    }else{
      break
    }
    
    top10_cm <- top10_cm %>%
      arrange(cluster, rev(avg_log2FC))
    top10_cm$label <- paste0("(Cluster ", top10_cm$cluster, ", ", "avg LFC = ", round(top10_cm$avg_log2FC,2), ")")
    top10_cm$duplicated <- 0
    
    for(i in 1:nrow(top10_cm)){
      gene <- top10_cm[i,"gene"]
      if(gene %in% top10_cm$gene[c(i+1: nrow(top10_cm))]){
        loc = which(top10_cm$gene %in% gene)[-1]
        duplicated <- top10_cm[top10_cm$gene %in% gene,]
        top10_cm[i, "label"] <- paste0(duplicated$label, collapse = "\n")
        top10_cm[i, "duplicated"] <- 1
        top10_cm[loc, "duplicated"] <- 2
      }
    }
    
    top10_cm$label <- paste(top10_cm$gene_symbol, top10_cm$label)
    top10_cm[top10_cm$duplicated == 0,]$label <- paste0(top10_cm[top10_cm$duplicated == 0,]$label, "\n") 
    top10_cm <- top10_cm[top10_cm$duplicated != 2,]
    
    label <- top10_cm$label
    names(label) <- top10_cm$gene
    
    top10_cm$caption <- paste0(top10_cm$gene_symbol, " (Cluster ", top10_cm$cluster, ")")
    caption_cm <- paste("Not in scale data: ", paste(top10_cm[!top10_cm$gene %in% scale_genes, "caption"], 
                                                     collapse = ", "))
    
    mtx <- seurat_obj[["RNA"]]$counts
    mtx <- mtx[top10_cm$gene,]
    df_matrix <- as.data.frame(mtx)
    df_matrix$gene <- row.names(df_matrix)
    df_long <- pivot_longer(df_matrix, 
                            cols = -gene, 
                            names_to = "X", 
                            values_to = "Value")
    
    annotations_cm <- annotations_df[colnames(df_matrix), c("X", cm)]
    df_joined <- left_join(df_long, annotations_cm, by = "X")
    df_joined$`Not Expressed` <- df_joined$Value == 0
    
    proportion_zeros <- df_joined %>%
      group_by(gene, .dots = cm) %>%
      count(`Not Expressed`) %>%
      mutate(zero_prop = n/sum(n))
    
    if(exclude_zeros){
      data <- FetchData(seurat_obj,
                        vars = top10_cm$gene,
                        layer = "data")
      
      df <- as.data.frame(data)
      df$Spot <- rownames(df)
      long_df <- pivot_longer(df, 
                              cols = -Spot, 
                              names_to = "Gene", 
                              values_to = "Value")
      long_df[long_df == 0] <- NA
      long_df <- cbind(long_df, seurat_obj[[cm]][long_df$Spot, ])
      colnames(long_df)[4] <- "cluster"
      
      p_list_no_zeros <- list()
      for(gene in top10_cm$gene){
        p <- ggplot(long_df[long_df$Gene==gene,],
                    aes(x = cluster, y = Value, fill = cluster)) +
          geom_violin() +
          geom_jitter(size = 0.1, size = 1)+
          scale_fill_manual(values = colors_cluster) +
          labs(y = "Expression Level", x = "Identity", fill = NULL, title = gene) 
        p_list_no_zeros[[gene]] <- p
      }
    }
    
    p_list <- VlnPlot(seurat_obj, 
                      features = top10_cm$gene, 
                      combine = F,
                      cols = colors_cluster)
    
    
    
    for (i in 1:length(p_list)) {
      p <- p_list[[i]]
      title <- p$labels$title
      new_title <- gsub(" \\(", "\n\\(", label[names(label) == title])
      
      p_zeros <- ggplot(proportion_zeros[proportion_zeros$gene == title, ],
                        aes_string(x = cm, y = "zero_prop",  fill = "`Not Expressed`")) + 
        geom_col(position = position_stack(), width = 0.8) + 
        geom_text(aes(label= paste0(n, " (", round(zero_prop*100,1), "%)")), 
                  position = position_stack(vjust = 0.5), size = (2.5 / (length(colors_cluster) / 2))) +
        labs(title = new_title, fill = NULL, x = "Cluster Identity", y = "Proportion")+
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
      
      # p_zeros <- ggplot(proportion_zeros[proportion_zeros$gene == title, ],
      #                   aes_string(x = cm, y = "zero_prop",  fill = cm)) +
      #   geom_bar(stat = "identity", color = "black", width = 0.8) +
      #   scale_fill_manual(values = colors_cluster) +
      #   labs(title = new_title, fill = NULL, x = "Identity")+
      #   theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) 
      # 
      
      if(exclude_zeros){
        p_new <- p + ggtitle("Including Zeros") + labs(x = "Cluster Identity")
        p_new_no_zeros <- p_list_no_zeros[[i]] + ggtitle("Excluding Zeros") 
        p_list[[i]] <- p_zeros/p_new/p_new_no_zeros +
          plot_layout(heights = c(1, 2, 2))
      }else{
        p_new <- p + ggtitle(NULL) + labs(x = "Cluster Identity")
        p_list[[i]] <- p_zeros/p_new +
          plot_layout(heights = c(1, 2))
      }
    }
    
    pdf(file.path(cluster_dir, paste0(fn_prefix, "_", cm, "_", fn_suffix, ".pdf")), 
        width = pdf_width, height = pdf_height)
    
    n <- n_col*n_row
    page <- length(p_list) %/% n
    remainder <- length(p_list) %% n
    
    if(page > 0){
      for (i in 1:page) {
        p <- wrap_plots(p_list[(n*(i-1)+1):(n*i)]) + 
          plot_layout(ncol=n_col, nrow = n_row) +
          plot_annotation(cm, theme=theme(plot.title=element_text(hjust=0.5, size = 22, face = "bold")))
        
        print(p)
      }
    }
    
    if(remainder > 0){
      p <- wrap_plots(p_list[(n*page+1):length(p_list)]) + 
        plot_layout(ncol=n_col, nrow = n_row) +
        plot_annotation(cm, theme=theme(plot.title=element_text(hjust=0.5, size = 22, face = "bold")))
      print(p)
    }
    
    dev.off()
    
    rm(p)
    gc()
  }
}


visualize_markers_umap <- function(cluster_methods, 
                                   seurat_obj, 
                                   markers_list, 
                                   annotations_df, 
                                   n_row = 5, 
                                   n_col = 10, 
                                   fn_prefix = "umap_top10_DE_markers", 
                                   fn_suffix = "filt_zeros_all",
                                   pdf_width = 35, 
                                   pdf_height = 20){
  for (cm in cluster_methods){
    print(cm)
    
    m <- strsplit(cm, "_")[[1]][1]
    norm <- strsplit(m, ".", fixed = TRUE)[[1]][3]
    
    Idents(object = seurat_obj) <- cm
    
    top10_cm <- markers_list[[cm]]
    if(nrow(top10_cm)==0) next
    
    if(norm == "lognorm"){
      DefaultAssay(seurat_obj) <- "RNA"
      scale_genes <- rownames(seurat_obj@assays$RNA$scale.data)
    }else if(norm == "sct"){
      DefaultAssay(seurat_obj) <- "SCT"
      scale_genes <- rownames(seurat_obj@assays$SCT$scale.data)
    }else{
      break
    }
    
    top10_cm <- top10_cm %>%
      arrange(cluster, rev(avg_log2FC))
    top10_cm$label <- paste0("(Cluster ", top10_cm$cluster, ", ", "avg LFC = ", round(top10_cm$avg_log2FC,2), ")")
    top10_cm$duplicated <- 0
    
    for(i in 1:nrow(top10_cm)){
      gene <- top10_cm[i,"gene"]
      if(gene %in% top10_cm$gene[c(i+1: nrow(top10_cm))]){
        loc = which(top10_cm$gene %in% gene)[-1]
        duplicated <- top10_cm[top10_cm$gene %in% gene,]
        top10_cm[i, "label"] <- paste0(duplicated$label, collapse = "\n")
        top10_cm[i, "duplicated"] <- 1
        top10_cm[loc, "duplicated"] <- 2
      }
    }
    
    top10_cm$label <- paste(top10_cm$gene_symbol, top10_cm$label)
    top10_cm[top10_cm$duplicated == 0,]$label <- paste0(top10_cm[top10_cm$duplicated == 0,]$label, "\n") 
    top10_cm <- top10_cm[top10_cm$duplicated != 2,]
    
    label <- top10_cm$label
    names(label) <- top10_cm$gene
    
    top10_cm$caption <- paste0(top10_cm$gene_symbol, " (Cluster ", top10_cm$cluster, ")")
    caption_cm <- paste("Not in scale data: ", paste(top10_cm[!top10_cm$gene %in% scale_genes, "caption"], 
                                                     collapse = ", "))
    
    p_list <- FeaturePlot(seurat_obj, 
                          reduction = paste0("umap.", strsplit(cm, "_")[[1]][1]),
                          features = top10_cm$gene, combine = F)
    
    for (i in 1:length(p_list)) {
      p <- p_list[[i]]
      title <- p$labels$title
      new_title <- gsub(" \\(", "\n\\(", label[names(label) == title])
      p_new <- p + ggtitle(new_title) +
        theme(plot.title = element_text(size = 12))
      p_list[[i]] <- p_new
    }
    
    pdf(file.path(cluster_dir, paste0(fn_prefix, "_", cm, "_", fn_suffix, ".pdf")), 
        width = pdf_width, height = pdf_height)
    
    n <- n_col*n_row
    page <- length(p_list) %/% n
    remainder <- length(p_list) %% n
    
    if(page > 0){
      for (i in 1:page) {
        p <- wrap_plots(p_list[(n*(i-1)+1):(n*i)]) + 
          plot_layout(ncol=n_col, nrow = n_row) +
          plot_annotation(cm, theme=theme(plot.title=element_text(hjust=0.5, size = 22, face = "bold")))
        
        print(p)
      }
    }
    
    if(remainder > 0){
      p <- wrap_plots(p_list[(n*page+1):length(p_list)]) + 
        plot_layout(ncol=n_col, nrow = n_row) +
        plot_annotation(cm, theme=theme(plot.title=element_text(hjust=0.5, size = 22, face = "bold")))
      print(p)
    }
    
    dev.off()
    
    rm(p)
    gc()
  }
}


visualize_clusters_comparison_heatmap <- function(annotation_df, 
                                                  colnameList.1,
                                                  colnameList.2,
                                                  name.1 = "Unfiltered", 
                                                  name.2 = "Filtered",
                                                  fn_prefix = "clusters_compare",
                                                  fn_suffix,
                                                  pdf_width = 8,
                                                  pdf_height = 8){
  if(length(colnameList.1) != length(colnameList.2)){
    break
  }
  
  pdf(file.path(cluster_dir, paste0(fn_prefix, "_", fn_suffix, ".pdf")), 
      width = pdf_width,
      height = pdf_height)
  
  
  for (i in 1:length(colnameList.1)) {
    
    setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, 
                                                             height=0.9, 
                                                             name="vp", 
                                                             just=c("right","top"))), 
            action="prepend")
    pheatmap::pheatmap(as.matrix(table(annotation_df[, colnameList.1[i]], 
                                       annotation_df[, colnameList.2[i]])),
                       display_numbers = TRUE,
                       number_format = "%.f",
                       color = colorRampPalette(RColorBrewer::brewer.pal(n = 5, 
                                                                         name ="Reds"))(100),
                       cluster_row = FALSE,
                       cluster_col = FALSE,
                       fontsize = 14,
                       main = paste0(colnameList.1[i], "\nvs\n",colnameList.2[i]))
    setHook("grid.newpage", NULL, "replace")
    grid.text(name.1, x=-0.07, rot=90, gp=gpar(fontsize = 16))
    grid.text(name.2, y=-0.07, gp=gpar(fontsize = 16))
  }
  
  dev.off()
}

visualize_clusters_comparison_spotplot <- function(annotation_df_list,
                                                   colnameList,
                                                   locList = list(c("row", "col"),
                                                                  c("row", "col")),
                                                   titles,
                                                   pdf_width,
                                                   pdf_height,
                                                   fn_prefix = "spotplot",
                                                   fn_suffix,
                                                   fn_names,
                                                   color_list=NULL){
  
  for(i in 1:length(fn_names)){
    labels = unlist(lapply(colnameList, "[", i))
    
    # Assign colors to labels
    colors_cluster_list <- list()
    for(j in 1:length(labels)){
      colors_cluster <- color_list[[j]]
      if(is.null(colors_cluster)){
        values <- unique(annotation_df_list[[j]][, labels[j]])
        colors_cluster <- scales::hue_pal()(length(values))
        colors_cluster <- setNames(colors_cluster, levels(values))
        colors_cluster_list[[j]] <- colors_cluster
      }else{
        colors_cluster_list[[j]] <- colors_cluster
      }
    }
    
    # Create plots
    p_list <- list()
    for (sid in sids) {
      p_list_s <- list()
      for (j in 1:length(labels)) {
        annotation_df <- annotation_df_list[[j]]
        annotation_df_s <- annotation_df[annotation_df$sid == sid, ]
        p <- visualize_spots(mapping_df = annotation_df_s,
                             label = labels[j], x_loc = locList[[j]][1], 
                             y_loc = locList[[j]][2], title = "sid", 
                             strata = "sid", 
                             pointSize = 2, 
                             colors = colors_cluster_list[[j]]) + 
          ggtitle(titles[j]) +
          theme(plot.title = element_text(size = 15, face = "bold"))
        
        p_list_s[[j]] <- p
      }
      
      p_list[[sid]] <- wrap_plots(p_list_s, ncol = 1)
    }
    
    
    pdf(file.path(cluster_dir, paste0(fn_prefix, "_", fn_names[i], "_", fn_suffix, ".pdf")),
        width = pdf_width, height = pdf_height)
    for (pid in pids) {
      p <- wrap_plots(p_list[startsWith(names(p_list), pid)]) +
        plot_layout(nrow = 1) +
        plot_annotation(pid, theme=theme(plot.title=element_text(hjust=0.5, size = 22, face = "bold")))
      print(p)
    }
    dev.off()
  }
}



map_gene_symbol <- function(x){ # x is a vector of gene ids
  for(i in 1:length(x)){
    if(x[i] %in% matched_symbols$id){
      symbols <- na.omit(unique(unlist(matched_symbols[matched_symbols$id==x[i], 3:6])))
      symbol <- ifelse(length(symbols)>1, paste(symbols, collapse = "/"), symbols)
    }else{
      symbol = x[i]
    }
    x[i] <- symbol
  }
  x
}

