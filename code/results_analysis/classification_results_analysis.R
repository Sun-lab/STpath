rm(list = ls())
library(knitr)
library(stringr)
library(ggplot2)
library(patchwork)
library(car)
library(yarrr)
library(gridExtra)
library(tidyverse)
library(kableExtra)
library(ggbeeswarm)
library(caret)
library(ggplot2)     
library(grid)
library(gridExtra)           
library(likert)
library(multiROC)
library(cvms)
library(RColorBrewer)
theme_set(theme_bw())


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


###### TASK: Wu et al. (CV) ----

# path for accuracy and loss plots
metrics_dir = "/Users/zhiningsui/GitHub/STpath/output/Wu_2021/classification/CV/eval_Wu_CV_11_Classification"
# path for ROC with testing dataset
ROC_dir = "/Users/zhiningsui/GitHub/STpath/output/Wu_2021/classification/CV/roc_Wu_CV_11_Classification"

df_files = paste0(c("training_", "testing_", "validation_"), "Wu_CV_11_Classification", ".csv")
df = NULL
for(f1 in df_files){
  s1 <- gsub(".csv", "",f1)
  dataset <- str_split(s1, pattern = "_")[[1]][1]
  p1 <- read.csv(file.path("/Users/zhiningsui/GitHub/STpath/output/Wu_2021/classification/CV/", f1))[,-1]
  p1$dataset <- dataset
  df <- rbind(df, p1)
}

df$dataset <- factor(df$dataset, levels = c("training", "validation", "testing"))

create_levels <- function(df){
  colnames(df)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                         "dropout_rate", "dense_layer_size")
  # colnames(df)[7:ncol(df)] <- gsub("\\.", " ", colnames(df)[7:ncol(df)]) 
  rownames(df) <- NULL
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128", "256", "512"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.01", "0.001", "0.0001"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))
  
  return(df)
}

#* Section 2 ----

val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

#* -- fig: val_eval ----
p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dense_layer_size)) +
  geom_beeswarm(cex = 4) +
  labs(subtitle = "Validation Loss", 
       x = "Batch size",
       y = "Categorical Crossentropy") +
  facet_wrap(~sample, scales = "free") +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_color_discrete("Dense layer units")

ggplot(data = val_eval,
       aes(x = batch_size, y = val_loss,
           color = dense_layer_size,
           shape = dropout_rate)) +
  geom_beeswarm(cex = 4) +
  labs(subtitle = "Validation Loss", 
       x = "Batch size",
       y = "Categorical Crossentropy") +
  facet_wrap(~learning_rate) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_color_discrete("Dense layer units")

ggplot(data = val_eval,
       aes(x = batch_size, y = accuracy,
           color = dense_layer_size)) +
  geom_beeswarm(cex = 4) +
  labs(subtitle = "Validation Metrics", 
       x = "Batch size", 
       y = "Accuracy") +
  facet_wrap(~sample) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_color_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = accuracy,
                 color = dense_layer_size)) +
  geom_beeswarm(cex = 4) +
  labs(subtitle = "Validation Metrics", 
       x = "Batch size", 
       y = "Accuracy") +
  facet_wrap(~sample, scales = "free") +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_color_discrete("Dense layer units")

mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), 
                              p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("chap09/NN_classification_Wu_per_sample_val_eval.png", p, width = 11, height = 4)

loss_lowest <- val_eval %>%
  filter(val_loss == min(val_loss, na.rm=TRUE))

pred_files <- loss_lowest[,1:7] %>%
  unite(name, sep="_")
pred_files$name <- paste0("test_pred_", pred_files$name, ".csv")

#* Section 3 ----

for(f1 in pred_files$name){
  pred = NULL
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1), check.names = F)
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
  colnames(pred)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                           "dropout_rate", "dense_layer_size")
  pred$label.pred <- factor(pred$label.pred + 1)
  class <- colnames(pred)[10:ncol(pred)]
  pred$Prediction <- class[pred$label.pred]
  
  pred_true <- merge(df, pred, by = "X")
  pred_true <- pred_true[,c("X", "Classification", "Prediction")]
  conf_mat <- confusion_matrix(targets = pred_true$Classification,
                               predictions = pred_true$Prediction)
  
  
  sample <- unique(pred_true$sid)
  
  p <- plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]], 
                             add_sums = TRUE,
                             sums_settings = sum_tile_settings(
                               palette = "Oranges",
                               label = "Total",
                               tc_tile_border_color = "black"
                             ),
                             add_arrows = T,
                             arrow_size = 0.08,
                             arrow_nudge_from_text = 0.1,
                             add_zero_shading = T,
                             # diag_percentages_only = T,
                             rotate_y_text = F,
                             counts_on_top = T,
                             place_x_axis_above = F,
                             font_counts = list("color" = "red", "size" = ((length(class)+1)/2), "nudge_y" = -0.15, "fontface" = "bold"),
                             font_normalized = list("size" = ((length(class)+1)/2), "nudge_y" = 0.25, "fontface" = "bold"),
                             font_row_percentages = list("size" = ((length(class)+1)/3), "alpha" = 0.6),
                             font_col_percentages = list("size" = ((length(class)+1)/3), "alpha" = 0.6),
                             theme_fn = theme_bw) +
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(color = "black", vjust = 0.5, hjust=1),
          plot.title = element_text(face = "bold", size = (2*(length(class)+1))),
          axis.title = element_text(size = (2*(length(class)+1))),
          axis.text = element_text(size = (1.7*(length(class)+1)))) +
    labs(x = "Classification", y = "Prediction",
         title = paste0("Confusion Matrix of ", sample)) +
    scale_x_discrete(
      position = "bottom",
      limits = rev(levels(conf_mat$`Confusion Matrix`$Target))
    )
  
  p
  
  ggsave(paste0(sample, "_cm.png"), p, width = (1.1*(length(class)+1)), height = (1.1*(length(class)+1)))
}


#* -- fig: ROC_test ----

classes = c("Average", class)
myColors <- c("black", "#FF0000", "#008000", "#0000FF")
names(myColors) <- classes
myLinetypes <- c("solid", rep("dashed", 3))
names(myLinetypes) <- classes

pred = NULL
s1 <- gsub(".csv", "",f1)
params_1 <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
p1 <- read.csv(file.path(metrics_dir, f1), check.names = F)
params <- params_1[rep(1, each = nrow(p1)),]
pred <- rbind(pred, cbind(params, p1))
colnames(pred)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                         "dropout_rate", "dense_layer_size")
pred$label.pred <- factor(pred$label.pred)
pred$Prediction <- pred$label.pred

pred_true <- merge(df, pred, by = "X")

roc_curve <- roc(pred_true$Classification, pred_true$Stroma)
plot(roc_curve, col = "blue")
auc_value <- auc(roc_curve)
text(0.3, 0.1, labels = paste("AUC =", round(auc_value, 2)), col = "red")





y_scores <- pred_true[,(ncol(pred_true)-3):(ncol(pred_true)-1)]
colnames(y_scores) <- paste0(colnames(y_scores), "_pred_NN")

y_onehot <- fastDummies::dummy_cols(pred_true[,"Classification"])
colnames(y_onehot) <- c('drop', gsub(".data_", "", colnames(y_onehot)[-1]))
y_onehot <- subset(y_onehot, select = -c(drop))
class = colnames(y_onehot)
colnames(y_onehot) <- paste0(colnames(y_onehot),"_true")
z = cbind(y_onehot, y_scores)

roc <- multi_roc(z)
auc <- as.data.frame(roc$AUC$NN) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column()
colnames(auc) <- c("class", "AUC")
AUC_df <- cbind(params_1, auc)

plot_roc_df <- plot_roc_data(roc)
plot_roc_df <- plot_roc_df[plot_roc_df$Group != "Macro", ] # remove Macro
plot_roc_df$Group <- ifelse(plot_roc_df$Group == "Micro", "Average", plot_roc_df$Group)
plot_roc_df$Legend <- ifelse(plot_roc_df$Group == "Average", paste0(plot_roc_df$Group, ' (AUC = ', round(plot_roc_df$AUC,3), ')'),
                             paste0("Cluster ", plot_roc_df$Group, " vs Others", ' (AUC = ', round(plot_roc_df$AUC,3), ')'))


plot_roc_df$Legend <- factor(plot_roc_df$Legend, 
                             levels = sort(unique(plot_roc_df$Legend))) 
plot_roc_df$Group <- factor(plot_roc_df$Group, levels = c("Average", class))
plot_roc_df$Type <- ifelse(plot_roc_df$Group %in% c("Average"), "1", "0")

group_order <- levels(plot_roc_df$Group)
names(group_order) <- levels(plot_roc_df$Legend)
colScale <- scale_colour_manual(name = "Class (AUC)", 
                                values = myColors,
                                breaks = group_order,
                                labels = names(group_order))

linetype_order <- levels(plot_roc_df$Group)
names(linetype_order) <- levels(plot_roc_df$Legend)

linetypeScale <- scale_linetype_manual(name = "Class (AUC)", 
                                       values = myLinetypes,
                                       breaks = linetype_order,
                                       labels = names(linetype_order))

p_AUC <- ggplot(plot_roc_df[plot_roc_df$Group == "Average", ], 
            aes(x = 1-Specificity, y=Sensitivity)) + 
  geom_path(aes(color = Group, linetype = Group), size = 0.7) + 
  labs(x = "False Positive Rate",
       y = "True Positive Rate",
       title = "") +
  theme_bw() + 
  geom_text("AUC")
  theme(plot.title = element_text(hjust = 0.5, size = 17), 
        legend.justification=c(1, 0), 
        legend.position="none", 
        legend.title=element_blank(), 
        legend.text = element_text(size = 12),
        legend.key.width = unit(1, 'cm'),
        legend.background = element_rect(fill="transparent", 
                                         size=0.5, 
                                         linetype="solid", 
                                         colour ="black"),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 15))  +
  colScale +
  linetypeScale

ggsave(paste0(title, ".png"), p, width = 5, height = 5)


