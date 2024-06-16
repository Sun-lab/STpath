library(knitr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(tidyverse)
library(kableExtra)
library(ggbeeswarm)
library(caret)
library(likert)
library(multiROC)
library(cvms)
library(RColorBrewer)
library(Seurat)
library(gtable)
library(reshape2)
library(cowplot)
library(fastDummies)
theme_set(theme_cowplot())

source("visualization_helper.R")

# Function to create factor levels
create_levels <- function(df, is_wu = TRUE){
  if (is_wu) {
    colnames(df)[1:7] <- c("sample", "image", "optimizer", "batch_size", "learning_rate", 
                           "dropout_rate", "dense_layer_size")
    df$optimizer <- factor(df$optimizer, levels = c("Adam", "SGD"))
    df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
    df$learning_rate <- factor(df$learning_rate, levels = c("0.0001"))
    df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0"))
    df$dense_layer_size <- factor(df$dense_layer_size, levels = c("256", "512"))
  } else {
    colnames(df)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                           "dropout_rate", "dense_layer_size")
    df$batch_size <- factor(df$batch_size, levels = c("32", "64", "128", "256", "512"))
    df$learning_rate <- factor(df$learning_rate, levels = c("0.1", "0.001", "0.0001"))
    df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
    df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))
  }
  rownames(df) <- NULL
  return(df)
}

# Function to load and process datasets
load_and_process_datasets <- function(path, pattern, is_wu = TRUE) {
  files = list.files(path, pattern = pattern)
  data = NULL
  for (f in files) {
    dataset_name <- str_split(gsub(".csv", "", f), pattern = "_")[[1]][2]
    df <- read.csv(file.path(path, f))[, -1]
    df$dataset <- dataset_name
    df <- df[, c("X", "patientid", "dataset", "subtype", "Classification")]
    data <- rbind(data, df)
  }
  data$dataset <- factor(data$dataset, levels = c("training", "validation", "testing"))
  return(data)
}

# Load Wu et al. datasets
data_dir_wu = "../../output/Wu_2021/classification"
eval_dir_wu = "../../output/Wu_2021/classification/Prediction_multiclass_per_sample/"
val_dir_wu = "../../output/Wu_2021/classification/val_multiclass_per_sample/"


df_wu <- load_and_process_datasets(data_dir_wu, pattern = ".csv")

# Evaluation for validation datasets
val_files_wu = list.files(val_dir_wu, pattern = "metrics_")
val_eval_wu = NULL
for (f in val_files_wu) {
  params <- as.data.frame(t(str_split(gsub(".csv", "", f), pattern = "_")[[1]][-1]))
  p <- read.csv(file.path(val_dir_wu, f))
  val_eval_wu <- rbind(val_eval_wu, cbind(params, t(p)))
}

val_eval_wu <- create_levels(val_eval_wu)
colnames(val_eval_wu)[8:9] <- c("loss", "accuracy")

loss_lowest_wu <- val_eval_wu %>%
  group_by(sample) %>%
  filter(loss == min(loss, na.rm = TRUE))



files_acc = list.files(val_dir, pattern = "acc_")
files_loss = list.files(val_dir, pattern = "loss_")



###### TASK: Wu et al. ----
# path for accuracy and loss plots
eval_dir = "/Users/zhiningsui/GitHub/DL_ST_BRCA/output/Wu_2021/Classification/Metrics_multiclass_per_sample/"
# path for the prediction of testing dataset
eval_dir = "/Users/zhiningsui/GitHub/DL_ST_BRCA/output/Wu_2021/Classification/Prediction_multiclass_per_sample/"
# path for evaluation with validation dataset
val_dir = "/Users/zhiningsui/GitHub/DL_ST_BRCA/output/Wu_2021/Classification/val_multiclass_per_sample/"
# path for ROC with testing dataset
ROC_dir = "/Users/zhiningsui/GitHub/DL_ST_BRCA/output/Wu_2021/Classification/ROC_multiclass_per_sample/"

create_levels <- function(df){
  colnames(df)[1:7] <- c("sample", "image", "optimizer", "batch_size", "learning_rate", 
                         "dropout_rate", "dense_layer_size")
  # colnames(df)[7:ncol(df)] <- gsub("\\.", " ", colnames(df)[7:ncol(df)]) 
  rownames(df) <- NULL
  df$optimizer <- factor(df$optimizer, levels = c("Adam", "SGD"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("256", "512"))
  
  return(df)
}

df_files = list.files("../../output/Wu_2021/classification", pattern=".csv")
df = NULL
for(f1 in df_files){
  s1 <- gsub(".csv", "",f1)
  dataset <- str_split(s1, pattern = "_")[[1]][2]
  p1 <- read.csv(file.path("../../output/Wu_2021/classification/", f1))[,-1]
  p1$dataset <- dataset
  p1 <- p1[, c("X", "patientid", "dataset", "subtype", "Classification")]
  df <- rbind(df, p1)
}

df$dataset <- factor(df$dataset, levels = c("training", "validation", "testing"))
table(df$Classification, df$dataset, df$patientid)

#* val_eval ----
val_files = list.files(val_dir, pattern="metrics_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][-1]))
  p1 <- read.csv(file.path(val_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, t(p1)))
}

val_eval <- create_levels(val_eval)
colnames(val_eval)[8:9] <- c("loss", "accuracy")

loss_lowest <- val_eval %>%
  group_by(sample) %>%
  filter(loss == min(loss, na.rm=TRUE))

#* -- tab: loss_lowest ----
kbl(loss_lowest[,c(1, 8:9, 4,7)], digits = 3, format = "latex",
    caption = "Batch size and hidden dense layer units of the neutral network that gave the lowest validation loss at the optimal epoch of each sample.",
    col.names = c("Sample ID", "Validation Loss", "Validation Accuracy",
                  "Batch size", "Dense layer Units"),
    escape = FALSE) 

pred_files <- loss_lowest[,1:7] %>%
  unite(name, sep="_")
pred_files$name <- paste0("pred_", pred_files$name, ".csv")

#* -- fig: ROC_test ----
test_files = list.files("../output/Wu_2021/Classification/", pattern="_testing.csv")
test_df = NULL
for(f1 in test_files){
  s1 <- gsub(".csv", "",f1)
  dataset <- str_split(s1, pattern = "_")[[1]][2]
  p1 <- read.csv(file.path("../output/Wu_2021/Classification/", f1))[,-1]
  p1$dataset <- dataset
  p1 <- p1[, c("X", "patientid", "dataset", "subtype", "Classification")]
  test_df <- rbind(test_df, p1)
}

classes = c("Macro-average", "Micro-average", unique(test_df$Classification))
myColors <- c("black", "black", c("#FF0000", "#008000", "#0000FF", "#00CC00", "#800080", "#ff7f00",
                                  "#1f78b4", "#FF00FF", "#8b8b8b", "#000080", "#1f78b4", "#8c564b"))
names(myColors) <- classes
myLinetypes <- c("dashed", "solid", rep("solid", length(unique(test_df$Classification))))
names(myLinetypes) <- classes

AUC_df <- data.frame()
for(f1 in pred_files$name){
  pred = NULL
  s1 <- gsub(".csv", "",f1)
  params_1 <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][-1]))
  p1 <- read.csv(file.path(eval_dir, f1), check.names = F)
  params <- params_1[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
  colnames(pred)[1:7] <- c("sample", "image", "optimizer", "batch_size", "learning_rate", 
                           "dropout_rate", "dense_layer_size")
  pred$label.pred <- factor(pred$label.pred + 1)
  class <- colnames(pred)[10:ncol(pred)]
  pred$Prediction <- class[pred$label.pred]
  
  pred_true <- merge(df, pred, by = "X")
  y_scores <- pred_true[,14:(ncol(pred_true)-1)]
  colnames(y_scores) <- paste0(colnames(y_scores), "_pred_NN")
  
  y_onehot <- fastDummies::dummy_cols(pred_true$Classification)
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
  AUC_df <- rbind(AUC_df, cbind(params_1, auc))
}

colnames(AUC_df)[1:7] <- c("sample", "image", "optimizer", "batch_size", "learning_rate", 
                           "dropout_rate", "dense_layer_size")

AUC_df$class <- gsub("\\.\\.\\.", " + ", AUC_df$class)
AUC_df$class <- gsub("\\.", " ", AUC_df$class)
AUC_df$class <- gsub("micro", "Average", AUC_df$class)
AUC_df$class <- gsub("macro", "Macro-average", AUC_df$class)

AUC_df$class <- factor(AUC_df$class, 
                       levels = c("Macro-average", "Average", "Stroma", 
                                  "Invasive cancer + stroma + lymphocytes", 
                                  "Invasive cancer + lymphocytes", 
                                  "Invasive cancer + stroma", "Lymphocytes",
                                  "Invasive cancer", "DCIS",
                                  "Normal + stroma + lymphocytes",
                                  "Normal glands + lymphocytes",
                                  "Stroma + adipose tissue",
                                  "Adipose tissue","Necrosis"))

p_AUC_Wu <- ggplot(AUC_df[!AUC_df$class %in% c("Macro-average"),], 
                          aes(x = class, y = AUC, color = sample)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.box.margin = margin(0,0,0,-250)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  labs(color = "Sample") + 
  scale_x_discrete(limits=rev) + 
  coord_flip() 

###### TASK: He et al. ----
fns <- c("He_integrated.cca.lognorm_clusters_0.1")
# path for accuracy and loss plots
eval_dir = paste0("/Users/zhiningsui/GitHub/STpath/output/He_2020/classification/eval_He_integrated.cca.lognorm_clusters_0.1/") 
# path for the prediction of testing dataset
eval_dir = paste0("../output/BRCA/classification/pred_", fns) 
# path for evaluation with validation dataset
val_dir = paste0("../output/BRCA/classification/val_", fns) 
# path for ROC with testing dataset
ROC_dir = paste0("../output/BRCA/classification/ROC_", fns) 

df_files = paste0(c("training_", "testing_", "validation_"), rep(fns[1],3), ".csv")
df = NULL
for(f1 in df_files){
  s1 <- gsub(".csv", "",f1)
  dataset <- str_split(s1, pattern = "_")[[1]][1]
  p1 <- read.csv(file.path("/Users/zhiningsui/GitHub/STpath/output/He_2020/classification/", f1))[,-1]
  p1$dataset <- dataset
  df <- rbind(df, p1)
}

df$dataset <- factor(df$dataset, levels = c("training", "validation", "testing"))

create_levels <- function(df){
  colnames(df)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                         "dropout_rate", "dense_layer_size")
  # colnames(df)[7:ncol(df)] <- gsub("\\.", " ", colnames(df)[7:ncol(df)]) 
  rownames(df) <- NULL
  df$batch_size <- factor(df$batch_size, levels = c("32", "64", "128", "256", "512"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.1", "0.001", "0.0001"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))
  
  return(df)
}

#* Section 2 ----

val_files = list.files(eval_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(eval_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

#* val_eval ----
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


loss_lowest <- val_eval %>%
  filter(val_loss == min(val_loss, na.rm=TRUE))

pred_files <- loss_lowest[,1:7] %>%
  unite(name, sep="_")
f1 <- paste0("test_pred_", pred_files$name, ".csv")
f1 <- "test_pred_ResNet50_o_Adam_256_0.001_0.0_256.csv"
#* -- fig: ROC_test ----
classes = c("Average", "0", "1", "2")
myColors <- c("black", "#FF0000", "#008000", "#0000FF")
names(myColors) <- classes
myLinetypes <- c("solid", rep("dashed", 3))
names(myLinetypes) <- classes

pred = NULL
s1 <- gsub(".csv", "",f1)
params_1 <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
p1 <- read.csv(file.path(eval_dir, f1), check.names = F)
params <- params_1[rep(1, each = nrow(p1)),]
pred <- rbind(pred, cbind(params, p1))
colnames(pred)[1:7] <- c("basemodel", "image", "optimizer", "batch_size", "learning_rate", 
                         "dropout_rate", "dense_layer_size")
pred$label.pred <- factor(pred$label.pred)
pred$Prediction <- pred$label.pred

pred_true <- merge(df, pred, by = "X")
y_scores <- pred_true[,(ncol(pred_true)-3):(ncol(pred_true)-1)]
colnames(y_scores) <- paste0(colnames(y_scores), "_pred_NN")

y_onehot <- fastDummies::dummy_cols(pred_true[,"integrated.cca.lognorm_clusters_0.1"])
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
colScale <- scale_colour_manual(name = NULL, 
                                values = myColors,
                                breaks = group_order,
                                labels = names(group_order))

linetype_order <- levels(plot_roc_df$Group)
names(linetype_order) <- levels(plot_roc_df$Legend)

linetypeScale <- scale_linetype_manual(name = NULL, 
                                       values = myLinetypes,
                                       breaks = linetype_order,
                                       labels = names(linetype_order))

p_AUC_He <- ggplot(plot_roc_df, 
                   aes(x = 1-Specificity, y=Sensitivity)) + 
  geom_path(aes(color = Group, linetype = Group), size = 0.7) + 
  labs(x = "False Positive Rate",
       y = "True Positive Rate",
       title = "") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 17), 
        legend.justification=c(1, 0), 
        legend.position=c(.98, .02), 
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



# ========== FILTERED (zeros all samples, 90%) ========== ----------------------

### directories ----
data_dir = "../../data/He_2020/"
cluster_dir = "../../output/He_2020/clustering/"

annotations <- read.csv(file.path(data_dir, "He_clusters_filt_zeros_all_90.csv")) %>%
  remove_rownames() %>%
  column_to_rownames("X") %>%
  mutate(across(where(is.integer), as.factor))

tab1 <- as.data.frame.matrix(table(annotations$pid, annotations$label))
tab2 <- as.data.frame.matrix(table(annotations$pid, annotations$integrated.cca.lognorm_clusters_0.1))

tab <- bind_cols(tab1, tab2)

library(kableExtra)
kbl(tab, col.names = c("Non", "Tumor", "Cluster 0", "Cluster 1", "Cluster 2"),
    format = "latex") %>%
  kable_classic() %>%
  add_header_above(c(" " = 1, "He et al." = 2, "Seurat v5" = 3))

seurat_clustered <- readRDS(file.path(cluster_dir, "clustered_st_obj_filt_zeros_all.rds"))
p_list <- visualize_umap(cluster_methods = "integrated.cca.lognorm_clusters_0.1",
                         seurat_obj = seurat_clustered,
                         fn_suffix = "integrated.cca.lognorm_clusters_0.1",
                         split.by.1 = "label",
                         n_col.1 = 1,
                         strip.text.size = 16,
                         axis.text.size = 12,
                         legend.text.size = 12,
                         axis.title.size = 14,
                         legend.point.size = 4,
                         tag.size = 14,
                         pdf_height = 6,
                         pdf_width = 7, 
                         plot_titles = "",
                         includes_number = T,
                         returns_plot = 1)

p_umap <- p_list[[1]]

p_AUC_Wu <- p_AUC_Wu + 
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16)) 

fig6 <- free(p_AUC_Wu) / ((p_umap | free(p_AUC_He)) & 
                            theme(axis.text = element_text(size = 14),
                                  axis.title = element_text(size = 16),
                                  legend.text = element_text(size = 14),
                                  legend.title = element_text(size = 16))) + 
  plot_annotation(tag_levels = 'A',
                  tag_prefix = "(",
                  tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 16, face = "bold")) 

ggsave("../../figure/paper_figure/Fig6.classification.jpg", fig6, width = 12, height = 11)
















fns <- c("He_integrated.cca.lognorm_clusters_0.1")
# path for accuracy and loss plots
eval_dir = paste0("/Users/zhiningsui/GitHub/STpath/output/Wu_2021/classification/eval_Wu_multiclass_per_sample/") 
# path for the prediction of testing dataset
eval_dir = paste0("../output/BRCA/classification/pred_", fns) 
# path for evaluation with validation dataset
val_dir = paste0("../output/BRCA/classification/val_", fns) 
# path for ROC with testing dataset
ROC_dir = paste0("../output/BRCA/classification/ROC_", fns) 











