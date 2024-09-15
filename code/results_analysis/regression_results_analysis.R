rm(list = ls())
library(knitr)
library(stringr)
library(ggplot2)
library(car)
library(yarrr)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggbeeswarm)
library(kableExtra)
theme_set(theme_bw())

set.seed(123)
colors = c(`B cells` = "#2a5fbd",
           CAFs = "#7ddffa",
           `Cancer Epithelial` = "#117d30",
           Endothelial = "#71f25a",
           Myeloid = "#ebc857",
           `Normal Epithelial` ="#D39200",
           Plasmablasts = "#F8766D",
           PVL = "#DB72FB",
           `T cells` = "#bd2a84")

colors_patho <- c(`Invasive Cancer` = "#C77CFF",
                  Lymphocyte = "#00BFC4",
                  Stroma = "#7CAE00",
                  Others = "#F8766D")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


###### TASK: evaluate_NN_regression_Wu_CV_1142243F ----
create_levels <- function(df){
  colnames(df)[1:7] <- c("base", "image", "optimizer", "batch_size", "learning_rate",
                         "dropout_rate", "dense_layer_size")
  rownames(df) <- NULL
  df$image <- factor(df$image, levels = c("o"))
  df$optimizer <- factor(df$optimizer, levels = c("Adam"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001", "0.001", "0.1"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))

  return(df)
}


# path for accuracy and loss plots
metrics_dir = "../../output/Wu_2021/regression/eval_Wu_CV_1142243F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial/"
# path for cell type proportions obtained by CARD
dvn_dir = "../../output/Wu_2021/regression/"

# load the Observed proportions
major_train <- read.csv("../../output/Wu_2021/regression/training_Wu_CV_1142243F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
major_val <- read.csv("../../output/Wu_2021/regression/validation_Wu_CV_1142243F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
major_test <- read.csv("../../output/Wu_2021/regression/testing_Wu_CV_1142243F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
colnames(major_train) <- gsub("\\.", " ", colnames(major_train))
colnames(major_val) <- gsub("\\.", " ", colnames(major_val))
colnames(major_test) <- gsub("\\.", " ", colnames(major_test))

celltypes <- colnames(major_train)[2:10]
celltypes.obs <- paste(celltypes, "obs", sep = ".")
celltypes.pred <- paste(celltypes, "pred", sep = ".")

major_train_obs <- major_train[, c("X", celltypes)]
major_train_long <- major_train_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Training data")

major_val_obs <- major_val[, c("X", celltypes)]
major_val_long <- major_val_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Validation data")

major_test_obs <- major_test[, c("X", celltypes)]
major_test_long <- major_test_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Testing data")

major_long <- bind_rows(major_train_long, major_val_long, major_test_long)

p <- ggplot(major_long,
       aes(x=factor(name), y=value, color=factor(name))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  scale_color_manual("Cell type", values = colors) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~factor(dataset,
                     levels = c("Training data", "Validation data", "Testing data")),
             scales = "free_x")

ggsave("NN_regression_Wu_CV_1142243F_composition.png", p, width = 6.5, height = 2.5)

val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Loss",
       x = "Batch size",
       y = "MSE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_mae,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Metrics",
       x = "Batch size",
       y = "MAE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")


mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                         p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("NN_regression_Wu_CV_1142243F_val_eval.png", p, width = 13, height = 5)

pred_files = list.files(metrics_dir, pattern="test_pred_")
pred = NULL
for(f1 in pred_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
}
pred <- create_levels(pred)
colnames(pred) <- gsub("\\.", " ", colnames(pred))


colnames(major_test_obs) <- c("X", celltypes.obs)

pred_true <- merge(pred, major_test_obs, by= "X")
pred_true <- pred_true %>%
  arrange(`Cancer Epithelial.obs`)

pred_true_long <- reshape(pred_true, direction =  "long",
                          varying = list(celltypes, celltypes.obs),
                          timevar = "celltype", times = celltypes,
                          v.names = c("proportion.pred", "proportion.obs"),
                          idvar = colnames(pred_true)[1:8])
rownames(pred_true_long) <- NULL

X_order = unlist(unique(pred_true_long[pred_true_long$celltype == "Cancer Epithelial","X"]))
pred_true_long$X <- factor(pred_true_long$X, levels=X_order)

slope <- pred_true_long %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(slope = map_dbl(model, ~summary(.)$coefficients[2,1]),
         model = NULL) %>%
  as.data.frame(.)

ggplot(data = slope,
       aes(x = batch_size, y = slope,
           color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3, alpha = 0.9, size = 1) +
  labs(x = "Batch size of training",
       y = "Slope of OLS",
       shape="Dense layer units", colour="Dropout proportion") +
  facet_grid(celltype~optimizer+learning_rate,
             labeller = labeller(learning_rate = c("0.0001" = "Learning rate = 0.0001",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.1" = "Learning rate = 0.01"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9))
ggsave("NN_regression_CV_1142243F_slope1.png", height = 12, width = 10)

slope_reduced <- slope[!((slope$optimizer == "Adam" & slope$learning_rate == "0.01")),]

ggplot(data = slope_reduced,
       aes(x = batch_size, y = slope,
           color = celltype)) +
  geom_point(size = 1.2) +
  labs(x = "Batch size",
       y = "Slope of OLS",
       shape="Learning rate", colour="Cell type") +
  facet_grid(optimizer+learning_rate~dropout_rate+dense_layer_size,
             labeller = labeller(learning_rate = c("0.01" = "Learning rate = 0.01",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.0001" = "Learning rate = 0.0001"),
                                 dropout_rate = c("0.0" = "No dropout layer",
                                                  "0.2" = "Dropout proportion = 0.2",
                                                  "0.5" = "Dropout proportion = 0.5"),
                                 dense_layer_size = c("0" = "No dense layer",
                                                      "256" = "Dense layer units = 256",
                                                      "512" = "Dense layer units = 512"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9)) +
  scale_color_manual("Cell type", values = colors) +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_CV_1142243F_slope2.png", height = 5, width = 11)

max_slope <- slope %>%
  group_by(celltype) %>%
  dplyr::filter(slope == max(slope, na.rm=TRUE)) %>%
  arrange(slope)

kbl(max_slope[,c(6:7, 1:5)], format = "latex",
    caption = "Parameters of the neutral network that gave the highest slope of best-fit line for each cell type in the training dataset.",
    col.names = c("Cell type", "Slope", "Optimizer", "Batch size",
                  "Learning rate", "Dropout Proportion", "Dense layer Units"),
    escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

tmp <- subset(pred_true_long,
              optimizer == "Adam" &
                batch_size == "64" & learning_rate == "0.0001" &
                dropout_rate == "0.0" & dense_layer_size == "256")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)

r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(celltype) %>%
  summarise(max_x = max(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
max$max_x <- 0.7*max$max_x
max$delta_y <- 0.1*max$range_y
r2_tmp <- merge(r2_tmp, max, by = "celltype")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = transparent("black", trans.val = .7), size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1.5, stroke=NA) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = max_x, y = (min_y+delta_y)),
            size = 5,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 64, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 256",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=23, face = "bold"),
        strip.text.x = element_text(size=17),
        legend.position="none", legend.box = "horizontal",
        legend.key.size = unit(1, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=16), axis.title.x = element_text(size=17),
        axis.title.y= element_text(size=19),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1142243F_scatterplot.png", height = 12, width = 12)


tmp_2 <- reshape(
  tmp,
  direction =  "long",
  varying = c("proportion.pred", "proportion.obs"),
  timevar = "type",
  times = c("pred", "obs"),
  v.names = c("proportion"),
  idvar = colnames(tmp)[1:9]
)

tmp_2 <- tmp_2 %>%
  group_by(celltype) %>%
  arrange(proportion, .by_group = TRUE)

X_order = unlist(unique(tmp_2[tmp_2$type == "obs" & tmp_2$celltype == "Cancer Epithelial","X"]))
tmp_2$X <- factor(tmp_2$X, levels=X_order)
tmp_2$celltype <- factor(tmp_2$celltype, levels = sort(celltypes))
tmp_2$celltype <- relevel(tmp_2$celltype, "Cancer Epithelial")

ggplot(tmp_2, aes(fill=factor(celltype, levels=rev(levels(celltype))),
                  color=factor(celltype, levels=rev(levels(celltype))),
                  y=proportion, x=X)) +
  geom_bar(position="stack", stat="identity", width = 1) +
  facet_wrap(~type,
             labeller = labeller(type = c("obs" = "Observed",
                                          "pred" = "Predicted"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  scale_fill_manual("Cell type", values = colors) +
  scale_color_manual("Cell type", values = colors) +
  labs(x = "Patches in Testing Data", y = "Proportion")
ggsave("NN_regression_Wu_CV_1142243F_obs_vs_pred_composition.png", height = 3.5, width = 10)

p_list = list()
for (i in celltypes) {
  tmp <- tmp_2[tmp_2$celltype ==i, ]

  ordered <- tmp %>%
    group_by(type) %>%
    arrange(proportion, .by_group = TRUE)

  X_order = unique(unlist(ordered$X))
  tmp$X <- factor(tmp$X, levels=X_order)

  p <- ggplot(tmp, aes(fill=factor(celltype), color=factor(celltype), y=proportion, x=X)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(strip.text.x = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18)) +
    scale_fill_manual("Cell type",
                      values=colors) +
    scale_color_manual("Cell type",
                       values=colors) +
    labs(title = i,
         x = "Patches in Testing Data",
         y = "Proportion") +
    facet_wrap(~type,
               labeller = labeller(type = c("obs" = "Observed",
                                            "pred" = "Predicted")))
  p_list[[i]] <- p
}

p <- grid.arrange(p_list[[1]],p_list[[2]],p_list[[3]],
                  p_list[[4]],p_list[[5]],p_list[[6]],
                  p_list[[7]],p_list[[8]],p_list[[9]],nrow = 3)
ggsave("NN_regression_Wu_CV_1142243F_obs_vs_pred_composition_per_type.png", p, height = 6, width = 15)

mu <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), mean)
mu$measure <- "mean"
med <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), median)
med$measure <- "median"

measure <- rbind(mu, med)

ggplot(tmp_2, aes(x = proportion, color = type)) +
  geom_histogram(fill = "white", alpha = 0.05, position = "dodge", bins = 70) +
  geom_vline(data=measure, aes(xintercept=proportion, color=type, linetype = measure))+
  scale_colour_discrete("Proportion", breaks=c("obs", "pred"),
                        labels=c("Observed (from CARD)",
                                 "Predicted (from trained network)")) +
  scale_linetype_discrete("Measure", breaks=c("mean", "median"),
                          labels = c("Mean", "Median")) +
  labs(title="Optimizer = Adam, Batch size = 32, Learning rate = 0.0001,
       Dropout layer proportion = 0.2, Dense layer units = 512",
       x="Proportion",
       y = "Count") +
  theme(legend.position="right", legend.box = "vertical",
        legend.key.size = unit(1, 'cm'), legend.title = element_text(size=19),
        plot.title = element_text(size=22, face = "bold"),
        strip.text.x = element_text(size=20),
        legend.text = element_text(size=17),
        axis.title.x = element_text(size=21),
        axis.title.y= element_text(size=21),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free") +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1142243F_obs_vs_pred_histogram.png", height = 11, width = 15)

major_patho <- data.frame(X = major_test_obs$X,
                          `Invasive.Cancer.obs` = major_test_obs$`Cancer Epithelial.obs`,
                          Stroma.obs = major_test_obs$Endothelial.obs + major_test_obs$PVL.obs + major_test_obs$CAFs.obs,
                          Lymphocyte.obs = major_test_obs$`B cells.obs` + major_test_obs$`T cells.obs` + major_test_obs$Plasmablasts.obs,
                          Others.obs = major_test_obs$`Normal Epithelial.obs` + major_test_obs$Myeloid.obs)

pred_patho <- data.frame(`Invasive.Cancer.pred` = pred$`Cancer Epithelial`,
                         Stroma.pred = pred$Endothelial + pred$PVL + pred$CAFs,
                         Lymphocyte.pred =  pred$`T cells` + pred$`B cells` + pred$Plasmablasts,
                         Others.pred = pred$`Normal Epithelial` + pred$Myeloid)

pred_patho <- cbind(pred[,c(1:8)], pred_patho)
pred_true_patho <- merge(pred_patho, major_patho, by= "X")

pred_true_patho_long <- reshape(pred_true_patho, direction =  "long",
                                varying = list(c("Invasive.Cancer.obs",
                                                 "Stroma.obs",
                                                 "Lymphocyte.obs",
                                                 "Others.obs"),
                                               c("Invasive.Cancer.pred",
                                                 "Stroma.pred",
                                                 "Lymphocyte.pred",
                                                 "Others.pred")),
                                timevar = "pathology",
                                times = c("Invasive Cancer", "Stroma", "Lymphocyte", "Others"),
                                v.names = c("proportion.obs", "proportion.pred"),
                                idvar = colnames(pred_true_patho)[1:8])
rownames(pred_true_patho_long) <- NULL

tmp <-  subset(pred_true_patho_long,
               optimizer == "Adam" &
                 batch_size == "64" & learning_rate == "0.0001" &
                 dropout_rate == "0.0" & dense_layer_size == "256")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, pathology) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(pathology) %>%
  summarise(range_x = max(proportion.obs)-min(proportion.obs),
            min_x = min(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
r2_tmp <- merge(r2_tmp, max, by = "pathology")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 64, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 256",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~pathology, nrow = 1, scales = "free") +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1142243F_scatterplot_patho.png", height = 3, width = 10)

tmp_nine <- tmp




###### TASK: evaluate_NN_regression_Wu_CV_1160920F ----
create_levels <- function(df){
  colnames(df)[1:7] <- c("base", "image", "optimizer", "batch_size", "learning_rate",
                         "dropout_rate", "dense_layer_size")
  rownames(df) <- NULL
  df$image <- factor(df$image, levels = c("o"))
  df$optimizer <- factor(df$optimizer, levels = c("Adam"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001", "0.001", "0.1"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))

  return(df)
}


# path for accuracy and loss plots
metrics_dir = "../../output/Wu_2021/regression/eval_Wu_CV_1160920F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial/"
# path for cell type proportions obtained by CARD
dvn_dir = "../../output/Wu_2021/regression/"

# load the Observed proportions
major_train <- read.csv("../../output/Wu_2021/regression/training_Wu_CV_1160920F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
major_val <- read.csv("../../output/Wu_2021/regression/validation_Wu_CV_1160920F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
major_test <- read.csv("../../output/Wu_2021/regression/testing_Wu_CV_1160920F_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv")[,-1]
colnames(major_train) <- gsub("\\.", " ", colnames(major_train))
colnames(major_val) <- gsub("\\.", " ", colnames(major_val))
colnames(major_test) <- gsub("\\.", " ", colnames(major_test))

celltypes <- colnames(major_train)[2:10]
celltypes.obs <- paste(celltypes, "obs", sep = ".")
celltypes.pred <- paste(celltypes, "pred", sep = ".")

major_train_obs <- major_train[, c("X", celltypes)]
major_train_long <- major_train_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Training data")

major_val_obs <- major_val[, c("X", celltypes)]
major_val_long <- major_val_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Validation data")

major_test_obs <- major_test[, c("X", celltypes)]
major_test_long <- major_test_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Testing data")

major_long <- bind_rows(major_train_long, major_val_long, major_test_long)

p <- ggplot(major_long,
            aes(x=factor(name), y=value, color=factor(name))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  scale_color_manual("Cell type", values = colors) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~factor(dataset,
                     levels = c("Training data", "Validation data", "Testing data")),
             scales = "free_x")

ggsave("NN_regression_Wu_CV_1160920F_composition.png", p, width = 6.5, height = 2.5)

val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Loss",
       x = "Batch size",
       y = "MSE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_mae,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Metrics",
       x = "Batch size",
       y = "MAE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")


mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("NN_regression_Wu_CV_1160920F_val_eval.png", p, width = 13, height = 5)

pred_files = list.files(metrics_dir, pattern="test_pred_")
pred = NULL
for(f1 in pred_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
}
pred <- create_levels(pred)
colnames(pred) <- gsub("\\.", " ", colnames(pred))


colnames(major_test_obs) <- c("X", celltypes.obs)

pred_true <- merge(pred, major_test_obs, by= "X")
pred_true <- pred_true %>%
  arrange(`Cancer Epithelial.obs`)

pred_true_long <- reshape(pred_true, direction =  "long",
                          varying = list(celltypes, celltypes.obs),
                          timevar = "celltype", times = celltypes,
                          v.names = c("proportion.pred", "proportion.obs"),
                          idvar = colnames(pred_true)[1:8])
rownames(pred_true_long) <- NULL

X_order = unlist(unique(pred_true_long[pred_true_long$celltype == "Cancer Epithelial","X"]))
pred_true_long$X <- factor(pred_true_long$X, levels=X_order)

slope <- pred_true_long %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(slope = map_dbl(model, ~summary(.)$coefficients[2,1]),
         model = NULL) %>%
  as.data.frame(.)

ggplot(data = slope,
       aes(x = batch_size, y = slope,
           color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3, alpha = 0.9, size = 1) +
  labs(x = "Batch size of training",
       y = "Slope of OLS",
       shape="Dense layer units", colour="Dropout proportion") +
  facet_grid(celltype~optimizer+learning_rate,
             labeller = labeller(learning_rate = c("0.0001" = "Learning rate = 0.0001",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.1" = "Learning rate = 0.01"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9))
ggsave("NN_regression_Wu_CV_1160920F_slope1.png", height = 12, width = 10)

slope_reduced <- slope[!((slope$optimizer == "Adam" & slope$learning_rate == "0.01")),]

ggplot(data = slope_reduced,
       aes(x = batch_size, y = slope,
           color = celltype)) +
  geom_point(size = 1.2) +
  labs(x = "Batch size",
       y = "Slope of OLS",
       shape="Learning rate", colour="Cell type") +
  facet_grid(optimizer+learning_rate~dropout_rate+dense_layer_size,
             labeller = labeller(learning_rate = c("0.01" = "Learning rate = 0.01",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.0001" = "Learning rate = 0.0001"),
                                 dropout_rate = c("0.0" = "No dropout layer",
                                                  "0.2" = "Dropout proportion = 0.2",
                                                  "0.5" = "Dropout proportion = 0.5"),
                                 dense_layer_size = c("0" = "No dense layer",
                                                      "256" = "Dense layer units = 256",
                                                      "512" = "Dense layer units = 512"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9)) +
  scale_color_manual("Cell type", values = colors) +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_slope2.png", height = 5, width = 11)

max_slope <- slope %>%
  group_by(celltype) %>%
  dplyr::filter(slope == max(slope, na.rm=TRUE)) %>%
  arrange(slope)

kbl(max_slope[,c(6:7, 1:5)], format = "latex",
    caption = "Parameters of the neutral network that gave the highest slope of best-fit line for each cell type in the training dataset.",
    col.names = c("Cell type", "Slope", "Optimizer", "Batch size",
                  "Learning rate", "Dropout Proportion", "Dense layer Units"),
    escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

tmp <- subset(pred_true_long,
              optimizer == "Adam" &
                batch_size == "16" & learning_rate == "0.0001" &
                dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)

r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(celltype) %>%
  summarise(max_x = max(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
max$max_x <- 0.7*max$max_x
max$delta_y <- 0.1*max$range_y
r2_tmp <- merge(r2_tmp, max, by = "celltype")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = transparent("black", trans.val = .7), size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1.5, stroke=NA) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = max_x, y = (min_y+delta_y)),
            size = 5,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 512",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=23, face = "bold"),
        strip.text.x = element_text(size=17),
        legend.position="none", legend.box = "horizontal",
        legend.key.size = unit(1, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=16), axis.title.x = element_text(size=17),
        axis.title.y= element_text(size=19),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1160920F_scatterplot.png", height = 12, width = 12)


tmp_2 <- reshape(
  tmp,
  direction =  "long",
  varying = c("proportion.pred", "proportion.obs"),
  timevar = "type",
  times = c("pred", "obs"),
  v.names = c("proportion"),
  idvar = colnames(tmp)[1:9]
)

tmp_2 <- tmp_2 %>%
  group_by(celltype) %>%
  arrange(proportion, .by_group = TRUE)

X_order = unlist(unique(tmp_2[tmp_2$type == "obs" & tmp_2$celltype == "Cancer Epithelial","X"]))
tmp_2$X <- factor(tmp_2$X, levels=X_order)
tmp_2$celltype <- factor(tmp_2$celltype, levels = sort(celltypes))
tmp_2$celltype <- relevel(tmp_2$celltype, "Cancer Epithelial")

ggplot(tmp_2, aes(fill=factor(celltype, levels=rev(levels(celltype))),
                  color=factor(celltype, levels=rev(levels(celltype))),
                  y=proportion, x=X)) +
  geom_bar(position="stack", stat="identity", width = 1) +
  facet_wrap(~type,
             labeller = labeller(type = c("obs" = "Observed",
                                          "pred" = "Predicted"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  scale_fill_manual("Cell type", values = colors) +
  scale_color_manual("Cell type", values = colors) +
  labs(x = "Patches in Testing Data", y = "Proportion")
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_composition.png", height = 3.5, width = 10)

p_list = list()
for (i in celltypes) {
  tmp <- tmp_2[tmp_2$celltype ==i, ]

  ordered <- tmp %>%
    group_by(type) %>%
    arrange(proportion, .by_group = TRUE)

  X_order = unique(unlist(ordered$X))
  tmp$X <- factor(tmp$X, levels=X_order)

  p <- ggplot(tmp, aes(fill=factor(celltype), color=factor(celltype), y=proportion, x=X)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(strip.text.x = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18)) +
    scale_fill_manual("Cell type",
                      values=colors) +
    scale_color_manual("Cell type",
                       values=colors) +
    labs(title = i,
         x = "Patches in Testing Data",
         y = "Proportion") +
    facet_wrap(~type,
               labeller = labeller(type = c("obs" = "Observed",
                                            "pred" = "Predicted")))
  p_list[[i]] <- p
}

p <- grid.arrange(p_list[[1]],p_list[[2]],p_list[[3]],
                  p_list[[4]],p_list[[5]],p_list[[6]],
                  p_list[[7]],p_list[[8]],p_list[[9]],nrow = 3)
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_composition_per_type.png", p, height = 6, width = 15)

mu <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), mean)
mu$measure <- "mean"
med <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), median)
med$measure <- "median"

measure <- rbind(mu, med)

ggplot(tmp_2, aes(x = proportion, color = type)) +
  geom_histogram(fill = "white", alpha = 0.05, position = "dodge", bins = 70) +
  geom_vline(data=measure, aes(xintercept=proportion, color=type, linetype = measure))+
  scale_colour_discrete("Proportion", breaks=c("obs", "pred"),
                        labels=c("Observed (from CARD)",
                                 "Predicted (from trained network)")) +
  scale_linetype_discrete("Measure", breaks=c("mean", "median"),
                          labels = c("Mean", "Median")) +
  labs(title="Optimizer = Adam, Batch size = 32, Learning rate = 0.0001,
       Dropout layer proportion = 0.2, Dense layer units = 512",
       x="Proportion",
       y = "Count") +
  theme(legend.position="right", legend.box = "vertical",
        legend.key.size = unit(1, 'cm'), legend.title = element_text(size=19),
        plot.title = element_text(size=22, face = "bold"),
        strip.text.x = element_text(size=20),
        legend.text = element_text(size=17),
        axis.title.x = element_text(size=21),
        axis.title.y= element_text(size=21),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free") +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_histogram.png", height = 11, width = 15)

major_patho <- data.frame(X = major_test_obs$X,
                          `Invasive.Cancer.obs` = major_test_obs$`Cancer Epithelial.obs`,
                          Stroma.obs = major_test_obs$Endothelial.obs + major_test_obs$PVL.obs + major_test_obs$CAFs.obs,
                          Lymphocyte.obs = major_test_obs$`B cells.obs` + major_test_obs$`T cells.obs` + major_test_obs$Plasmablasts.obs,
                          Others.obs = major_test_obs$`Normal Epithelial.obs` + major_test_obs$Myeloid.obs)

pred_patho <- data.frame(`Invasive.Cancer.pred` = pred$`Cancer Epithelial`,
                         Stroma.pred = pred$Endothelial + pred$PVL + pred$CAFs,
                         Lymphocyte.pred =  pred$`T cells` + pred$`B cells` + pred$Plasmablasts,
                         Others.pred = pred$`Normal Epithelial` + pred$Myeloid)

pred_patho <- cbind(pred[,c(1:8)], pred_patho)
pred_true_patho <- merge(pred_patho, major_patho, by= "X")

pred_true_patho_long <- reshape(pred_true_patho, direction =  "long",
                                varying = list(c("Invasive.Cancer.obs",
                                                 "Stroma.obs",
                                                 "Lymphocyte.obs",
                                                 "Others.obs"),
                                               c("Invasive.Cancer.pred",
                                                 "Stroma.pred",
                                                 "Lymphocyte.pred",
                                                 "Others.pred")),
                                timevar = "pathology",
                                times = c("Invasive Cancer", "Stroma", "Lymphocyte", "Others"),
                                v.names = c("proportion.obs", "proportion.pred"),
                                idvar = colnames(pred_true_patho)[1:8])
rownames(pred_true_patho_long) <- NULL

tmp <-  subset(pred_true_patho_long,
               optimizer == "Adam" &
                 batch_size == "16" & learning_rate == "0.0001" &
                 dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, pathology) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(pathology) %>%
  summarise(range_x = max(proportion.obs)-min(proportion.obs),
            min_x = min(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
r2_tmp <- merge(r2_tmp, max, by = "pathology")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 512",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~pathology, nrow = 1, scales = "free") +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_scatterplot_patho.png", height = 3, width = 10)

tmp_nine <- tmp


###### TASK: evaluate_NN_regression_CV_Wu_10x ----
create_levels <- function(df){
  colnames(df)[1:7] <- c("base", "image", "optimizer", "batch_size", "learning_rate",
                         "dropout_rate", "dense_layer_size")
  rownames(df) <- NULL
  df$image <- factor(df$image, levels = c("o"))
  df$optimizer <- factor(df$optimizer, levels = c("Adam"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001", "0.001", "0.1"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))

  return(df)
}


# path for accuracy and loss plots
metrics_dir = "../../output/CV/regression/eval_CV_Wu_10x_invasive.cancer_stroma_lymphocyte_others/"
# path for cell type proportions obtained by CARD
dvn_dir = "../../output/CV/regression/"

# load the Observed proportions
major_train <- read.csv("../../output/CV/regression/training_CV_Wu_10x_invasive.cancer_stroma_lymphocyte_others.csv")[,-1]
major_val <- read.csv("../../output/CV/regression/validation_CV_Wu_10x_invasive.cancer_stroma_lymphocyte_others.csv")[,-1]
major_test <- read.csv("../../output/CV/regression/testing_CV_Wu_10x_invasive.cancer_stroma_lymphocyte_others.csv")[,-1]
colnames(major_train) <- gsub("\\.", " ", colnames(major_train))
colnames(major_val) <- gsub("\\.", " ", colnames(major_val))
colnames(major_test) <- gsub("\\.", " ", colnames(major_test))

celltypes <- colnames(major_train)[3:6]
celltypes.obs <- paste(celltypes, "obs", sep = ".")
celltypes.pred <- paste(celltypes, "pred", sep = ".")

major_train_obs <- major_train[, c("X", celltypes)]
major_train_long <- major_train_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Training data")

major_val_obs <- major_val[, c("X", celltypes)]
major_val_long <- major_val_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Validation data")

major_test_obs <- major_test[, c("X", "sid", celltypes)]
major_test_long <- major_test_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Testing data")

major_long <- bind_rows(major_train_long, major_val_long, major_test_long)

p <- ggplot(major_long,
            aes(x=factor(name), y=value, color=factor(name))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  scale_color_manual("Cell type", values = colors) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~factor(dataset,
                     levels = c("Training data", "Validation data", "Testing data")),
             scales = "free_x")

ggsave("NN_regression_CV_Wu_10x_composition.png", p, width = 6.5, height = 2.5)

val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Loss",
       x = "Batch size",
       y = "MSE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_mae,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Metrics",
       x = "Batch size",
       y = "MAE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")


mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("NN_regression_Wu_CV_1160920F_val_eval.png", p, width = 13, height = 5)

pred_files = list.files(metrics_dir, pattern="test_pred_")
pred = NULL
for(f1 in pred_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
}
pred <- create_levels(pred)
colnames(pred) <- gsub("\\.", " ", colnames(pred))

colnames(major_test_obs) <- c("X", "sid", celltypes.obs)

pred_true <- merge(pred, major_test_obs, by= "X")
pred_true <- pred_true %>%
  arrange(`invasive cancer.obs`)

pred_true_long <- reshape(pred_true, direction =  "long",
                          varying = list(celltypes, celltypes.obs),
                          timevar = "celltype", times = celltypes,
                          v.names = c("proportion.pred", "proportion.obs"),
                          idvar = colnames(pred_true)[c(1:8,13)])
rownames(pred_true_long) <- NULL

X_order = unlist(unique(pred_true_long[pred_true_long$celltype == "Cancer Epithelial","X"]))
pred_true_long$X <- factor(pred_true_long$X, levels=X_order)

slope <- pred_true_long %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(slope = map_dbl(model, ~summary(.)$coefficients[2,1]),
         model = NULL) %>%
  as.data.frame(.)

ggplot(data = slope,
       aes(x = batch_size, y = slope,
           color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3, alpha = 0.9, size = 1) +
  labs(x = "Batch size of training",
       y = "Slope of OLS",
       shape="Dense layer units", colour="Dropout proportion") +
  facet_grid(celltype~optimizer+learning_rate,
             labeller = labeller(learning_rate = c("0.0001" = "Learning rate = 0.0001",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.1" = "Learning rate = 0.01"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9))
ggsave("NN_regression_Wu_CV_1160920F_slope1.png", height = 12, width = 10)

slope_reduced <- slope[!((slope$optimizer == "Adam" & slope$learning_rate == "0.01")),]

ggplot(data = slope_reduced,
       aes(x = batch_size, y = slope,
           color = celltype)) +
  geom_point(size = 1.2) +
  labs(x = "Batch size",
       y = "Slope of OLS",
       shape="Learning rate", colour="Cell type") +
  facet_grid(optimizer+learning_rate~dropout_rate+dense_layer_size,
             labeller = labeller(learning_rate = c("0.01" = "Learning rate = 0.01",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.0001" = "Learning rate = 0.0001"),
                                 dropout_rate = c("0.0" = "No dropout layer",
                                                  "0.2" = "Dropout proportion = 0.2",
                                                  "0.5" = "Dropout proportion = 0.5"),
                                 dense_layer_size = c("0" = "No dense layer",
                                                      "256" = "Dense layer units = 256",
                                                      "512" = "Dense layer units = 512"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9)) +
  scale_color_manual("Cell type", values = colors) +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_slope2.png", height = 5, width = 11)

max_slope <- slope %>%
  group_by(celltype) %>%
  dplyr::filter(slope == max(slope, na.rm=TRUE)) %>%
  arrange(slope)

kbl(max_slope[,c(6:7, 1:5)], format = "latex",
    caption = "Parameters of the neutral network that gave the highest slope of best-fit line for each cell type in the training dataset.",
    col.names = c("Cell type", "Slope", "Optimizer", "Batch size",
                  "Learning rate", "Dropout Proportion", "Dense layer Units"),
    escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

tmp <- subset(pred_true_long,
              optimizer == "Adam" &
                batch_size == "16" & learning_rate == "0.0001" &
                dropout_rate == "0.0" & dense_layer_size == "256")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)

r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(celltype) %>%
  summarise(max_x = max(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
max$max_x <- 0.65*max$max_x
max$delta_y <- 0.1*max$range_y
r2_tmp <- merge(r2_tmp, max, by = "celltype")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  geom_point(shape = 21, aes(fill = sid), alpha = 0.5, size = 1.5, stroke=NA) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = max_x, y = (min_y+delta_y)),
            size = 5,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 256",
       x="Observed Proportion",
       y="Predicted Proportion",) +
  theme(plot.title = element_text(size=23, face = "bold"),
        strip.text.x = element_text(size=17),
        legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(1, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=16), axis.title.x = element_text(size=17),
        axis.title.y= element_text(size=19),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free") +
  ylim(c(0,1)) +
  xlim(c(0,1))
ggsave("NN_regression_Wu_CV_1160920F_scatterplot.png", height = 12, width = 12)


tmp_2 <- reshape(
  tmp,
  direction =  "long",
  varying = c("proportion.pred", "proportion.obs"),
  timevar = "type",
  times = c("pred", "obs"),
  v.names = c("proportion"),
  idvar = colnames(tmp)[1:9]
)

tmp_2 <- tmp_2 %>%
  group_by(celltype) %>%
  arrange(proportion, .by_group = TRUE)

X_order = unlist(unique(tmp_2[tmp_2$type == "obs" & tmp_2$celltype == "Cancer Epithelial","X"]))
tmp_2$X <- factor(tmp_2$X, levels=X_order)
tmp_2$celltype <- factor(tmp_2$celltype, levels = sort(celltypes))
tmp_2$celltype <- relevel(tmp_2$celltype, "Cancer Epithelial")

ggplot(tmp_2, aes(fill=factor(celltype, levels=rev(levels(celltype))),
                  color=factor(celltype, levels=rev(levels(celltype))),
                  y=proportion, x=X)) +
  geom_bar(position="stack", stat="identity", width = 1) +
  facet_wrap(~type,
             labeller = labeller(type = c("obs" = "Observed",
                                          "pred" = "Predicted"))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm'),
        strip.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  scale_fill_manual("Cell type", values = colors) +
  scale_color_manual("Cell type", values = colors) +
  labs(x = "Patches in Testing Data", y = "Proportion")
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_composition.png", height = 3.5, width = 10)

p_list = list()
for (i in celltypes) {
  tmp <- tmp_2[tmp_2$celltype ==i, ]

  ordered <- tmp %>%
    group_by(type) %>%
    arrange(proportion, .by_group = TRUE)

  X_order = unique(unlist(ordered$X))
  tmp$X <- factor(tmp$X, levels=X_order)

  p <- ggplot(tmp, aes(fill=factor(celltype), color=factor(celltype), y=proportion, x=X)) +
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(strip.text.x = element_text(size = 16),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18)) +
    scale_fill_manual("Cell type",
                      values=colors) +
    scale_color_manual("Cell type",
                       values=colors) +
    labs(title = i,
         x = "Patches in Testing Data",
         y = "Proportion") +
    facet_wrap(~type,
               labeller = labeller(type = c("obs" = "Observed",
                                            "pred" = "Predicted")))
  p_list[[i]] <- p
}

p <- grid.arrange(p_list[[1]],p_list[[2]],p_list[[3]],
                  p_list[[4]],p_list[[5]],p_list[[6]],
                  p_list[[7]],p_list[[8]],p_list[[9]],nrow = 3)
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_composition_per_type.png", p, height = 6, width = 15)

mu <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), mean)
mu$measure <- "mean"
med <- tmp_2 %>%
  group_by(celltype, type) %>%
  summarise_at(vars("proportion"), median)
med$measure <- "median"

measure <- rbind(mu, med)

ggplot(tmp_2, aes(x = proportion, color = type)) +
  geom_histogram(fill = "white", alpha = 0.05, position = "dodge", bins = 70) +
  geom_vline(data=measure, aes(xintercept=proportion, color=type, linetype = measure))+
  scale_colour_discrete("Proportion", breaks=c("obs", "pred"),
                        labels=c("Observed (from CARD)",
                                 "Predicted (from trained network)")) +
  scale_linetype_discrete("Measure", breaks=c("mean", "median"),
                          labels = c("Mean", "Median")) +
  labs(title="Optimizer = Adam, Batch size = 32, Learning rate = 0.0001,
       Dropout layer proportion = 0.2, Dense layer units = 512",
       x="Proportion",
       y = "Count") +
  theme(legend.position="right", legend.box = "vertical",
        legend.key.size = unit(1, 'cm'), legend.title = element_text(size=19),
        plot.title = element_text(size=22, face = "bold"),
        strip.text.x = element_text(size=20),
        legend.text = element_text(size=17),
        axis.title.x = element_text(size=21),
        axis.title.y= element_text(size=21),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free") +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1160920F_obs_vs_pred_histogram.png", height = 11, width = 15)

major_patho <- data.frame(X = major_test_obs$X,
                          `Invasive.Cancer.obs` = major_test_obs$`Cancer Epithelial.obs`,
                          Stroma.obs = major_test_obs$Endothelial.obs + major_test_obs$PVL.obs + major_test_obs$CAFs.obs,
                          Lymphocyte.obs = major_test_obs$`B cells.obs` + major_test_obs$`T cells.obs` + major_test_obs$Plasmablasts.obs,
                          Others.obs = major_test_obs$`Normal Epithelial.obs` + major_test_obs$Myeloid.obs)

pred_patho <- data.frame(`Invasive.Cancer.pred` = pred$`Cancer Epithelial`,
                         Stroma.pred = pred$Endothelial + pred$PVL + pred$CAFs,
                         Lymphocyte.pred =  pred$`T cells` + pred$`B cells` + pred$Plasmablasts,
                         Others.pred = pred$`Normal Epithelial` + pred$Myeloid)

pred_patho <- cbind(pred[,c(1:8)], pred_patho)
pred_true_patho <- merge(pred_patho, major_patho, by= "X")

pred_true_patho_long <- reshape(pred_true_patho, direction =  "long",
                                varying = list(c("Invasive.Cancer.obs",
                                                 "Stroma.obs",
                                                 "Lymphocyte.obs",
                                                 "Others.obs"),
                                               c("Invasive.Cancer.pred",
                                                 "Stroma.pred",
                                                 "Lymphocyte.pred",
                                                 "Others.pred")),
                                timevar = "pathology",
                                times = c("Invasive Cancer", "Stroma", "Lymphocyte", "Others"),
                                v.names = c("proportion.obs", "proportion.pred"),
                                idvar = colnames(pred_true_patho)[1:8])
rownames(pred_true_patho_long) <- NULL

tmp <-  subset(pred_true_patho_long,
               optimizer == "Adam" &
                 batch_size == "16" & learning_rate == "0.0001" &
                 dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, pathology) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(pathology) %>%
  summarise(range_x = max(proportion.obs)-min(proportion.obs),
            min_x = min(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
r2_tmp <- merge(r2_tmp, max, by = "pathology")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 512",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~pathology, nrow = 1, scales = "free") +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_scatterplot_patho.png", height = 3, width = 10)

tmp_nine <- tmp


###### TASK: evaluate_NN_regression_Wu_CV_11 ----
create_levels <- function(df){
  colnames(df)[1:7] <- c("base", "image", "optimizer", "batch_size", "learning_rate",
                         "dropout_rate", "dense_layer_size")
  rownames(df) <- NULL
  df$image <- factor(df$image, levels = c("o"))
  df$optimizer <- factor(df$optimizer, levels = c("Adam"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001", "0.001", "0.1"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))
  
  return(df)
}


# path for accuracy and loss plots
metrics_dir = "/Users/zhiningsui/GitHub/STpath/output/Wu_2021/regression/eval_Wu_CV_11_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial/"
# path for cell type proportions obtained by CARD
dvn_dir = "/Users/zhiningsui/GitHub/STpath/output/Wu_2021/regression/"

# load the Observed proportions
major_train <- read.csv(file.path(dvn_dir, "training_Wu_CV_11_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv"))[,-1]
major_val <- read.csv(file.path(dvn_dir, "validation_Wu_CV_11_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv"))[,-1]
major_test <- read.csv(file.path(dvn_dir, "testing_Wu_CV_11_Endothelial_CAFs_PVL_B.cells_T.cells_Myeloid_Normal.Epithelial_Plasmablasts_Cancer.Epithelial.csv"))[,-1]
colnames(major_train) <- gsub("\\.", " ", colnames(major_train))
colnames(major_val) <- gsub("\\.", " ", colnames(major_val))
colnames(major_test) <- gsub("\\.", " ", colnames(major_test))

celltypes <- colnames(major_train)[2:10]
celltypes.obs <- paste(celltypes, "obs", sep = ".")
celltypes.pred <- paste(celltypes, "pred", sep = ".")

major_train_obs <- major_train[, c("X", celltypes)]
major_train_long <- major_train_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Training data")

major_val_obs <- major_val[, c("X", celltypes)]
major_val_long <- major_val_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Validation data")

major_test_obs <- major_test[, c("X", "sid", celltypes)]
major_test_long <- major_test_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Testing data")

major_long <- bind_rows(major_train_long, major_val_long, major_test_long)

p <- ggplot(major_long,
            aes(x=factor(name), y=value, color=factor(name))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~factor(dataset,
                     levels = c("Training data", "Validation data", "Testing data")),
             scales = "free_x")

ggsave("NN_regression_CV_Wu_10x_composition.png", p, width = 6.5, height = 2.5)

colMeans(major_train[,2:10]) 
colMeans(major_test[,2:10]) 

val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Loss",
       x = "Batch size",
       y = "MSE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_mae,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Metrics",
       x = "Batch size",
       y = "MAE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")


mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("NN_regression_Wu_CV_1160920F_val_eval.png", p, width = 13, height = 5)

pred_files = list.files(metrics_dir, pattern="test_pred_")
pred = NULL
for(f1 in pred_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][3:9]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
}
pred <- create_levels(pred)
colnames(pred) <- gsub("\\.", " ", colnames(pred))

colnames(major_test_obs) <- c("X", "sid", celltypes.obs)

pred_true <- merge(pred, major_test_obs, by= "X")
pred_true <- pred_true %>%
  arrange(`invasive cancer.obs`)

pred_true_long <- reshape(pred_true, direction =  "long",
                          varying = list(celltypes, celltypes.obs),
                          timevar = "celltype", times = celltypes,
                          v.names = c("proportion.pred", "proportion.obs"),
                          idvar = colnames(pred_true)[c(1:8,18)])
rownames(pred_true_long) <- NULL

X_order = unlist(unique(pred_true_long[pred_true_long$celltype == "Cancer Epithelial","X"]))
pred_true_long$X <- factor(pred_true_long$X, levels=X_order)

slope <- pred_true_long %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(slope = map_dbl(model, ~summary(.)$coefficients[2,1]),
         model = NULL) %>%
  as.data.frame(.)

ggplot(data = slope,
       aes(x = batch_size, y = slope,
           color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3, alpha = 0.9, size = 1) +
  labs(x = "Batch size of training",
       y = "Slope of OLS",
       shape="Dense layer units", colour="Dropout proportion") +
  facet_grid(celltype~optimizer+learning_rate,
             labeller = labeller(learning_rate = c("0.0001" = "Learning rate = 0.0001",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.1" = "Learning rate = 0.01"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9))
ggsave("NN_regression_Wu_CV_1160920F_slope1.png", height = 12, width = 10)

slope_reduced <- slope[!((slope$optimizer == "Adam" & slope$learning_rate == "0.01")),]

ggplot(data = slope_reduced,
       aes(x = batch_size, y = slope,
           color = celltype)) +
  geom_point(size = 1.2) +
  labs(x = "Batch size",
       y = "Slope of OLS",
       shape="Learning rate", colour="Cell type") +
  facet_grid(optimizer+learning_rate~dropout_rate+dense_layer_size,
             labeller = labeller(learning_rate = c("0.01" = "Learning rate = 0.01",
                                                   "0.001" = "Learning rate = 0.001",
                                                   "0.0001" = "Learning rate = 0.0001"),
                                 dropout_rate = c("0.0" = "No dropout layer",
                                                  "0.2" = "Dropout proportion = 0.2",
                                                  "0.5" = "Dropout proportion = 0.5"),
                                 dense_layer_size = c("0" = "No dense layer",
                                                      "256" = "Dense layer units = 256",
                                                      "512" = "Dense layer units = 512"))) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(size=9)) +
  scale_color_manual("Cell type", values = colors) +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_slope2.png", height = 5, width = 11)

max_slope <- slope %>%
  group_by(celltype) %>%
  dplyr::filter(slope == max(slope, na.rm=TRUE)) %>%
  arrange(slope)

kbl(max_slope[,c(6:7, 1:5)], format = "latex",
    caption = "Parameters of the neutral network that gave the highest slope of best-fit line for each cell type in the training dataset.",
    col.names = c("Cell type", "Slope", "Optimizer", "Batch size",
                  "Learning rate", "Dropout Proportion", "Dense layer Units"),
    escape = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")

tmp <- subset(pred_true_long,
              optimizer == "Adam" &
                batch_size == "16" & learning_rate == "0.0001" &
                dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)

r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(celltype) %>%
  summarise(max_x = max(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
max$max_x <- 0.65*max$max_x
max$delta_y <- 0.1*max$range_y
r2_tmp <- merge(r2_tmp, max, by = "celltype")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  geom_point(shape = 21, aes(fill = sid), alpha = 0.5, size = 1.5, stroke=NA) +
  geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = max_x, y = (min_y+delta_y)),
            size = 5,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 256",
       x="Observed Proportion",
       y="Predicted Proportion",) +
  theme(plot.title = element_text(size=23, face = "bold"),
        strip.text.x = element_text(size=17),
        legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(1, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=16), axis.title.x = element_text(size=17),
        axis.title.y= element_text(size=19),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 3, scales = "free")
ggsave("NN_regression_Wu_CV_1160920F_scatterplot.png", height = 12, width = 12)


major_patho <- data.frame(X = major_test_obs$X,
                          `Invasive.Cancer.obs` = major_test_obs$`Cancer Epithelial.obs`,
                          Stroma.obs = major_test_obs$Endothelial.obs + major_test_obs$PVL.obs + major_test_obs$CAFs.obs,
                          Lymphocyte.obs = major_test_obs$`B cells.obs` + major_test_obs$`T cells.obs` + major_test_obs$Plasmablasts.obs,
                          Others.obs = major_test_obs$`Normal Epithelial.obs` + major_test_obs$Myeloid.obs)

pred_patho <- data.frame(`Invasive.Cancer.pred` = pred$`Cancer Epithelial`,
                         Stroma.pred = pred$Endothelial + pred$PVL + pred$CAFs,
                         Lymphocyte.pred =  pred$`T cells` + pred$`B cells` + pred$Plasmablasts,
                         Others.pred = pred$`Normal Epithelial` + pred$Myeloid)

pred_patho <- cbind(pred[,c(1:8)], pred_patho)
pred_true_patho <- merge(pred_patho, major_patho, by= "X")

pred_true_patho_long <- reshape(pred_true_patho, direction =  "long",
                                varying = list(c("Invasive.Cancer.obs",
                                                 "Stroma.obs",
                                                 "Lymphocyte.obs",
                                                 "Others.obs"),
                                               c("Invasive.Cancer.pred",
                                                 "Stroma.pred",
                                                 "Lymphocyte.pred",
                                                 "Others.pred")),
                                timevar = "pathology",
                                times = c("Invasive Cancer", "Stroma", "Lymphocyte", "Others"),
                                v.names = c("proportion.obs", "proportion.pred"),
                                idvar = colnames(pred_true_patho)[1:8])
rownames(pred_true_patho_long) <- NULL

tmp <-  subset(pred_true_patho_long,
               optimizer == "Adam" &
                 batch_size == "16" & learning_rate == "0.0001" &
                 dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, pathology) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(pathology) %>%
  summarise(range_x = max(proportion.obs)-min(proportion.obs),
            min_x = min(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
r2_tmp <- merge(r2_tmp, max, by = "pathology")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 512",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~pathology, nrow = 1, scales = "free") +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_scatterplot_patho.png", height = 3, width = 10)

tmp_nine <- tmp



###### TASK: evaluate_NN_regression_Janesick ----
create_levels <- function(df){
  colnames(df)[1:7] <- c("base", "image", "optimizer", "batch_size", "learning_rate",
                         "dropout_rate", "dense_layer_size")
  rownames(df) <- NULL
  df$image <- factor(df$image, levels = c("o"))
  df$optimizer <- factor(df$optimizer, levels = c("Adam"))
  df$batch_size <- factor(df$batch_size, levels = c("16", "32", "64", "128"))
  df$learning_rate <- factor(df$learning_rate, levels = c("0.0001", "0.001", "0.1"))
  df$dropout_rate <- factor(df$dropout_rate, levels = c("0.0", "0.2", "0.5"))
  df$dense_layer_size <- factor(df$dense_layer_size, levels = c("0", "256", "512"))
  
  return(df)
}


# path for accuracy and loss plots
metrics_dir = "/Users/zhiningsui/GitHub/STpath/output/Janesick_2023/regression/eval_Janesick_VisiumS1_Lymphocyte.xenium_Myeloid.xenium_Tumor.xenium_Stroma.xenium_Others.xenium/"
# path for cell type proportions obtained by CARD
dvn_dir = "/Users/zhiningsui/GitHub/STpath/output/Janesick_2023/regression/"

# load the Observed proportions
major_train <- read.csv(file.path(dvn_dir, "training_Janesick_VisiumS1_Lymphocyte.xenium_Myeloid.xenium_Tumor.xenium_Stroma.xenium_Others.xenium.csv"), check.names = F)[,-1]
major_val <- read.csv(file.path(dvn_dir, "validation_Janesick_VisiumS1_Lymphocyte.xenium_Myeloid.xenium_Tumor.xenium_Stroma.xenium_Others.xenium.csv"), check.names = F)[,-1]
major_test <- read.csv(file.path(dvn_dir, "testing_Janesick_VisiumS1_Lymphocyte.xenium_Myeloid.xenium_Tumor.xenium_Stroma.xenium_Others.xenium.csv"), check.names = F)[,-1]

celltypes <- colnames(major_train)[2:6]
celltypes.obs <- paste(celltypes, "obs", sep = ".")
celltypes.pred <- paste(celltypes, "pred", sep = ".")

major_train_obs <- major_train[, c("X", celltypes)]
major_train_long <- major_train_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Training data")

major_val_obs <- major_val[, c("X", celltypes)]
major_val_long <- major_val_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Validation data")

major_test_obs <- major_test[, c("X", celltypes)]
major_test_long <- major_test_obs %>%
  pivot_longer(cols = all_of(celltypes)) %>%
  mutate(dataset = "Testing data")

major_long <- bind_rows(major_train_long, major_val_long, major_test_long)

p <- ggplot(major_long,
            aes(x=factor(name), y=value, color=factor(name))) +
  geom_boxplot(outlier.alpha = 0.3, outlier.size = 1) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 10.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.4, "cm")) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~factor(dataset,
                     levels = c("Training data", "Validation data", "Testing data")),
             scales = "free_x")

ggsave("NN_regression_CV_Wu_10x_composition.png", p, width = 6.5, height = 2.5)


val_files = list.files(metrics_dir, pattern="metrics_history_")
val_eval = NULL
for(f1 in val_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][4:10]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  val_eval <- rbind(val_eval, cbind(params, p1[nrow(p1),]))
}

val_eval <- create_levels(val_eval)

p1 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_loss,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Loss",
       x = "Batch size",
       y = "MSE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")

p2 <- ggplot(data = val_eval,
             aes(x = batch_size, y = val_mae,
                 color = dropout_rate, shape = dense_layer_size)) +
  geom_beeswarm(cex = 3.5) +
  labs(subtitle = "Validation Metrics",
       x = "Batch size",
       y = "MAE") +
  ggh4x::facet_grid2(optimizer~learning_rate,
                     scales = "free_y", independent = "y",
                     labeller = labeller(dropout_rate = c("0.0" = "No Dropout Layer",
                                                          "0.2" = "Dropout Proportion = 0.2",
                                                          "0.5" = "Dropout Proportion = 0.5"),
                                         learning_rate = c("0.1" = "Learning Rate = 0.01",
                                                           "0.001" = "Learning Rate = 0.001",
                                                           "0.0001" = "Learning Rate = 0.0001"))) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  scale_colour_discrete("Dropout proportion") +
  scale_shape_discrete("Dense layer units")


mylegend <- g_legend(p1)
p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                              p2 + theme(legend.position="none"), nrow = 1), mylegend, heights=c(13, 1))
ggsave("NN_regression_Wu_CV_1160920F_val_eval.png", p, width = 13, height = 5)

pred_files = list.files(metrics_dir, pattern="test_pred_")
pred = NULL
for(f1 in pred_files){
  s1 <- gsub(".csv", "",f1)
  params <- as.data.frame(t(str_split(s1, pattern = "_")[[1]][4:10]))
  p1 <- read.csv(file.path(metrics_dir, f1))
  params <- params[rep(1, each = nrow(p1)),]
  pred <- rbind(pred, cbind(params, p1))
}
pred <- create_levels(pred)

colnames(major_test_obs) <- c("X", celltypes.obs)

pred_true <- merge(pred, major_test_obs, by= "X")

pred_true_long <- reshape(pred_true, direction =  "long",
                          varying = list(celltypes, celltypes.obs),
                          timevar = "celltype", times = celltypes,
                          v.names = c("proportion.pred", "proportion.obs"),
                          idvar = colnames(pred_true)[c(1:8)])
rownames(pred_true_long) <- NULL

tmp <- subset(pred_true_long,
              optimizer == "Adam" &
                batch_size == "32" & learning_rate == "0.0001" &
                dropout_rate == "0.0" & dense_layer_size == "256")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, celltype) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)

r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(celltype) %>%
  summarise(max_x = max(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
max$max_x <- 0.65*max$max_x
max$delta_y <- 0.1*max$range_y
r2_tmp <- merge(r2_tmp, max, by = "celltype")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1.5, stroke=NA) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = max_x, y = (min_y+delta_y)),
            size = 5,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 32, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 256",
       x="Observed Proportion",
       y="Predicted Proportion",) +
  theme(plot.title = element_text(size=23, face = "bold"),
        strip.text.x = element_text(size=17),
        legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(1, 'cm'), legend.title = element_blank(),
        legend.text = element_text(size=16), axis.title.x = element_text(size=17),
        axis.title.y= element_text(size=19),
        axis.text = element_text(size=15, color = "black")) +
  facet_wrap(~celltype, nrow = 1, scales = "free")

ggsave("NN_regression_Wu_CV_1160920F_scatterplot.png", height = 12, width = 12)


major_patho <- data.frame(X = major_test_obs$X,
                          `Invasive.Cancer.obs` = major_test_obs$`Cancer Epithelial.obs`,
                          Stroma.obs = major_test_obs$Endothelial.obs + major_test_obs$PVL.obs + major_test_obs$CAFs.obs,
                          Lymphocyte.obs = major_test_obs$`B cells.obs` + major_test_obs$`T cells.obs` + major_test_obs$Plasmablasts.obs,
                          Others.obs = major_test_obs$`Normal Epithelial.obs` + major_test_obs$Myeloid.obs)

pred_patho <- data.frame(`Invasive.Cancer.pred` = pred$`Cancer Epithelial`,
                         Stroma.pred = pred$Endothelial + pred$PVL + pred$CAFs,
                         Lymphocyte.pred =  pred$`T cells` + pred$`B cells` + pred$Plasmablasts,
                         Others.pred = pred$`Normal Epithelial` + pred$Myeloid)

pred_patho <- cbind(pred[,c(1:8)], pred_patho)
pred_true_patho <- merge(pred_patho, major_patho, by= "X")

pred_true_patho_long <- reshape(pred_true_patho, direction =  "long",
                                varying = list(c("Invasive.Cancer.obs",
                                                 "Stroma.obs",
                                                 "Lymphocyte.obs",
                                                 "Others.obs"),
                                               c("Invasive.Cancer.pred",
                                                 "Stroma.pred",
                                                 "Lymphocyte.pred",
                                                 "Others.pred")),
                                timevar = "pathology",
                                times = c("Invasive Cancer", "Stroma", "Lymphocyte", "Others"),
                                v.names = c("proportion.obs", "proportion.pred"),
                                idvar = colnames(pred_true_patho)[1:8])
rownames(pred_true_patho_long) <- NULL

tmp <-  subset(pred_true_patho_long,
               optimizer == "Adam" &
                 batch_size == "16" & learning_rate == "0.0001" &
                 dropout_rate == "0.0" & dense_layer_size == "512")

r2_tmp <- tmp %>%
  group_by(optimizer, batch_size, learning_rate, dropout_rate, dense_layer_size, pathology) %>%
  summarise(model = list(lm(proportion.pred ~ proportion.obs))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2_tmp$r2.adj <- sprintf("italic(R^2) == %.3f", r2_tmp$r2.adj)
max <- tmp %>%
  group_by(pathology) %>%
  summarise(range_x = max(proportion.obs)-min(proportion.obs),
            min_x = min(proportion.obs),
            range_y = max(proportion.pred)-min(proportion.pred),
            min_y = min(proportion.pred))
r2_tmp <- merge(r2_tmp, max, by = "pathology")

ggplot(tmp, aes(x = proportion.obs, y = proportion.pred)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2_tmp,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(title="Optimizer = Adam, Batch size = 16, Learning rate = 0.0001,
       Dropout layer proportion = 0.0, Dense layer units = 512",
       x="Observed Proportion",
       y="Predicted Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~pathology, nrow = 1, scales = "free") +
  guides(colour = guide_legend(nrow = 1))
ggsave("NN_regression_Wu_CV_1160920F_scatterplot_patho.png", height = 3, width = 10)

tmp_nine <- tmp





files = list.files()
times_df <- data.frame()
time_strings <- NULL
i=1
for (f1 in files) {
  print(i)
  p1 <- readLines(f1)
  time <- p1[str_detect(p1, "/step")]
  time <- time[str_detect(time, "val_loss")][-1]
  times.epoch <- NULL
  times.step <- NULL
  
  for (j in 1:length(time)) {
    substrings <- unlist(strsplit(time[j], " \\[|\\] | - "))
    steps <- unlist(strsplit(substrings[1], "/"))[1]
    t1 <- gsub("- ", "", substrings[3])
    
    time.per.epoch <- strsplit(t1, " ")[[1]][1]
    time.per.epoch <- as.numeric(sub("s.*", "", time.per.epoch))
    
    time.per.step <- strsplit(t1, " ")[[1]][2]
    # Check if milliseconds exist in the string
    
    total_seconds <- seconds + milliseconds / 1000
    
    if (!grepl("ms", time.per.step)) {
      time.per.step <- as.numeric(sub("s.*", "", time.per.step)) * 1000
    } else {
      time.per.step <- as.numeric(sub("([0-9]+)ms.*", "\\1", time.per.step))
    }
    times.epoch <- c(times.epoch, time.per.epoch)
    times.step <- c(times.step, time.per.step)
  }

  times_df <- rbind(times_df, c(strsplit(f1, "_")[[1]][c(2, 7:10)], steps, 
                                mean(times.epoch),
                                sd(times.epoch),
                                mean(times.step),
                                sd(times.step)))
  i=i+1
}

names(times_df) <- c("sample", "batch", "lr", "drop", "dense", "step", "time.per.epoch", "sd.time.per.epoch",  "time.per.step", "sd.time.per.step")
times_df$batch <- as.integer(times_df$batch)
times_df$step <- as.integer(times_df$step)
times_df$time.per.epoch = as.numeric(times_df$time.per.epoch)
times_df$sd.time.per.epoch = as.numeric(times_df$sd.time.per.epoch)
times_df$time.per.step = as.numeric(times_df$time.per.step)
times_df$size <- times_df$batch * times_df$step

times_df$batch <- factor(times_df$batch, levels = sort(as.numeric(unique(times_df$batch))))

times_df <- na.omit(times_df)
times_df <- times_df[times_df$time.per.epoch < 200,]

mean_df <- times_df %>%
  group_by(batch) %>%
  summarise(mean = mean(time.per.epoch))
mean_df$mean <- round(mean_df$mean, 1)
mean_df$label <- paste0(mean_df$mean, "s")

ggplot(times_df, aes(x = batch, y = time.per.epoch)) + 
  geom_boxplot() +
  labs(x = "Batch Size", y = "Time per Epoch in Each Task", 
       title = "Training Data Size: 8512 to 10544")+
  geom_text(data = mean_df, aes(y = mean + 20,
                                label = label),
            color = "red") +
  theme_bw()

ggplot(times_df, aes(x = batch, y = sd.time.per.epoch)) + 
  geom_boxplot() +
  labs(x = "Batch Size", y = "Standard Deviation of Time per Epoch in Each Task", 
       title = "Training Data Size: 8512 to 10544")


