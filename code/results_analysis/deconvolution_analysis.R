rm(list = ls())
library(knitr)
library(tidyverse)
library(stringr)
library(car)
library(yarrr)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(grid)
library(gridExtra)
library(reshape2)
library(kableExtra)
library(styler)
library(scales)
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


########### READ CARD RESULTS FOR ALL DATASETS ############

## Janesick et al
meta_dir = "../../data/Janesick_2023/"
dvn_dir = "../../output/Janesick_2023/deconvolution/"

## 10x Genomics
meta_dir = "/Users/zhiningsui/GitHub/st2image_data/10X/data/Visium_Human_Breast/analysis/clustering/graphclust/"
dvn_dir = "/Users/zhiningsui/GitHub/st2image_data/10X/data/deconvolution"
major_10X_FF = read.csv(file.path(dvn_dir, "Proportion_Visium_Human_Breast.csv"))
meta_dir = "/Users/zhiningsui/GitHub/st2image_data/10X/data/Visium_FFPE_Human_Breast_Cancer/analysis/clustering/graphclust/"
dvn_dir = "/Users/zhiningsui/GitHub/st2image_data/10X/data/deconvolution"
major_10X_FFPE = read.csv(file.path(dvn_dir, "Proportion_Visium_FFPE_Human_Breast_Cancer.csv"))

## We et al
meta_dir = "../../st2image_data/Wu_2021/data/metadata/"
dvn_dir = "../../output/Wu_2021/deconvolution/"


# Janesick et al ----------------------------------------------------------

dvn_prop_files = list.files(dvn_dir, pattern="_Annotation.csv")
Janesick_prop = NULL
for(f1 in dvn_prop_files){
  p1 <- read.csv(file.path(dvn_dir, f1), check.names = F)
  colnames(p1)[1] <- "X" 
  s1 <- gsub("_Annotation.csv", "",f1)
  s2 <- gsub("Proportion_","",s1)
  p1$X = paste(s2, p1$X, sep="_")
  Janesick_prop <- rbind(Janesick_prop, p1)
}
celltypes <- sort(colnames(Janesick_prop)[2:(ncol(Janesick_prop)-1)])
Janesick_prop <- Janesick_prop[, c("X", "sid", celltypes)]
Janesick_prop$max <- apply(Janesick_prop[,celltypes], 1, max)
Janesick_prop$wmax <- apply(Janesick_prop[,celltypes], 1, which.max)
Janesick_prop$max.celltype <- celltypes[Janesick_prop$wmax]

number_max_prop <- data.frame(celltype = celltypes)
t1 <- as.data.frame(table(celltypes[Janesick_prop$wmax]))
number_max_prop <- merge(number_max_prop, t1, by.x = "celltype", by.y = "Var1", all = T)
write.csv(number_max_prop, "number_of_max_celltype.csv")

Janesick_long <- reshape(Janesick_prop, direction =  "long", 
                         varying = celltypes,
                         timevar = "celltype", 
                         times = celltypes,
                         v.names = c("proportion"), 
                         idvar = c("X", "sid"))
rownames(Janesick_long) <- NULL
Janesick_long$celltype <- factor(Janesick_long$celltype, levels = celltypes)
samples <- sort(unique(Janesick_plot$sid))

ggplot(Janesick_long, 
       aes(x=celltype, y=proportion, color=celltype)) + 
  geom_boxplot(width = 0.5, outlier.alpha = 0.3, outlier.size = 0.8) + 
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 11)) +
  # scale_fill_manual("Cell type", values = colors) +
  # scale_color_manual("Cell type", values = colors) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~sid, nrow = 2) +
  guides(color=guide_legend(nrow=1)) +
  coord_flip()
ggsave("boxplot_proportions.png", width = 10, height = 5)

########## Consistency with True Proportion #########

Janesick_prop_xenium <- read.csv(file.path(meta_dir, "xenium_spotbinned_celltypes.tsv"),
                               sep = "\t", check.names = F)
colnames(Janesick_prop_xenium)[1] <- "X"
Janesick_prop_xenium$X <- paste0("VisiumS1_", Janesick_prop_xenium$X)
colnames(Janesick_prop_xenium) <- gsub("\\_", " ", colnames(Janesick_prop_xenium))
celltypes.xenium <- sort(colnames(Janesick_prop_xenium)[-1])
Janesick_prop_xenium <- Janesick_prop_xenium[, c("X", celltypes.xenium)]

sum(Janesick_prop_xenium$Unlabeled > 0)

colnames(Janesick_prop_xenium) <- c("X", paste0(celltypes.xenium, ".xenium"))
Janesick_prop_card <- Janesick_prop[c("X", celltypes)]
Janesick_prop_card$Unlabeled <- 0
colnames(Janesick_prop_card) <- c("X", paste0(celltypes.xenium, ".card"))

Janesick_prop_xenium_card <- merge(Janesick_prop_xenium, Janesick_prop_card, by = "X")
Janesick_prop_xenium_card_long <- reshape(Janesick_prop_xenium_card, direction =  "long",
                                       varying = list(paste0(celltypes.xenium, ".xenium"),
                                                      paste0(celltypes.xenium, ".card")),
                                timevar = "celltype",
                                times = celltypes.xenium,
                                v.names = c("proportion.xenium", "proportion.card"),
                                idvar = "X")
rownames(Janesick_prop_xenium_card_long) <- NULL

r2 <- Janesick_prop_xenium_card_long %>%
  group_by(celltype) %>%
  summarise(model = list(lm(proportion.card ~ proportion.xenium))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2$r2.adj <- sprintf("italic(R^2) == %.3f", r2$r2.adj)
max <- Janesick_prop_xenium_card_long %>%
  group_by(celltype) %>%
  summarise(range_x = max(proportion.xenium)-min(proportion.xenium),
            min_x = min(proportion.xenium),
            range_y = max(proportion.card)-min(proportion.card),
            min_y = min(proportion.card))
r2 <- merge(r2, max, by = "celltype")

ggplot(Janesick_prop_xenium_card_long, aes(x = proportion.xenium, y = proportion.card)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = xenium,
            hjust = 0,
            color = "red") +
  labs(title = "Assign 0 to Unlabeled in CARD proportion",
       x="xenium Proportion",
       y="CARD Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~celltype, nrow = 4, scales = "free") +
  guides(colour = guide_legend(nrow = 1))

ggsave("xenium_vs_card_prop_1.png", width = 12, height = 8)

Janesick_prop_xenium_card_nounlabeled <- Janesick_prop_xenium_card[Janesick_prop_xenium_card$Unlabeled.xenium == 0,]

Janesick_prop_xenium_card_nounlabeled_long <- reshape(Janesick_prop_xenium_card_nounlabeled, direction =  "long",
                                       varying = list(paste0(celltypes.xenium, ".xenium"),
                                                      paste0(celltypes.xenium, ".card")),
                                       timevar = "celltype",
                                       times = celltypes.xenium,
                                       v.names = c("proportion.xenium", "proportion.card"),
                                       idvar = "X")
rownames(Janesick_prop_xenium_card_nounlabeled_long) <- NULL
Janesick_prop_xenium_card_nounlabeled_long <- Janesick_prop_xenium_card_nounlabeled_long[Janesick_prop_xenium_card_nounlabeled_long$celltype != "Unlabeled", ]
r2 <- Janesick_prop_xenium_card_nounlabeled_long %>%
  group_by(celltype) %>%
  summarise(model = list(lm(proportion.card ~ proportion.xenium))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2$r2.adj <- sprintf("italic(R^2) == %.3f", r2$r2.adj)
max <- Janesick_prop_xenium_card_nounlabeled_long %>%
  group_by(celltype) %>%
  summarise(range_x = max(proportion.xenium)-min(proportion.xenium),
            min_x = min(proportion.xenium),
            range_y = max(proportion.card)-min(proportion.card),
            min_y = min(proportion.card))
r2 <- merge(r2, max, by = "celltype")

ggplot(Janesick_prop_xenium_card_nounlabeled_long, aes(x = proportion.xenium, y = proportion.card)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = xenium,
            hjust = 0,
            color = "red") +
  labs(title = "Exclude barcodes with nonzero Unlabeled in xenium proportion",
       x="xenium Proportion",
       y="CARD Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~celltype, nrow = 4, scales = "free") +
  guides(colour = guide_legend(nrow = 1))

ggsave("xenium_vs_card_prop_2.png", width = 12, height = 8)





Janesick_prop_xenium <- read.csv(file.path(meta_dir, "xenium_spotbinned_celltypes.tsv"),
                               sep = "\t", check.names = F)
colnames(Janesick_prop_xenium)[1] <- "X"
Janesick_prop_xenium$X <- paste0("VisiumS1_", Janesick_prop_xenium$X)
colnames(Janesick_prop_xenium) <- gsub("\\_", " ", colnames(Janesick_prop_xenium))

ggplot(Janesick_prop_xenium, aes(x = "", y=Unlabeled)) +
  geom_boxplot() +
  geom_jitter()

Janesick_prop_xenium_normalized <- Janesick_prop_xenium[Janesick_prop_xenium$Unlabeled < 0.3, -which(names(Janesick_prop_xenium) %in% c("Unlabeled"))]  
Janesick_prop_xenium_normalized[,-1] <- Janesick_prop_xenium_normalized[,-1] / rowSums(Janesick_prop_xenium_normalized[,-1])

Janesick_prop_xenium_normalized$DCIS <- Janesick_prop_xenium_normalized$`DCIS 1` + Janesick_prop_xenium_normalized$`DCIS 2`
Janesick_prop_xenium_normalized$DCs <- Janesick_prop_xenium_normalized$`IRF7+ DCs` + Janesick_prop_xenium_normalized$`LAMP3+ DCs`
Janesick_prop_xenium_normalized$Macrophages <- Janesick_prop_xenium_normalized$`Macrophages 1` + Janesick_prop_xenium_normalized$`Macrophages 2`
Janesick_prop_xenium_normalized$Myoepi <- Janesick_prop_xenium_normalized$`Myoepi ACTA2+` + Janesick_prop_xenium_normalized$`Myoepi KRT15+`

Janesick_prop_xenium_normalized <- Janesick_prop_xenium_normalized[, !(names(Janesick_prop_xenium_normalized) %in% c("DCIS 1", "DCIS 2",
                                                                                                                     "IRF7+ DCs", "LAMP3+ DCs",
                                                                                                                     "Macrophages 1", "Macrophages 2",
                                                                                                                     "Myoepi ACTA2+", "Myoepi KRT15+"))]
celltypes_collapsed <- sort(colnames(Janesick_prop_xenium_normalized)[-1])
Janesick_prop_xenium_normalized <- Janesick_prop_xenium_normalized[, c("X", celltypes_collapsed)]
colnames(Janesick_prop_xenium_normalized) <- c("X", paste0(celltypes_collapsed, ".xenium"))

Janesick_prop_card <- Janesick_prop[c("X", celltypes)]
Janesick_prop_card$DCIS <- Janesick_prop_card$`DCIS 1` + Janesick_prop_card$`DCIS 2`
Janesick_prop_card$DCs <- Janesick_prop_card$`IRF7+ DCs` + Janesick_prop_card$`LAMP3+ DCs`
Janesick_prop_card$Macrophages <- Janesick_prop_card$`Macrophages 1` + Janesick_prop_card$`Macrophages 2`
Janesick_prop_card$Myoepi <- Janesick_prop_card$`Myoepi ACTA2+` + Janesick_prop_card$`Myoepi KRT15+`
Janesick_prop_card <- Janesick_prop_card[, !(names(Janesick_prop_card) %in% c("DCIS 1", "DCIS 2",
                                                                              "IRF7+ DCs", "LAMP3+ DCs",
                                                                              "Macrophages 1", "Macrophages 2",
                                                                              "Myoepi ACTA2+", "Myoepi KRT15+"))]
Janesick_prop_card <- Janesick_prop_card[, c("X", celltypes_collapsed)]
colnames(Janesick_prop_card) <- c("X", paste0(celltypes_collapsed, ".card"))

Janesick_prop_xenium_card <- merge(Janesick_prop_xenium_normalized, Janesick_prop_card, by = "X")
Janesick_prop_xenium_card_long <- reshape(Janesick_prop_xenium_card, direction =  "long",
                                       varying = list(paste0(celltypes_collapsed, ".xenium"),
                                                      paste0(celltypes_collapsed, ".card")),
                                       timevar = "celltype",
                                       times = celltypes_collapsed,
                                       v.names = c("proportion.xenium", "proportion.card"),
                                       idvar = "X")
rownames(Janesick_prop_xenium_card_long) <- NULL



# Janesick_prop_card_collapsed <- data.frame(X = Janesick_prop_card$X,
#                                            Lymphocyte.card = Janesick_prop_card$`B Cells.card` + Janesick_prop_card$`CD4+ T Cells.card` + Janesick_prop_card$`CD8+ T Cells.card`,
#                                            Myeloid.card = Janesick_prop_card$DCs.card + Janesick_prop_card$`Mast Cells.card` + Janesick_prop_card$Macrophages.card,
#                                            Tumor.card = Janesick_prop_card$`Invasive Tumor.card` + Janesick_prop_card$`Prolif Invasive Tumor.card` + Janesick_prop_card$DCIS.card,
#                                            Stroma.card = Janesick_prop_card$Stromal.card + Janesick_prop_card$Endothelial.card + Janesick_prop_card$`Perivascular-Like.card` + Janesick_prop_card$Myoepi.card,
#                                            Others.card = Janesick_prop_card$`Stromal & T Cell Hybrid.card` + Janesick_prop_card$`T Cell & Tumor Hybrid.card`)
# 
# Janesick_prop_xenium_collapsed <- data.frame(X = Janesick_prop_xenium_normalized$X,
#                                            Lymphocyte.xenium = Janesick_prop_xenium_normalized$`B Cells.xenium` + Janesick_prop_xenium_normalized$`CD4+ T Cells.xenium` + Janesick_prop_xenium_normalized$`CD8+ T Cells.xenium`,
#                                            Myeloid.xenium = Janesick_prop_xenium_normalized$DCs.xenium + Janesick_prop_xenium_normalized$`Mast Cells.xenium` + Janesick_prop_xenium_normalized$Macrophages.xenium,
#                                            Tumor.xenium = Janesick_prop_xenium_normalized$`Invasive Tumor.xenium` + Janesick_prop_xenium_normalized$`Prolif Invasive Tumor.xenium` + Janesick_prop_xenium_normalized$DCIS.xenium,
#                                            Stroma.xenium = Janesick_prop_xenium_normalized$Stromal.xenium + Janesick_prop_xenium_normalized$Endothelial.xenium + Janesick_prop_xenium_normalized$`Perivascular-Like.xenium` + Janesick_prop_xenium_normalized$Myoepi.xenium,
#                                            Others.xenium = Janesick_prop_xenium_normalized$`Stromal & T Cell Hybrid.xenium` + Janesick_prop_xenium_normalized$`T Cell & Tumor Hybrid.xenium`)
# 
Janesick_prop_card_collapsed <- data.frame(X = Janesick_prop_card$X,
                                           Lymphocyte.card = Janesick_prop_card$`B Cells.card` + Janesick_prop_card$`CD4+ T Cells.card` + Janesick_prop_card$`CD8+ T Cells.card`,
                                           Myeloid.card = Janesick_prop_card$DCs.card + Janesick_prop_card$`Mast Cells.card` + Janesick_prop_card$Macrophages.card,
                                           Tumor.card = Janesick_prop_card$`Invasive Tumor.card` + Janesick_prop_card$`Prolif Invasive Tumor.card` + Janesick_prop_card$DCIS.card,
                                           Stroma.card = Janesick_prop_card$Stromal.card + Janesick_prop_card$Endothelial.card + Janesick_prop_card$`Perivascular-Like.card`,
                                           Myoepi.card = Janesick_prop_card$Myoepi.card,
                                           Others.card = Janesick_prop_card$`Stromal & T Cell Hybrid.card` + Janesick_prop_card$`T Cell & Tumor Hybrid.card`)

Janesick_prop_xenium_collapsed <- data.frame(X = Janesick_prop_xenium_normalized$X,
                                             Lymphocyte.xenium = Janesick_prop_xenium_normalized$`B Cells.xenium` + Janesick_prop_xenium_normalized$`CD4+ T Cells.xenium` + Janesick_prop_xenium_normalized$`CD8+ T Cells.xenium`,
                                             Myeloid.xenium = Janesick_prop_xenium_normalized$DCs.xenium + Janesick_prop_xenium_normalized$`Mast Cells.xenium` + Janesick_prop_xenium_normalized$Macrophages.xenium,
                                             Tumor.xenium = Janesick_prop_xenium_normalized$`Invasive Tumor.xenium` + Janesick_prop_xenium_normalized$`Prolif Invasive Tumor.xenium` + Janesick_prop_xenium_normalized$DCIS.xenium,
                                             Stroma.xenium = Janesick_prop_xenium_normalized$Stromal.xenium + Janesick_prop_xenium_normalized$Endothelial.xenium + Janesick_prop_xenium_normalized$`Perivascular-Like.xenium`,
                                             Myoepi.xenium = Janesick_prop_xenium_normalized$Myoepi.xenium,
                                             Others.xenium = Janesick_prop_xenium_normalized$`Stromal & T Cell Hybrid.xenium` + Janesick_prop_xenium_normalized$`T Cell & Tumor Hybrid.xenium`)


Janesick_prop_xenium_card_collapsed <- merge(Janesick_prop_xenium_collapsed, Janesick_prop_card_collapsed, by = "X")
Janesick_prop_xenium_card_collapsed_long <- reshape(Janesick_prop_xenium_card_collapsed, direction =  "long",
                                          varying = list(paste0(c("Lymphocyte", "Myeloid", "Tumor", "Stroma", "Myoepi", "Others"), ".xenium"),
                                                         paste0(c("Lymphocyte", "Myeloid", "Tumor", "Stroma", "Myoepi", "Others"), ".card")),
                                          timevar = "celltype",
                                          times = c("Lymphocyte", "Myeloid", "Tumor", "Stroma", "Myoepi", "Others"),
                                          v.names = c("proportion.xenium", "proportion.card"),
                                          idvar = "X")
rownames(Janesick_prop_xenium_card_collapsed_long) <- NULL

r2 <- Janesick_prop_xenium_card_collapsed_long %>%
  group_by(celltype) %>%
  summarise(model = list(lm(proportion.card ~ proportion.xenium))) %>%
  mutate(r2.adj = map_dbl(model, ~summary(.)$adj.r.squared),
         model = NULL) %>%
  as.data.frame(.)
r2$r2.adj <- sprintf("italic(R^2) == %.3f", r2$r2.adj)
max <- Janesick_prop_xenium_card_collapsed_long %>%
  group_by(celltype) %>%
  summarise(range_x = max(proportion.xenium)-min(proportion.xenium),
            min_x = min(proportion.xenium),
            range_y = max(proportion.card)-min(proportion.card),
            min_y = min(proportion.card))
r2 <- merge(r2, max, by = "celltype")

ggplot(Janesick_prop_xenium_card_collapsed_long, aes(x = proportion.xenium, y = proportion.card)) +
  # geom_point(col = "black", alpha = 0.3, size = 0.7) +
  geom_point(shape = 21, fill = "black", alpha = 0.5, size = 1, stroke=NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 0.7) +
  # geom_smooth(method="lm", color = "blue", fill = "#6d93fc", size = 0.7) +
  geom_text(data = r2,
            aes(label = r2.adj, x = (min_x+0.7*range_x), y = (min_y+0.05*range_y)),
            size = 3,
            parse = TRUE,
            hjust = 0,
            color = "red") +
  labs(x="Xenium Proportion",
       y="CARD Proportion") +
  theme(plot.title = element_text(size=14, face = "bold"),
        strip.text.x = element_text(size = 12),
        legend.position="none", legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=10), axis.title.x = element_text(size=10),
        axis.title.y= element_text(size=10),
        axis.text = element_text(color = "black")) +
  facet_wrap(~celltype, nrow = 2, scales = "free") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_bw()

ggsave("Xenium_vs_CARD_prop_collapsed_for_visium.png", width = 8, height = 4)

Janesick_prop_xenium_card_collapsed$img_path <- paste0("/fh/scratch/delete90/sun_w/zsui/st2image_data/Janesick_2023/output/patch/GSM7782699_", 
                                                       Janesick_prop_xenium_card_collapsed$X)
write.csv(Janesick_prop_xenium_card_collapsed, "Janesick_prop_xenium_card_collapsed_6ct.csv")

########## Consistency with Annotation #########
metadata_files = list.files(meta_dir, pattern=".csv")
metadata = NULL
for(f1 in metadata_files){
  s1 <- gsub("Barcode_Type_", "",f1)
  s2 <- gsub(".csv", "", s1)
  if (s2 %in% samples){
    p1 <- read.csv(file.path(meta_dir, f1))
    colnames(p1)[1] <- "X"
    p1$X = paste(s2, p1$X, sep="_")
    metadata <- rbind(metadata, p1)
  }
}

table(metadata$Cluster, metadata$Annotation)
Janesick_plot <- merge(Janesick_prop, metadata, by = "X")

Janesick_plot_long <- melt(Janesick_plot, id.vars=c('Annotation', 'Cluster'), 
                           measure.vars=celltypes)

ggplot(data = Janesick_plot_long, 
       aes(x = Annotation, y = value, color = Annotation)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.5) + 
  theme(legend.position="none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        strip.text = element_text(size=8)) + 
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_wrap(~variable, scales = "free")
ggsave("consistency_boxplot_1.png", width = 15, height = 7)

ggplot(data = Janesick_plot_long, 
       aes(x = variable, y = value, color = variable)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.4) +
  xlab("Pathology annotation") + 
  ylab("Proportion") +
  coord_flip() +
  facet_wrap(~Annotation, nrow=3, scales = "free")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        strip.text = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position="none", legend.box = "horizontal") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

ggsave("consistency_boxplot_2.png", width = 25, height = 10)



# We et al ----------------------------------------------------------------

# metadata_files = list.files(meta_dir, pattern="_metadata.csv")
# metadata = NULL
# for(f1 in metadata_files){
#   p1 <- read.csv(file.path(meta_dir, f1))
#   s1 <- gsub("_metadata.csv", "",f1)
#   p1$X = paste(s1, p1$X, sep="_")
#   metadata <- rbind(metadata, p1)
# }
dvn_prop_files = list.files(dvn_dir, pattern="_celltype_major.csv")
major_Wu_full = NULL
for(f1 in dvn_prop_files){
  p1 <- read.csv(file.path(dvn_dir, f1))
  s1 <- gsub("_celltype_major.csv", "",f1)
  s2 <- gsub("Proportion_","",s1)
  p1$X = paste(s2, p1$X, sep="_")
  major_Wu_full <- rbind(major_Wu_full, p1)
}
major_10X_FF <- major_10X_FF[, c("X", sort(colnames(major_10X_FF)[2:10]))]
major_10X_FFPE <- major_10X_FFPE[, c("X", sort(colnames(major_10X_FFPE)[2:10]))]
major_Wu_full <- major_Wu_full[, c("X", sort(colnames(major_Wu_full)[2:10]))]

major_prop_FF = major_10X_FF
rownames(major_prop_FF) = major_prop_FF$X
major_prop_FF = major_prop_FF[,-1]
major_prop_FF$max <- apply(major_prop_FF, 1, max)
major_prop_FF$wmax <- apply(major_prop_FF, 1, which.max)
major_prop_FF$sample <- "10X Genomics: Fresh Frozen"
colnames(major_prop_FF)[1:9] <- gsub("\\.", " ", colnames(major_prop_FF)[1:9])

major_prop_FFPE = major_10X_FFPE
rownames(major_prop_FFPE) = major_prop_FFPE$X
major_prop_FFPE = major_prop_FFPE[,-1]
major_prop_FFPE$max <- apply(major_prop_FFPE, 1, max)
major_prop_FFPE$wmax <- apply(major_prop_FFPE, 1, which.max)
major_prop_FFPE$sample <- "10X Genomics: FFPE"
colnames(major_prop_FFPE)[1:9] <- gsub("\\.", " ", colnames(major_prop_FFPE)[1:9])

major_prop_Wu_full = major_Wu_full
rownames(major_prop_Wu_full) = major_prop_Wu_full$X
major_prop_Wu_full = major_prop_Wu_full[,-1]
major_prop_Wu_full$max <- apply(major_prop_Wu_full, 1, max)
major_prop_Wu_full$wmax <- apply(major_prop_Wu_full[,1:9], 1, which.max)
major_prop_Wu_full$sample <- sapply(str_split(rownames(major_prop_Wu_full), "_"), "[[", 1)
major_prop_Wu_full$sample <- paste0("Wu et al.: ", major_prop_Wu_full$sample)
colnames(major_prop_Wu_full)[1:9] <- gsub("\\.", " ", colnames(major_prop_Wu_full)[1:9])

number_max_prop <- data.frame(celltype = colnames(major_prop_FF)[1:9])
t1 <- as.data.frame(table(colnames(major_prop_FF)[1:9][major_prop_FF$wmax]))
number_max_prop <- merge(number_max_prop, t1, by.x = "celltype", by.y = "Var1", all = T)
t2 <- as.data.frame(table(colnames(major_prop_FFPE)[1:9][major_prop_FFPE$wmax]))
number_max_prop <- merge(number_max_prop, t2, by.x = "celltype", by.y = "Var1", all = T)
t3 <- as.data.frame.matrix(table(colnames(major_prop_Wu_full)[1:9][major_prop_Wu_full$wmax], major_prop_Wu_full$sample))
t3$Var1 <- rownames(t3)
number_max_prop <- merge(number_max_prop, t3, by.x = "celltype", by.y = "Var1", all = T)
number_max_prop[is.na(number_max_prop)] <- 0
write.csv(number_max_prop, "number_of_max_celltype.csv")

major_prop_all <- bind_rows(major_prop_FF, major_prop_FFPE, major_prop_Wu_full)
major_prop_all$max.celltype <- colnames(major_prop_all)[1:9][major_prop_all$wmax]


major_prop_all %>% 
  group_by(sample) %>%
  mutate(med = median(as.numeric(max))) %>%
  ggplot(aes(x = max, fill = max.celltype)) +
  geom_histogram(bins = 100, boundary = 0, closed = "left") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  geom_vline(aes(xintercept= med), colour='red', linetype = "dashed") +
  facet_wrap(~sample, scales = "free", nrow = 2) +
  scale_fill_manual("Cell type",
                    values=colors) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm")) +
  labs(x = "Patch-wise Maximum Proportion", y = "Count") 
ggsave("hist_max_proportions.png", width = 10, height = 4)

celltypes <- colnames(major_prop_all)[1:9]
major_plot <- major_prop_all
major_plot$X <- rownames(major_plot)

major_long <- reshape(major_plot, direction =  "long", 
                      varying = celltypes,
                      timevar = "celltype", 
                      times = celltypes,
                      v.names = c("proportion"), 
                      idvar = colnames(major_plot)[12:14])
rownames(major_long) <- NULL
major_long$celltype <- factor(major_long$celltype, levels = celltypes)
samples <- sort(unique(major_plot$sample))
max_celltypes <- c("CAFs", "Cancer Epithelial", "Cancer Epithelial", "Plasmablasts", 
                   "Cancer Epithelial", "CAFs", "B cells", "Cancer Epithelial")
p_list = list()
for(i in 1:length(samples)){
  major_tmp <- major_plot[major_plot$sample == samples[i], ]
  major_tmp <- major_tmp %>%
    arrange(across(max_celltypes[i]))
  X_order <- unlist(major_tmp$X)
  
  tmp <- major_long[major_long$sample == samples[i], ]
  tmp$X <- factor(tmp$X, levels=X_order)
  tmp$celltype <- relevel(tmp$celltype, max_celltypes[i])
  tmp$celltype <- factor(tmp$celltype, levels=rev(levels(tmp$celltype)))
  
  p <- ggplot(tmp, aes(fill=factor(celltype), color=factor(celltype), y=proportion, x=X)) + 
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 9),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14)) +
    scale_fill_manual("Cell type",
                      values=colors) +
    scale_color_manual("Cell type",
                       values=colors) +
    labs(title = samples[i],
         x = "Patches",
         y = "Observed Proportion") +
    guides(fill=guide_legend(nrow=1),
           color=guide_legend(nrow = 1))
  
  mylegend <- g_legend(p)
  p_list[[i]] <- p
} 

p <- grid.arrange(arrangeGrob(p_list[[1]] + theme(legend.position="none"),
                              p_list[[2]] + theme(legend.position="none"),
                              p_list[[3]] + theme(legend.position="none"),
                              p_list[[4]] + theme(legend.position="none"),
                              p_list[[5]] + theme(legend.position="none"),
                              p_list[[6]] + theme(legend.position="none"),
                              p_list[[7]] + theme(legend.position="none"),
                              p_list[[8]] + theme(legend.position="none"), 
                              nrow = 4), mylegend, heights=c(8, 1))

ggsave("barplot_prop.png",p,  width = 16, height = 10)

ggplot(major_long, 
       aes(x=celltype, y=proportion, color=celltype)) + 
  geom_boxplot(width = 0.5, outlier.alpha = 0.3, outlier.size = 0.8) + 
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual("Cell type", values = colors) +
  scale_color_manual("Cell type", values = colors) +
  labs(y = "Proportion", x = "Cell Type") +
  facet_wrap(~sample, nrow = 2) +
  guides(color=guide_legend(nrow=1))


ggsave("boxplot_proportions.png", width = 10, height = 3.8)


########### Collapsed cell types #######
major_patho_FF <- data.frame(sample = "10X Genomics: Fresh Frozen",
                             X = major_10X_FF$X,
                             `Invasive Cancer` = major_10X_FF$Cancer.Epithelial,
                             Stroma = major_10X_FF$Endothelial + major_10X_FF$PVL + major_10X_FF$CAFs,
                          Lymphocyte = major_10X_FF$T.cells + major_10X_FF$B.cells + major_10X_FF$Plasmablasts,
                          Others = major_10X_FF$Myeloid + major_10X_FF$Normal.Epithelial)

major_patho_FFPE <- data.frame(sample = "10X Genomics: FFPE",
                               X = major_10X_FFPE$X,
                               `Invasive Cancer` = major_10X_FFPE$Cancer.Epithelial,
                             Stroma = major_10X_FFPE$Endothelial + major_10X_FFPE$PVL + major_10X_FFPE$CAFs,
                             Lymphocyte = major_10X_FFPE$T.cells + major_10X_FFPE$B.cells + major_10X_FFPE$Plasmablasts,
                             Others = major_10X_FFPE$Myeloid + major_10X_FFPE$Normal.Epithelial)

major_patho_Wu <- data.frame(sample = paste0("Wu et al.: ", sapply(str_split(rownames(major_prop_Wu_full), "_"), "[[", 1)),
                             X = major_Wu_full$X,
                             `Invasive Cancer` = major_Wu_full$Cancer.Epithelial,
                               Stroma = major_Wu_full$Endothelial + major_Wu_full$PVL + major_Wu_full$CAFs,
                               Lymphocyte = major_Wu_full$T.cells + major_Wu_full$B.cells + major_Wu_full$Plasmablasts,
                               Others = major_Wu_full$Myeloid + major_Wu_full$Normal.Epithelial)
major_patho <- bind_rows(major_patho_FF, major_patho_FFPE, major_patho_Wu)
colnames(major_patho)[3:6] <- c("Invasive Cancer", "Stroma", "Lymphocyte", "Others")
major_patho$max <- apply(major_patho[, 3:6], 1, max)
major_patho$wmax <- apply(major_patho[, 3:6], 1, which.max)
major_patho$max.celltype <- colnames(major_patho)[3:6][major_patho$wmax]

number_max_patho <- data.frame(celltype = colnames(major_patho)[3:6])
number_max_patho <- as.data.frame(table(colnames(major_patho)[3:6][major_patho$wmax], major_patho$sample))
number_max_patho <- number_max_patho %>%
  pivot_wider(names_from = Var2, values_from = Freq)
write.csv(number_max_patho, "number_of_max_celltype_patho.csv")
major_patho %>% 
  group_by(sample) %>%
  mutate(med = median(as.numeric(max))) %>%
  ggplot(aes(x = max, fill = max.celltype)) +
  geom_histogram(bins = 100, boundary = 0, closed = "left") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2)) +
  geom_vline(aes(xintercept= med), colour='red', linetype = "dashed") +
  facet_wrap(~sample, scales = "free", nrow = 2) +
  scale_fill_manual("Collapsed Cell type",
                    values=colors_patho) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.5, "cm")) +
  labs(x = "Patch-wise Maximum Proportion", y = "Count")

ggsave("hist_max_proportions_patho.png", width = 10, height = 3.8)


celltypes <- colnames(major_patho)[3:6]
major_plot <- major_patho
major_long <- reshape(major_plot, direction =  "long", 
                      varying = celltypes,
                      timevar = "celltype", 
                      times = celltypes,
                      v.names = c("proportion"), 
                      idvar = colnames(major_plot)[1:2])
rownames(major_long) <- NULL
major_long$celltype <- factor(major_long$celltype, levels = celltypes)
samples <- sort(unique(major_plot$sample))
max_celltypes <- c("Others", "Invasive Cancer", "Lymphocyte", "Lymphocyte", 
                   "Invasive Cancer", "Stroma", "Lymphocyte", "Invasive Cancer")
p_list = list()
for(i in 1:length(samples)){
  major_tmp <- major_plot[major_plot$sample == samples[i], ]
  major_tmp <- major_tmp %>%
    arrange(across(max_celltypes[i]))
  X_order <- unlist(major_tmp$X)
  
  tmp <- major_long[major_long$sample == samples[i], ]
  tmp$X <- factor(tmp$X, levels=X_order)
  tmp$celltype <- relevel(tmp$celltype, max_celltypes[i])
  tmp$celltype <- factor(tmp$celltype, levels=rev(levels(tmp$celltype)))
  
  p <- ggplot(tmp, aes(fill=factor(celltype), color=factor(celltype), y=proportion, x=X)) + 
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(plot.title = element_text(size = 13),
          axis.title = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 9),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14)) +
    scale_fill_manual("Collapsed Cell type",
                      values=colors_patho) +
    scale_color_manual("Collapsed Cell type",
                       values=colors_patho) +
    labs(title = samples[i],
         x = "Patches",
         y = "Observed Proportion")
  mylegend <- g_legend(p)
  p_list[[i]] <- p
} 

p <- grid.arrange(arrangeGrob(p_list[[1]] + theme(legend.position="none"),
                              p_list[[2]] + theme(legend.position="none"),
                              p_list[[3]] + theme(legend.position="none"),
                              p_list[[4]] + theme(legend.position="none"),
                              p_list[[5]] + theme(legend.position="none"),
                              p_list[[6]] + theme(legend.position="none"),
                              p_list[[7]] + theme(legend.position="none"),
                              p_list[[8]] + theme(legend.position="none"), 
                              nrow = 4), mylegend, heights=c(8, 1))
ggsave("barplot_prop_patho.png", p, width = 16, height = 10)


ggplot(major_long, 
       aes(x=celltype, y=proportion, color=celltype)) + 
  geom_boxplot(width = 0.5, outlier.alpha = 0.3, outlier.size = 0.8) + 
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual("Collapsed Cell type", values = colors_patho) +
  scale_color_manual("Collapsed Cell type", values = colors_patho) +
  labs(y = "Proportion", x = "Collapsed Cell Type") +
  facet_wrap(~sample, nrow = 2)

ggsave("boxplot_proportions_patho.png", width = 10, height = 3.8)

########## Consistency #########
meta_dir = "../../st2image_data/Wu_2021/data/metadata/"
metadata_files = list.files(meta_dir, pattern="_metadata.csv")
metadata = NULL
for(f1 in metadata_files){
  p1 <- read.csv(file.path(meta_dir, f1))
  s1 <- gsub("_metadata.csv", "",f1)
  p1$X = paste(s1, p1$X, sep="_")
  metadata <- rbind(metadata, p1)
}
metadata$Classification = ifelse(metadata$Classification == "", NA, metadata$Classification)
colnames(major_Wu_full)[2:10] <- gsub("\\.", " ", colnames(major_Wu_full)[2:10])

major_plot <- merge(major_Wu_full, metadata[, c("X", "patientid", "Classification")], by = "X")
#groups <- droplevels(counts_1$Var1[1:10])
#major_full_plot <- major_full_plot[major_full_plot$Classification %in% groups,]

major_plot_long <- melt(major_plot, id.vars=c('Classification', 'patientid'), 
                        measure.vars=colnames(major_plot)[2:10])

ggplot(data = major_plot_long, 
       aes(x = Classification, y = value, color = Classification)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.5) + 
  theme(legend.position="none",
        axis.title = element_text(size = 8),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        strip.text = element_text(size=8)) + 
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_grid(patientid~variable, scales = "free")
ggsave("consistency_boxplot_1.png", width = 13, height = 7)

ggplot(data = major_plot_long, 
       aes(x = Classification, y = value, color = variable)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.4) +
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_wrap(~patientid, nrow=3, scales = "free")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        strip.text = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position="bottom", legend.box = "horizontal") +
  scale_color_manual("Cell Type", values = colors) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))

ggsave("consistency_boxplot_2.png", width = 16, height = 25)


### Collapsed. 
colnames(major_patho_Wu) <- gsub("\\.", " ", colnames(major_patho_Wu))

major_patho_plot <- merge(major_patho_Wu, metadata[, c("X", "patientid", "Classification")], by = "X")
#groups <- droplevels(counts_1$Var1[1:10])
#major_full_plot <- major_full_plot[major_full_plot$Classification %in% groups,]

major_patho_plot_long <- melt(major_patho_plot, id.vars=c('Classification', 'patientid', 'X'), 
                        measure.vars=colnames(major_patho_plot)[3:6])

ggplot(data = major_patho_plot_long, 
       aes(x = Classification, y = value, color = Classification)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.5) + 
  theme(legend.position="none",
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        strip.text = element_text(size=10)) + 
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_grid(patientid~variable, scales = "free")

ggsave("consistency_boxplot_patho_1.png", width = 13, height = 7.5)

ggplot(data = major_patho_plot_long, 
       aes(x = Classification, y = value, color = variable)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.4, width = 0.7) +
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_wrap(~patientid, nrow=2, scales = "free")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        strip.text = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position="bottom", legend.box = "horizontal") +
  scale_color_manual("Collapsed Cell Type", values = colors_patho) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))

ggsave("consistency_boxplot_patho_2.png", width = 16, height = 13)


exclude_annotation <- c("NA", "Uncertain", "Artefact")
major_patho_plot_long_red <- major_patho_plot_long[!(major_patho_plot_long$Classification %in% exclude_annotation), ]
major_patho_plot_long_red <- na.omit(major_patho_plot_long_red)

ggplot(data = major_patho_plot_long_red, 
       aes(x = Classification, y = value, color = Classification)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.5) + 
  theme(legend.position="none",
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        strip.text = element_text(size=10)) + 
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  ggh4x::facet_grid2(patientid~variable, scales = "free_y")

ggsave("consistency_boxplot_patho_reduced_1.png", width = 13, height = 7.5)


samples <- sort(unique(major_patho_plot_long_red$patientid))
p_list = list()
for(i in 1:length(samples)){
  p <- ggplot(data = major_patho_plot_long_red[major_patho_plot_long_red$patientid == samples[i], ], 
              aes(x = X, y = value, color = variable, fill = variable)) + 
    geom_bar(position="stack", stat="identity", width = 1) +
    theme(axis.title = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 9),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11),
          strip.text.x = element_text(size = 11)) +
    scale_fill_manual("Collapsed Cell type",
                      values=colors_patho) +
    scale_color_manual("Collapsed Cell type",
                       values=colors_patho) +
    labs(title = samples[i],
         x = "Patches",
         y = "Observed Proportion") +
    guides(fill=guide_legend(nrow=1),
           color=guide_legend(nrow = 1)) +
    facet_wrap(~Classification, scales = "free", ncol = 3)
  
  mylegend <- g_legend(p)
  p_list[[i]] <- p
} 

p <- grid.arrange(arrangeGrob(p_list[[1]] + theme(legend.position="none"),
                              p_list[[2]] + theme(legend.position="none"),
                              p_list[[3]] + theme(legend.position="none"),
                              p_list[[4]] + theme(legend.position="none"),
                              p_list[[5]] + theme(legend.position="none"),
                              p_list[[6]] + theme(legend.position="none"),
                              nrow = 3), mylegend, heights=c(8, 1))
p_list[[6]]

ggsave("boxplot_per_annotation_6.png",p_list[[6]],  width = 10, height = 4)

counts <- as.data.frame.matrix(table(major_patho_plot_long$Classification, major_patho_plot_long$patientid))
counts <- counts != 0
indicator <- rowSums(counts)
include_annotation <- names(indicator[which(indicator != 1)])
major_patho_plot_long$Classification <- gsub(" \\+", "\n\\+", major_patho_plot_long$Classification)
major_patho_plot_long$Classification <- gsub("Cancer trapped in lymphocyte aggregation", "Cancer trapped\nin lymphocyte\naggregation", major_patho_plot_long$Classification)

ggplot(data = major_patho_plot_long[major_patho_plot_long$Classification %in% include_annotation, ], 
       aes(x = patientid, y = value, color = variable)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.4, width = 0.5) +
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_wrap(~Classification, nrow=2, scales = "free")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 9, color = "black"), 
        axis.text.y = element_text(size = 9, color = "black"),
        strip.text = element_text(size=9),
        legend.text = element_text(size=9),
        legend.position="bottom", 
        legend.box = "horizontal") +
  scale_color_manual("Collapsed Cell Type", values = colors_patho) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_y_continuous(limits = c(0,1))
ggsave("consistency_across_sample.png", width = 10, height = 8)

  
major_patho_plot_long_red$Classification <- gsub(" \\+", "\n\\+", major_patho_plot_long_red$Classification)
major_patho_plot_long_red$Classification <- gsub("Cancer trapped in lymphocyte aggregation", "Cancer trapped\nin lymphocyte\naggregation", major_patho_plot_long_red$Classification)

ggplot(data = major_patho_plot_long_red, 
       aes(x = Classification, y = value, color = variable)) + 
  geom_boxplot(outlier.alpha = .5, outlier.size = 0.4, width = 0.5) +
  xlab("Pathology annotation") + 
  ylab("Proportion") + 
  coord_flip() +
  facet_wrap(~patientid, nrow=2, scales = "free")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"), 
        axis.text.y = element_text(size = 15, color = "black"),
        strip.text = element_text(size=17),
        legend.text = element_text(size=15),
        legend.position="bottom", legend.box = "horizontal") +
  scale_color_manual("Collapsed Cell Type", values = colors_patho) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_y_continuous(limits = c(0,1))+
  scale_x_discrete(expand=c(0.15,0))
ggsave("consistency_boxplot_patho_reduced_2.png", width = 16, height = 11)


