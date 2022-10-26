options(stringsAsFactors = FALSE)

###Load libraries
library(ggplot2)
library(cowplot)
library(stringr)

###Read data
Us <- readRDS("../../output/gtex_mr_mash_analysis/data_driven_matrices/output/fold_1.ted_unconstrained.rds")$U
U <- Us$tPCA$mat
#U <- U/max(U)
Uc <- cov2cor(U)

###Order tissues
brain_tissues <- rownames(Uc)[stringr::str_detect(rownames(Uc), "Brain_")]
brain_cereb_tissues <- brain_tissues[stringr::str_detect(brain_tissues, "Brain_Cer")]
brain_other_tissues <- brain_tissues[!(brain_tissues %in% brain_cereb_tissues)]
cell_tissues <- rownames(Uc)[stringr::str_detect(rownames(Uc), "Cells_")]
testis <- "Testis"
whole_blood <- "Whole_Blood"
muscle_skeletal <- "Muscle_Skeletal"
other_tissues <- rownames(Uc)[!(rownames(Uc) %in% c(brain_tissues,
                                                    cell_tissues,
                                                    testis,
                                                    whole_blood,
                                                    muscle_skeletal,
                                                    "Kidney_Cortex"))]
ordered_tissues <- c(brain_cereb_tissues, brain_other_tissues, other_tissues, muscle_skeletal, cell_tissues, whole_blood, testis)

###Prepare the data for plotting
Uc <- Uc[ordered_tissues, ordered_tissues]
Uc[upper.tri(Uc,diag = FALSE)] <- NA
r <- nrow(Uc)
colors <- c("#4575b4","#74add1","#abd9e9","#e0f3f8",
            "#fee090","#fdae61","#f46d43","#d73027")
dat <- data.frame(x = rep(1:r,each = r),
                  y = rep(1:r,times = r),
                  u = as.vector(Uc))
dat <- transform(dat,
                 x = factor(x,1:r),
                 y = factor(y,r:1),
                 u = cut(u,breaks = 8))
tissues <- rownames(Uc)
levels(dat$x) <- tissues
levels(dat$y) <- tissues[r:1]

###Plot
p <- ggplot(dat,aes(x = x,y = y,fill = u)) +
  geom_tile() +
  scale_fill_manual(values = colors, na.translate = FALSE) +
  labs(x = "",y = "") +
  theme_cowplot(font_size = 7) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  guides(fill=guide_legend(title=""))


ggsave("../../analysis/paper_figures/gtex_tPCA_heatmap_fold_1.pdf",p,
       height = 7.2,width = 8)


