# some population structure based sanity checks on masked_regions
library(tidyverse)
library(adegenet)
library(FactoMineR)
library(ggrepel)

# get a vector of populations in the order of samples
get_loc_from_samples <- function(sample_list){
  loc <- rep(NA, length(sample_list))
  loc[grep("ba", sample_list)] <- "Southern Eurasian"
  loc[grep("ca", sample_list)] <- "Southern Eurasian"
  loc[grep("cr", sample_list)] <- "Western Eurasian"
  loc[grep("ka", sample_list)] <- "Eastern Eurasian"
  loc[grep("ki", sample_list)] <- "Western Eurasian"
  loc[grep("la", sample_list)] <- "Western Eurasian"
  loc[grep("no", sample_list)] <- "Western Eurasian"
  loc[grep("po", sample_list)] <- "Western Eurasian"
  loc[grep("og", sample_list)] <- "Eastern Eurasian"
  loc[grep("to", sample_list)] <- "Eastern Eurasian"
  loc[grep("tu", sample_list)] <- "Eastern Eurasian"
  loc[grep("ur", sample_list)] <- "Western Eurasian"
  loc[grep("vl", sample_list)] <- "Eastern Eurasian"
  loc[grep("ya", sample_list)] <- "Eastern Eurasian"
  loc[grep("sm", sample_list)] <- "Iberian"
  return(loc)
}

get_color_from_samples <- function(sample_list){
  # get a vector of colors in the order of samples
  col <- rep(NA, length(sample_list))
  col[grep("ba", sample_list)] <- "#A035AF"
  col[grep("ca", sample_list)] <- "#B8860b"
  col[grep("cr", sample_list)] <- "#CAB2D6"
  col[grep("ka", sample_list)] <- "#FDBF6F"
  col[grep("ki", sample_list)] <- "#440154FF"
  col[grep("la", sample_list)] <- "#B2DF8A"
  col[grep("no", sample_list)] <- "#3B528BFF"
  col[grep("po", sample_list)] <- "#21908CFF"
  col[grep("og", sample_list)] <- "#FDBF6F"
  col[grep("to", sample_list)] <- "#FDBF6F"
  col[grep("tu", sample_list)] <- "#FF7F00"
  col[grep("ur", sample_list)] <- "#0F4909"
  col[grep("vl", sample_list)] <- "#FB9A99"
  col[grep("ya", sample_list)] <- "#E31A1C"
  col[grep("sm", sample_list)] <- "#FFFF77"
  return(col)
}

# read plink raw data into a data.frame of the genotypes
get_gt_df_from_raw <- function(raw_table){
  gl <- read.PLINK(raw_table)
  raw_df <- data.frame(as.matrix(gl))
  return(raw_df)
}

# build a data.frame of samples (PCs, populations, colors)
# from a PCA object
build_df_from_pca <- function(pca_obj){
  pca_df <- as.data.frame(pca_obj$ind$coord)
  pca_df$population <- get_loc_from_samples(row.names(pca_df))
  pca_df$color <- get_color_from_samples(row.names(pca_df))
  return(pca_df)
}

# produce a pca biplot from pca data.frame and variance percentages
plot_pca <- function(pca_df, variance_percents, x = 1, y = 2){
  pca_plot <- ggplot(data = pca_df, 
                     aes(x = pca_df[,x], y = pca_df[,y],
                         label=rownames(pca_df))
                     ) +
    geom_point(fill = pca_df[,"color"], shape = 21, size = 2.5) +
    theme_bw() +
    theme(panel.background = element_blank(),
          # panel.grid = element_blank(),
          plot.background = element_blank(),
          #axis.text.x = element_blank(),
          #axis.ticks.x=element_blank()
    ) +
    # xlim(min = min_lim, max = max_lim) +
    # ylim(min = min_lim, max = max_lim) +
    
    xlab(paste0("PC", x, " - ", round(variance_percents[x], 2), "%")) +
    ylab(paste0("PC", y, " - ", round(variance_percents[y], 2), "%")) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9)) + geom_label_repel(max.overlaps = 100)
  
  return(pca_plot)
}

## LOAD DATA ##

# load masked_regions df
gt_df <- get_gt_df_from_raw("data/demographic_inference/lpa-eel.masked_regions_only.raw")

# Remove SNPs with missing data
all_df <- gt_df[, colSums(is.na(gt_df)) == 0]

# list of lpa individuals and lpa only df
lpa_samples <- rownames(all_df)[grep("lp", rownames(all_df))]
lpa_df <- all_df[rownames(all_df) %in% lpa_samples, ]
# filter out non-variants
different_columns <- sapply(lpa_df, function(x) length(unique(x)) > 1)
lpa_df <- lpa_df[, different_columns]

# list of eel individuals and eel only df
eel_samples <- rownames(all_df)[grep("ll", rownames(all_df))]
eel_df <- all_df[rownames(all_df) %in% eel_samples, ]
# filter out non-variants
different_columns <- sapply(eel_df, function(x) length(unique(x)) > 1)
eel_df <- eel_df[, different_columns]

## PCA AND PLOT ##

# all samples - pca
all_pca <- FactoMineR::PCA(all_df, graph = F)
# all samples - pca dataframe
all_pca_df <- build_df_from_pca(pca_obj = all_pca)
# all samples - variance percents
all_pca_variance_percents <- all_pca$eig[,2]
# all samples - pca plot
all_pca_plot <- plot_pca(pca_df = all_pca_df,
                         variance_percents = all_pca_variance_percents)

ggsave(filename = paste0("plots/demographic_inference/sanity/lpa-eel.all_pca.pdf"),
       plot = all_pca_plot, width = 300, height = 160, units = "mm")


# lpa samples - pca
lpa_pca <- FactoMineR::PCA(lpa_df, graph = F)
# lpa samples - pca dataframe
lpa_pca_df <- build_df_from_pca(pca_obj = lpa_pca)
# lpa samples - variance percents
lpa_pca_variance_percents <- lpa_pca$eig[,2]
# lpa samples - pca plot
lpa_pca_plot <- plot_pca(pca_df = lpa_pca_df,
                         variance_percents = lpa_pca_variance_percents)

ggsave(filename = paste0("plots/demographic_inference/sanity/lpa-eel.lpa_pca.pdf"),
       plot = lpa_pca_plot, width = 300, height = 160, units = "mm")

# eel samples - pca
eel_pca <- FactoMineR::PCA(eel_df, graph = F)
# eel samples - pca dataframe
eel_pca_df <- build_df_from_pca(pca_obj = eel_pca)
# eel samples - variance percents
eel_pca_variance_percents <- eel_pca$eig[,2]
# eel samples - pca plot
eel_pca_plot <- plot_pca(pca_df = eel_pca_df,
                         variance_percents = eel_pca_variance_percents)

ggsave(filename = paste0("plots/demographic_inference/sanity/lpa-eel.eel_pca.pdf"),
       plot = eel_pca_plot, width = 300, height = 160, units = "mm")

