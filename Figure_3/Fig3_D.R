library(tidyverse)
library(circlize)
library(ComplexHeatmap)



### Heatmap of VDJ usage

df = read.csv("Data/TCR_usage.csv")

df$Genomic_subtype = factor(df$Genomic_subtype, levels = c("HOXA", "LMO2_LYL1", "TLX1", "TLX3", "TAL1", "T-other"))
df$Flow_subtype = factor(df$Flow_subtype, levels = c("ETP", "nonETP"))
df$Response = factor(df$Response, levels = c("Responsive", "nonResponsive"))
df$ZBTB16_status = factor(df$ZBTB16_status, levels = c("ZBTB16_neg", "ZBTB16_pos"))

df = df %>% arrange(ZBTB16_status, TRB, TRA, TRD, TRG, Genomic_subtype, Flow_subtype)
mat = as.matrix(df[, c("TRD", "TRG", "TRB", "TRA")])
mat = t(mat)
colnames(mat) = df$Patient_ID

bottom_ha = HeatmapAnnotation(
  which = "column",
  'Response' = df$Response,
  col = list(
    'Response' = c('Responsive' = '#1D71B8', 'nonResponsive' = '#EA0017')
  ),
  gp = gpar(col = "black"),
  gap = unit(1, "mm"),
  annotation_name_side = "right",
  annotation_name_gp= gpar(fontsize = 10),
  annotation_legend_param = list(border = "black")
)

top_ha = HeatmapAnnotation(
  which = "column",
  'Genomic_subtype' = anno_text(df$Genomic_subtype, location = 0, just = "left", gp = gpar(fontsize = 8)),
  'Flow_subtype' = df$Flow_subtype,
  col = list('Flow_subtype' = c('nonETP' = '#FFFFFF', 'ETP' = 'grey30')),
  gp = gpar(col = "black"),
  gap = unit(1, "mm"),
  annotation_name_side = "right",
  annotation_name_gp= gpar(fontsize = 10),
  annotation_legend_param = list(border = "black")
)

hm = Heatmap(
  mat,
  col = c('1' = 'grey30', '0' = '#FFFFFF'),
  name = "TCR expressed",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  row_names_side = "right",
  column_split = df$ZBTB16_status,
  column_title = NULL,
  rect_gp = gpar(col = "black"),
  bottom_annotation = bottom_ha,
  top_annotation = top_ha,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(at = c(1, 0), labels = c('Yes', 'No'), border = "black")
)
pdf("Plots/Fig3D_VDJ_usage.pdf", width = 5.8, height = 2.4)
hm = draw(
  hm,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
)
dev.off()


