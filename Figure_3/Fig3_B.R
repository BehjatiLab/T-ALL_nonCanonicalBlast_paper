library(tidyverse)



### Figure 3B: Normal-to-leukaemia transcriptome comparison by logistic regression

# Load data
similarity_matrix = read.csv(file = paste0("Data/LogisticRegression.csv"), row.names = 1, check.names = FALSE)
similarity_matrix = as.matrix(similarity_matrix)

# Ordering of rows and columns
celltypes_order = c('ZBTB16_neg', 'ZBTB16_pos', 'SCP')
cell_groupings_order = c('HSC_MPP', 'CYCLING_MPP', 'LMPP_MLP', 'DN(early)_T', 'DN(P)_T', 'DN(Q)_T', 'DP(P)_T', 'DP(Q)_T', 'ABT(ENTRY)', 'CD8+T', 'CD4+T', 'TREG', 'CYCLING_T', 'CD8AA', 'TYPE_1_INNATE_T', 'TYPE_3_INNATE_T', 'ILC2', 'ILC3', 'CYCLING_ILC', 'NK', 'CYCLING_NK', 'SCP')

# Extract cell groupings
cell_groupings = rownames(similarity_matrix) %>% sub(pattern = "::.*", replacement = "")

# Convert to discrete values
mtx = similarity_matrix
mtx_discrete = matrix(nrow = nrow(mtx), ncol = ncol(mtx))
mtx_discrete[mtx < 0.3] = "a"
mtx_discrete[mtx >= 0.3 & mtx < 0.4] = "b"
mtx_discrete[mtx >= 0.4 & mtx <= 0.6] = "c"
mtx_discrete[mtx > 0.6 & mtx <= 0.7] = "d"
mtx_discrete[mtx > 0.7] = "e"
rm(mtx)

# Add column names and convert mtx_discrete into a data frame
colnames(mtx_discrete) = colnames(similarity_matrix)
mtx_discrete = as.data.frame(mtx_discrete)

# Add cell_groupings
mtx_discrete$cell_groupings = cell_groupings

# Pivot longer
mtx_discrete = mtx_discrete %>%
  pivot_longer(-c(cell_groupings), names_to = "ref_celltypes", values_to = "similarity")

# Set order for similarity, ref_celltypes and cell_groupings
mtx_discrete = mtx_discrete %>%
  mutate(similarity = factor(similarity, levels = c("e", "d", "c", "b", "a")),
         ref_celltypes = factor(ref_celltypes, levels = celltypes_order),
         cell_groupings = factor(cell_groupings, levels = cell_groupings_order))

# Plot stacked barplot
p = ggplot(mtx_discrete, aes(y = 0, fill = similarity)) +
  geom_bar(stat = "count", position = "fill", color = "black", width = 0.1) +
  scale_x_reverse() +
  scale_fill_manual(values = c("a" = "#7F7F7F", "b" = "#B3B3B3", "c" = "#FFFFFF", "d" = "#77BDAD", "e" = "#309E88"),
                    labels = c("a"= "0.0 - <0.3", "b" = "0.3 - <0.4", "c" = "0.4 - 0.6", "d" = ">0.6 - 0.7", "e" = ">0.7 - 1.0")) +
  facet_grid(cell_groupings ~ ref_celltypes, switch = "x") +
  theme_void() +
  theme(strip.text.x = element_text(margin = margin(5, 5, 5, 5), angle = 90, vjust = 0.5, hjust = 1, size = 10),
        strip.text.y = element_text(margin = margin(5, 5, 5, 5), vjust = 0.5, hjust = 0, size = 10),
        legend.position = "bottom",
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"))
ggsave("Plots/Fig3B_LogisticRegression.pdf", p, height = 8, width = 8)