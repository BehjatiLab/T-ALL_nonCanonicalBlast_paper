library(tidyverse)
library(ggpubr)
library(ggrastr)



### Phased BAF for P058_D0_WGS and P058_D28_WGS

# Load data from P058_D0_WGS_BAF and P058_D28_WGS_BAF
df_D0 = read.csv("Data/P058_D0_WGS_BAF.csv")
df_D28 = read.csv("Data/P058_D28_WGS_BAF.csv")

# Transfer phasing from P058_D0_WGS to P058_D28_WGS
match_index = match(x = df_D28$CHR_POS_REF_ALT, table = df_D0$CHR_POS_REF_ALT)
df_D28$altIsMum = df_D0$altIsMum[match_index]

list_of_df = list(df_D0, df_D28)
list_of_sample_ID = c("P058_D0_WGS", "P058_D28_WGS")

for (i in seq_along(list_of_df)) {
  
  list_of_CHR = c('chr9', 'chr17')
  list_of_breakpoints = c(36706115, 57793447)
  
  for (j in seq_along(list_of_CHR)) {
    
    df = list_of_df[[i]][list_of_df[[i]]$CHR == list_of_CHR[j], ]
    
    # Plot WGS BAF
    p = ggplot(df, aes(x = POS, y = WGS_BAF)) +
      rasterise(geom_point(aes(color = altIsMum), size = 0.1), dpi = 1200) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = list_of_breakpoints[j], linetype = "dashed") +
      scale_x_continuous(breaks = c(list_of_breakpoints[j]), labels = round(c(list_of_breakpoints[j])/1e6, digits = 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
      scale_color_manual(values = c("TRUE" = "#F77F11", "FALSE" = "#000000"), na.value = "#D3D3D3") +
      ylab("BAF") +
      xlab("Genomic position (Mb)") +
      theme_bw() +
      theme(
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    ggexport(p, filename = paste0("Plots/Fig4C_WGS_BAFphased_", list_of_sample_ID[i], "_", list_of_CHR[j], ".pdf"), width = 4, height = 2.6)
  }
  
}



### scRNA BAF plots for P058_D0 and P058_D28

list_of_sample_ID = c('P058_D0', 'P058_D28')

for (i in seq_along(list_of_sample_ID)) {

  # Load scRNA BAF (all chromosomes)
  df_all = read.csv(paste0("Data/", list_of_sample_ID[i], "_scRNA_BAF.csv"))
  
  # Subset for leukaemia cells and SNPs with at least 10 reads
  df_all = df_all[df_all$cellID == "Leuk" & df_all$totCount >= 10, ]
  
  # Calcuate BAF of each SNP using counts aggregated across cells (scRNA_BAF)
  df_all$scRNA_BAF = df_all$altCount / df_all$totCount
  
  # Transfer phasing from P058_D0_WGS
  df_all$CHR_POS_REF_ALT = paste(df_all$CHR, df_all$POS, df_all$REF, df_all$ALT, sep = "_")
  match_index = match(x = df_all$CHR_POS_REF_ALT, table = df_D0$CHR_POS_REF_ALT)
  df_all$altIsMum = df_D0$altIsMum[match_index]
  
  list_of_CHR = c('chr9', 'chr17')
  list_of_breakpoints = c(36706115, 57793447)
  
  for (j in seq_along(list_of_CHR)) {
    
    df = df_all[df_all$CHR == list_of_CHR[j], ]
    
    start_chr = min(df$POS)
    end_chr = max(df$POS)
    
    # Plot scRNA BAF
    p = ggplot(df, aes(x = POS, y = scRNA_BAF)) +
      geom_point(aes(color = altIsMum), size = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = list_of_breakpoints[j], linetype = "dashed") +
      scale_x_continuous(breaks = c(list_of_breakpoints[j]), labels = round(c(list_of_breakpoints[j])/1e6, digits = 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"), limits = c(0, 1)) +
      scale_color_manual(values = c("TRUE" = "#F77F11", "FALSE" = "#000000"), na.value = "#D3D3D3") +
      ylab("BAF") +
      xlab("Genomic position (Mb)") +
      theme_bw() +
      theme(
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    ggexport(p, filename = paste0("Plots/Fig4D_scRNA_BAFphased_", list_of_sample_ID[i], "_", list_of_CHR[j], ".pdf"), width = 4, height = 2.6)
  }
  
}