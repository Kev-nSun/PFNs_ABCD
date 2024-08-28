library(ggcorrplot)

corr <- read.csv('F1_F2_GP_corrs_BASELINE.csv') #read in corrs
names <- corr[,1]
corr <- corr[-c(1)]
rownames(corr) <- names
colnames(corr) <- names

ggcorrplot(corr, show.diag = FALSE, hc.order = TRUE, outline.color = "white", lab = TRUE)
ggsave("Corrs_clus_sq_nodiag_BASELINE.tiff", dpi = 600, width = 8, height = 7, units = "cm")

ggcorrplot(corr, hc.order = TRUE, type = "lower", outline.color = "white", lab = TRUE)
ggsave("Corrs_clus_BASELINE.tiff", dpi = 600, width = 7, height = 6, units = "cm")

p_vals <- c(2.87e-49, 7.87e-20, 0.048)
p_vals_FDR <- p.adjust(p_vals,"fdr") #FDR corrected p vals
write.csv(p_vals_FDR, file = "F1_F2_pfac_corr_p_vals_FDR.csv") #save out FDR p vals

p_vals_LME <- c(4.41e-48, 7.42e-18, 0.060)
p_vals_LME_FDR <- p.adjust(p_vals_LME,"fdr") #FDR corrected p vals
write.csv(p_vals_LME_FDR, file = "F1_F2_pfac_LME_p_vals_FDR.csv") #save out FDR p vals 
