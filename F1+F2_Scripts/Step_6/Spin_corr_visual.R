library(ggcorrplot)

spin_corr <- read.csv('Spin_corr_mat_BASELINE-FULL_3x3.csv') #read in correlation matrix
names <- spin_corr[,1] #set variable names
spin_corr <- spin_corr[-c(1)] #remove variable name column
rownames(spin_corr) <- names
colnames(spin_corr) <- names

ggcorrplot(spin_corr, show.diag = FALSE, hc.order = FALSE, outline.color = "white", lab = TRUE)
ggsave("Spin_Results/Spin_corrs_BASELINE-FULL_3x3.tiff", width = 9, height = 7.5, dpi = 600, units = "cm")