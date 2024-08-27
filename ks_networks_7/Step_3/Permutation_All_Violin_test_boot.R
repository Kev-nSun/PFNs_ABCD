
library(R.matlab)
library(ggplot2)
library(visreg)
library(PupillometryR)

# Actual data for A and B
PB1_mat_A <- readMat('../../Step_2/2_Server_All/BASELINE/baseline_results_031324/all_network/General_PB1_acc_testA.mat')
PB1_mat_B <- readMat('../../Step_2/2_Server_All/BASELINE/baseline_results_031324/all_network/General_PB1_acc_testB.mat')
F1_mat_A <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_1_acc_testA.mat')
F1_mat_B <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_1_acc_testB.mat')
F2_mat_A <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_2_acc_testA.mat')
F2_mat_B <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_040524/all_network/PRS_2_acc_testB.mat')

# Load in actual A corr values
GP_A <- PB1_mat_A$PB.test
F1_A <- F1_mat_A$PB.test
F2_A <- F2_mat_A$PB.test
# Load in actual B corr values
GP_B <- PB1_mat_B$PB.test
F1_B <- F1_mat_B$PB.test
F2_B <- F2_mat_B$PB.test
# Create dataframe with column for actual pred acc
data = data.frame(Corr_Actual = c(as.numeric(matrix(GP_A, 1, 1000)), as.numeric(matrix(GP_B, 1, 1000)), as.numeric(matrix(F1_A, 1, 1000)), as.numeric(matrix(F1_B, 1, 1000)), as.numeric(matrix(F2_A, 1, 1000)), as.numeric(matrix(F2_B, 1, 1000))))

# Load in permutation null data for A and B
GP_null_A <- list()
GP_null_B <- list()
F1_null_A <- list()
F1_null_B <- list()
F2_null_A <- list()
F2_null_B <- list()

for (phenotype in c("PB1","PRS_1","PRS_2")) {
  for (grp in c('A','B')) {
    for (dir in c('A','B','C','D','E','F','G','H','I','J')) {
      if (phenotype == "PB1") {
        dirin <- '../../Step_2/2_Server_All/BASELINE/'
        GP_mat<- readMat(paste0(dirin,'results_null_1000_031424/results_matchedsamples_null_1000_031424', dir, '/all_network/General_PB1_acc_null_test', grp, '.mat'))
        GP_nulls<- t(GP_mat$PB.test)
        if (grp == 'A') {
          GP_null_A <- append(GP_null_A,GP_nulls)
        } else {
          GP_null_B <- append(GP_null_B,GP_nulls)
        }
      } else if (phenotype == "PRS_1") {
        dirin <- '../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/'
        F1_mat<- readMat(paste0(dirin,'results_null_1000_040824/results_matchedsamples_null_1000_040824', dir, '/PRS_1_acc_null_test', grp, '.mat'))
        F1_nulls<- t(F1_mat$PB.test)
        if (grp == 'A') {
          F1_null_A <- append(F1_null_A,F1_nulls)
        } else {
          F1_null_B <- append(F1_null_B,F1_nulls)
        }
      } else {
        dirin <- '../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/'
        F2_mat<- readMat(paste0(dirin,'results_null_1000_040824/results_matchedsamples_null_1000_040824', dir, '/PRS_2_acc_null_test', grp, '.mat'))
        F2_nulls<- t(F2_mat$PB.test)
        if (grp == 'A') {
          F2_null_A <- append(F2_null_A,F2_nulls)
        } else {
          F2_null_B <- append(F2_null_B,F2_nulls)
        }
      }
    }
  }
}

# Add permutation null data column
data$Corr_Nulls = c(GP_null_A,GP_null_B,F1_null_A,F1_null_B,F2_null_A,F2_null_B)

# Establish plot groups
data$group = as.factor(c(matrix(0.5, 1, 1000), matrix(1.5, 1, 1000), matrix(2.5, 1, 1000), matrix(3.5, 1, 1000), matrix(4.5, 1, 1000), matrix(5.5, 1, 1000)))

# Ensure corr values are numeric
data$Corr_Actual = as.numeric(data$Corr_Actual)
data$Corr_Nulls = as.numeric(data$Corr_Nulls)

# Permutation Null scatter plots 
fig <- ggplot(data=data, aes(y=Corr_Actual, x=group, fill=group)) +
  geom_jitter(aes(y=Corr_Nulls), color='gray', width=.15, size=.25, alpha=0.8) + #null scatter
  geom_boxplot(aes(y=Corr_Nulls),color="black", width=.22, alpha = 0) + #null boxplot
  geom_point(aes(y=Corr_Actual, size=6, color=group, fill=group)) + #actual corr value
  scale_color_manual(values = c("#FFC000","#FFC000", "#CF7172","#CF7172", "#867DB1","#867DB1")) +
  scale_fill_manual(values = c("#FFC000","#FFC000", "#CF7172","#CF7172", "#867DB1","#867DB1")) +
  labs(x = "", y = expression(paste("Correlation (", italic("r"), ")"))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 27, color = 'black'),
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.title=element_text(size = 27)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  scale_y_continuous(breaks = c(0, 0.1)) + 
  scale_x_discrete(labels = c("P-factor Discovery","P-factor Replication","PRS-F1 Discovery","PRS-F1 Replication","PRS-F2 Discovery","PRS-F2 Replication")) 
fig + coord_cartesian(ylim = c(-0.08, 0.15)) #set consistent axes 
ggsave('Corr_perm_test_A_and_B_scatter_gray.png', width = 20, height = 16, dpi = 1000, units = "cm")

# Permutation Null violin plots 
fig <- ggplot(data=data, aes(y=Corr_Actual, x=group, fill=group)) +
  geom_flat_violin(aes(y=Corr_Nulls,x=group), fill='gray', alpha=.65, width=1) + #null dist violin
  geom_linerange(aes(ymin=-0.20, ymax=0.20), color='black', linewidth=0.3) + #"x-axis" lines for dist violins
  geom_point(aes(y=Corr_Actual, size=5, color=group, fill=group)) + #actual corr value
  scale_color_manual(values = c("#FFC000","#FFC000", "#CF7172","#CF7172", "#867DB1","#867DB1")) +
  scale_fill_manual(values = c("#FFC000","#FFC000", "#CF7172","#CF7172", "#867DB1","#867DB1")) +
  labs(x = "", y = expression(paste("Correlation (", italic("r"), ")"))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 27, color = 'black'),
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.title=element_text(size = 27)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) + 
  scale_x_discrete(labels = c("P-factor Discovery","P-factor Replication","PRS-F1 Discovery","PRS-F1 Replication","PRS-F2 Discovery","PRS-F2 Replication")) 
fig + coord_cartesian(ylim = c(-0.08, 0.15)) #set consistent axes 
ggsave('Corr_perm_test_A_and_B_violin.png', width = 20, height = 16, dpi = 1000, units = "cm")



# Actual corrs for pooled data between A and B
PB1_mat_pool <- readMat('../../Step_2/2_Server_All/BASELINE/results_pooled_050224/General_PB1_acc_pooled.mat')
F1_mat_pool <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_pooled_050224/PRS_1_acc_pooled.mat')
F2_mat_pool <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_pooled_050224/PRS_2_acc_pooled.mat')

# Load in actual pooled corr values and load them 1000x each
GP_pool <- PB1_mat_pool$acc.pooled
F1_pool <- F1_mat_pool$acc.pooled
F2_pool <- F2_mat_pool$acc.pooled
# Create dataframe with column for actual pred acc pooled
data_pooled = data.frame(Corr_Actual = c(as.numeric(matrix(GP_pool, 1, 100)),as.numeric(matrix(F1_pool, 1, 100)),as.numeric(matrix(F2_pool, 1, 100))))

# 2F-CV boot data (pooled)
PB1_mat_boot <- readMat('../../Step_2/2_Server_All/BASELINE/results_boot_040524/results_matchedsamples_boot_040524/all_network/General_PB1_acc_boot.mat')
F1_mat_boot <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_boot_042424/all_network/PRS_1_acc_boot.mat')
F2_mat_boot <- readMat('../../../../Alexander-Bloch/F1+F2_Scripts+files/Step_2-Ridge_reg/2_Server_All/FULL/results_boot_042424/all_network/PRS_2_acc_boot.mat')

# load boot corrs
GP_boot <- t(PB1_mat_boot$PB.test)
F1_boot <- t(F1_mat_boot$PB.test)
F2_boot <- t(F2_mat_boot$PB.test)
data_pooled$Corr_Boot = c(GP_boot,F1_boot,F2_boot) #load into dataframe

Mean_boot <- c(mean(GP_boot), mean(F1_boot), mean(F2_boot))
write.csv(Mean_boot, file = paste0("Mean_boot_corr.csv")) #save out means
Median_boot <- c(median(GP_boot), median(F1_boot), median(F2_boot))
write.csv(Median_boot, file = paste0("Median_boot_corr.csv")) #save out means

# Establish plot groups
data_pooled$group = as.factor(c(matrix(0.5, 1, 100), matrix(1.5, 1, 100), matrix(2.5, 1, 100)))
        
# 2F-CV boot violin plots
fig <- ggplot(data=data_pooled, aes(y=Corr_Boot, x=group, fill=group)) + #null dist
#  geom_flat_violin(#position=position_nudge(x=.2, y=0), 
#                   fill='gray', alpha=.8, width=1) + #boot null dist violin
  geom_flat_violin(aes(y=Corr_Boot,x=group,fill=group), #position=position_nudge(x=.2, y=0), 
                   alpha=.65, width=1) + #bootstrap dist violin
  geom_linerange(aes(ymin=-0.20, ymax=0.20), color='black', linewidth=0.3) + #"x-axis" lines for dist violins
  geom_hline(yintercept=0, linetype="dashed", color='gray', size=0.5) +
#  geom_jitter(aes(y=Corr_Null, color=group, fill=group), width=.05, size=.5, alpha=0.8) + #null scatter
#  geom_boxplot(color="black", width=.15, alpha = 0) + #null boxplot
#  geom_point(aes(y=Corr_Actual, size=6, color=group, fill=group)) + #actual corr value
  scale_color_manual(values = c("#FFC000", "#CF7172", "#867DB1")) +
  scale_fill_manual(values = c("#FFC000", "#CF7172", "#867DB1")) +
  labs(x = "", y = expression(paste("Correlation (", italic("r"), ")"))) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 27, color = 'black'),
	axis.text.x = element_text(size = 27, color = 'black'),
        axis.title=element_text(size = 27)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  scale_y_continuous(breaks = c(-0.1, 0, 0.1)) + 
  scale_x_discrete(labels = c("P-factor", "PRS-F1", "PRS-F2")) 
fig + coord_cartesian(ylim = c(-0.08, 0.15)) #set consistent axes 
ggsave('Corr_boot_violin.png', width = 12, height = 16, dpi = 600, units = "cm")

