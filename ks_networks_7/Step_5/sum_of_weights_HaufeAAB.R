
library(gifti)
library(R.matlab)
library(RNifti)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/workbench')
library(ggplot2)

dirin <- "../../Step_5/BASELINE/Weight_Haufe_031324/";
dirout <- "sum_of_weights_HaufeAAB/"
dir.create(dirout)

catg <- 'General_PB1'
for (grp in c('A','B')) {
  
  infile <- paste0(dirin, catg, '_', grp, '_Weight_haufetrans_all','.mat')
  
  # read in weights_by_vertex.mat
  readin <- readMat(infile);
  read_weights <- readin$PB.test
  
  #format Haufe weights into a matrix of 17 PFNs
  all_weights <- matrix(0,59412,17)
  for (vert in c(1:59412)) {
    for (PFN in c(1:17)) {
      all_weights[vert,PFN] <- read_weights[vert+(PFN-1)*59412]
    }
  }
  
  if (grp == 'A') {
    all_weights_A <- matrix(0,59412,17)
    all_weights_A <- all_weights
  }
  if (grp == 'B') {
    all_weights_B <- matrix(0,59412,17)
    all_weights_B <- all_weights
  }
}

#average Haufe weights between groups A and B
all_weights_avg <- matrix(0,59412,17) 
for (vert in c(1:59412)) {
  for (PFN in c(1:17)) {
    all_weights_avg[vert,PFN] <- sum(all_weights_A[vert,PFN]+all_weights_B[vert,PFN])/2
  }
}

sum_weights <- matrix(0,17,4)
#matrix of 4 values:
#1. NET- sum preserving sign
#2. ABS- sum of abs values
#3. POS- sum of only positive values
#4. NEG- sum of only negative values
# NET = POS + NEG
# ABS = POS + |NEG|
for (PFN in c(1:17)){
  for (vert in c(1:59412)) {
    sum_weights[PFN,1] <- sum_weights[PFN,1] + all_weights_avg[vert,PFN] #NET
    sum_weights[PFN,2] <- sum_weights[PFN,2] + abs(all_weights_avg[vert,PFN]) #ABS
    if (all_weights_avg[vert,PFN] > 0) {
      sum_weights[PFN,3] <- sum_weights[PFN,3] + all_weights_avg[vert,PFN] #POS
    } else {
      sum_weights[PFN,4] <- sum_weights[PFN,4] + all_weights_avg[vert,PFN] #NEG
    }
  }
}
#write.csv(sum_weights, file = paste0(dirout,"Sum_of_weights_HaufeAAB.csv")) #save out values


PFN_size <- read.csv('../PFN_vert_count.csv') #read in PFN soft parcel vert count
sum_weights_norm <- matrix(0,17,4)
for (PFN in c (1:17)) {
  for (sum in c (1:4)) {
    sum_weights_norm[PFN,sum] <- sum_weights[PFN,sum]/PFN_size[PFN,2] #normalized by PFN size (PFN vert count/avg of PFN vert count)
  }
}
#write.csv(sum_weights_norm, file = paste0(dirout,"Sum_of_weights_vert_count_norm_HaufeAAB.csv")) #save out normed values

#NORMED
# Using 7 colors scheme for bar plot
# bar plot for NET sum of weights
NET_sorted = rank(sum_weights_norm[,1]);
NET_inAscend <- c()
for (t in c(1:17)) { NET_inAscend <- append(NET_inAscend,which(NET_sorted == t))}

data <- data.frame(NET = as.numeric(sum_weights_norm[,1])) 
data$NET_Rank <- as.numeric(NET_sorted)
# BorderColor1 is for colors of 1-17 network numbers(to label the network numbers in a 1-17 order)
BorderColor1 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                  "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                  "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                  "#4E31A8", "#F5BA2E");

# BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
BorderColor2 <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#7499C2", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E443FF", "#FFFFFF",
                  "#FFFFFF", "#FFFFFF", "#7499C2", "#FFFFFF", "#F5BA2E", 
                  "#4E31A8", "#FFFFFF");
# Change the color order
BorderColor <- BorderColor1[NET_inAscend]
BorderColor3 <- BorderColor2[NET_inAscend]

# the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
LineType <- c("dashed", "dashed", "dashed", "solid", "dashed",
              "dashed", "dashed", "dashed", "solid", "dashed",
              "dashed", "dashed", "solid", "dashed", "solid",
              "solid", "dashed");

Fig <- ggplot(data, aes(NET_Rank, NET)) +
  geom_bar(stat = "identity", fill=BorderColor3[NET_sorted], 
           colour = BorderColor[NET_sorted], linetype = LineType, width = 0.8) + 
  labs(x = "Networks", y = paste0("Sum of Weights")) + theme_classic() +
  theme(axis.text.x = ggtext::element_markdown(size= 27, color = BorderColor), 
        axis.text.y = ggtext::element_markdown(size= 33, color = "black"), 
        axis.title=ggtext::element_markdown(size = 33)) +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = as.character(NET_inAscend),limits = as.character(NET_sorted)) 
Fig
ggsave(paste0(dirout,"/NET_pfac_b_vert_norm.tiff"), width = 19, height = 16, dpi = 600, units = "cm")


# bar plot for ABS sum of weights
ABS_sorted = rank(sum_weights_norm[,2]);
ABS_inAscend <- c()
for (t in c(1:17)) { ABS_inAscend <- append(ABS_inAscend,which(ABS_sorted == t))}

data <- data.frame(ABS = as.numeric(sum_weights_norm[,2])) 
#normalize by max
data$ABS <- data$ABS/max(data$ABS)
data$ABS_Rank <- as.numeric(ABS_sorted)
# BorderColor1 is for colors of 1-17 network numbers(to label the network numbers in a 1-17 order)
BorderColor1 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                  "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                  "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                  "#4E31A8", "#F5BA2E");

# BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
BorderColor2 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                  "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                  "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                  "#4E31A8", "#F5BA2E");
# Change the color order
BorderColor <- BorderColor1[ABS_inAscend]
BorderColor3 <- BorderColor2[ABS_inAscend]

# the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
LineType <- c("solid", "solid", "solid", "solid", "solid",
              "solid", "solid", "solid", "solid", "solid",
              "solid", "solid", "solid", "solid", "solid",
              "solid", "solid");

Fig <- ggplot(data, aes(ABS_Rank, ABS)) +
  geom_bar(stat = "identity", fill=BorderColor3[ABS_sorted], 
           colour = BorderColor[ABS_sorted], linetype = LineType, width = 0.8) + 
  labs(x = "Networks", y = paste0("Feature Importance")) + theme_classic() +
  theme(axis.text.x = ggtext::element_markdown(size= 27, color = BorderColor), 
        axis.text.y = ggtext::element_markdown(size= 33, color = "black"), 
        axis.title=ggtext::element_markdown(size = 33)) +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = as.character(ABS_inAscend),limits = as.character(ABS_sorted)) 
Fig
ggsave(paste0(dirout,"/ABS_pfac_b_vert_norm_FEAT_IMP_NORMED_1.tiff"), width = 19, height = 16, dpi = 600, units = "cm")


#SUM OF WEIGHTS IN A AND B SEPARATELY
for (grp in c('A','B')) {
  if (grp == 'A') {
    all_weights <- all_weights_A
  } else if (grp == 'B') {
    all_weights <- all_weights_B
  }
  
  sum_weights <- matrix(0,17,4)
  #matrix of 4 values:
  #1. NET- sum preserving sign
  #2. ABS- sum of abs values
  #3. POS- sum of only positive values
  #4. NEG- sum of only negative values
  # NET = POS + NEG
  # ABS = POS + |NEG|
  for (PFN in c(1:17)){
    for (vert in c(1:59412)) {
      sum_weights[PFN,1] <- sum_weights[PFN,1] + all_weights[vert,PFN] #NET
      sum_weights[PFN,2] <- sum_weights[PFN,2] + abs(all_weights[vert,PFN]) #ABS
      if (all_weights[vert,PFN] > 0) {
        sum_weights[PFN,3] <- sum_weights[PFN,3] + all_weights[vert,PFN] #POS
      } else {
        sum_weights[PFN,4] <- sum_weights[PFN,4] + all_weights[vert,PFN] #NEG
      }
    }
  }
  write.csv(sum_weights, file = paste0(dirout,"Sum_of_weights_",grp,"_HaufeAAB.csv")) #save out values
  
  PFN_size <- read.csv('../PFN_vert_count.csv') #read in PFN soft parcel vert count
  sum_weights_norm <- matrix(0,17,4)
  for (PFN in c (1:17)) {
    for (sum in c (1:4)) {
      sum_weights_norm[PFN,sum] <- sum_weights[PFN,sum]/PFN_size[PFN,2] #normalized by PFN size (PFN vert count/avg of PFN vert count)
    }
  }
  write.csv(sum_weights_norm, file = paste0(dirout,"Sum_of_weights_",grp,"_vert_count_norm_HaufeAAB.csv")) #save out normed values
  
  #NORMED NET
  # Using 7 colors scheme for bar plot
  # bar plot for NET sum of weights
  NET_sorted = rank(sum_weights_norm[,1]);
  NET_inAscend <- c()
  for (t in c(1:17)) { NET_inAscend <- append(NET_inAscend,which(NET_sorted == t))}
  
  data <- data.frame(NET = as.numeric(sum_weights_norm[,1])) 
  data$NET_Rank <- as.numeric(NET_sorted)
  # BorderColor1 is for colors of 1-17 network numbers(to label the network numbers in a 1-17 order)
  BorderColor1 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                    "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                    "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                    "#4E31A8", "#F5BA2E");
  if (grp == 'A') {  #can change based on sig testing of grp A (9, 15, 16)
    # BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
    BorderColor2 <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                      "#FFFFFF", "#FFFFFF", "#FFFFFF", "#E443FF", "#FFFFFF",
                      "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#F5BA2E", 
                      "#4E31A8", "#FFFFFF");
  } else if (grp == 'B') {  #can change based on sig testing of grp B (16)
    # BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
    BorderColor2 <- c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                      "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF",
                      "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF", 
                      "#4E31A8", "#FFFFFF");
  }
  # Change the color order
  BorderColor <- BorderColor1[NET_inAscend]
  BorderColor3 <- BorderColor2[NET_inAscend]
  if (grp == 'A') {  #can change based on sig testing of grp A (9, 15, 16)
    # the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
    LineType <- c("dashed", "dashed", "dashed", "dashed", "dashed",
                  "dashed", "dashed", "dashed", "solid", "dashed",
                  "dashed", "dashed", "dashed", "dashed", "solid",
                  "solid", "dashed");
  } else if (grp == 'B') {  #can change based on sig testing of grp B (16)
    # the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
    LineType <- c("dashed", "dashed", "dashed", "dashed", "dashed",
                  "dashed", "dashed", "dashed", "dashed", "dashed",
                  "dashed", "dashed", "dashed", "dashed", "dashed",
                  "solid", "dashed");
  }
  Fig <- ggplot(data, aes(NET_Rank, NET)) +
    geom_bar(stat = "identity", fill=BorderColor3[NET_sorted], 
             colour = BorderColor[NET_sorted], linetype = LineType, width = 0.8) + 
    labs(x = "Networks", y = paste0("Sum of Weights")) + theme_classic() +
    theme(axis.text.x = ggtext::element_markdown(size= 27, color = BorderColor), 
          axis.text.y = ggtext::element_markdown(size= 33, color = "black"), 
          axis.title=ggtext::element_markdown(size = 33)) +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = as.character(NET_inAscend),limits = as.character(NET_sorted)) 
  Fig
  ggsave(paste0(dirout,"/NET_pfac_sig_",grp,"_vert_norm.tiff"), width = 19, height = 16, dpi = 600, units = "cm")
  
  
  #NORMED ABS
  # Using 7 colors scheme for bar plot
  # bar plot for ABS sum of weights
  ABS_sorted = rank(sum_weights_norm[,2]);
  ABS_inAscend <- c()
  for (t in c(1:17)) { ABS_inAscend <- append(ABS_inAscend,which(ABS_sorted == t))}
  
  data <- data.frame(ABS = as.numeric(sum_weights_norm[,2])) 
  data$ABS_Rank <- as.numeric(ABS_sorted)
  # BorderColor1 is for colors of 1-17 network numbers(to label the network numbers in a 1-17 order)
  BorderColor1 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                    "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                    "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                    "#4E31A8", "#F5BA2E");
  if (grp == 'A') {  #can change based on sig testing of grp A (all)
    # BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
    BorderColor2 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                      "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                      "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                      "#4E31A8", "#F5BA2E");
  } else if (grp == 'B') {  #can change based on sig testing of grp B (all)
    # BorderColor2 is for colors of bar fill 1-17 networks (in a 1-17 order). If the bar is nofill, change the code to "#FFFFFF" for that PFN. ("black" is filled with black)
    BorderColor2 <- c("#E76178", "#7499C2", "#F5BA2E", "#7499C2", "#00A131",
                      "#AF33AD", "#E443FF", "#E76178", "#E443FF", "#AF33AD",
                      "#7499C2", "#E76178", "#7499C2", "#00A131", "#F5BA2E", 
                      "#4E31A8", "#F5BA2E");
  }
  # Change the color order
  BorderColor <- BorderColor1[ABS_inAscend]
  BorderColor3 <- BorderColor2[ABS_inAscend]
  if (grp == 'A') {  #can change based on sig testing of grp A (all)
    # the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
    LineType <- c("solid", "solid", "solid", "solid", "solid",
                  "solid", "solid", "solid", "solid", "solid",
                  "solid", "solid", "solid", "solid", "solid",
                  "solid", "solid");
  } else if (grp == 'B') {  #can change based on sig testing of grp B (all)
    # the line type for the 1-17 network bars (in a 1-17 order).  If a dashed line is needed for a bar, please change the value to "dashed".
    LineType <- c("solid", "solid", "solid", "solid", "solid",
                  "solid", "solid", "solid", "solid", "solid",
                  "solid", "solid", "solid", "solid", "solid",
                  "solid", "solid");
  }
  Fig <- ggplot(data, aes(ABS_Rank, ABS)) +
    geom_bar(stat = "identity", fill=BorderColor3[ABS_sorted], 
             colour = BorderColor[ABS_sorted], linetype = LineType, width = 0.8) + 
    labs(x = "Networks", y = paste0("Sum of Weights")) + theme_classic() +
    theme(axis.text.x = ggtext::element_markdown(size= 27, color = BorderColor), 
          axis.text.y = ggtext::element_markdown(size= 33, color = "black"), 
          axis.title=ggtext::element_markdown(size = 33)) +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = as.character(ABS_inAscend),limits = as.character(ABS_sorted)) 
  Fig
  ggsave(paste0(dirout,"/ABS_pfac_sig_",grp,"_vert_norm.tiff"), width = 19, height = 16, dpi = 600, units = "cm")
}
