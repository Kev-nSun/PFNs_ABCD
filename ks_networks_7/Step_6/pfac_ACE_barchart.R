library(ggplot2)

labels <- c("A","C","E")
ace_values <- c(0.54,0.20,0.26)
colors <- c("#636363","gray","white")
ACE_pfac <- data.frame(name=labels,value=ace_values)

bp <- barplot(height=ACE_pfac$value, 
        names=ACE_pfac$name, 
        col=colors,
        ylab="Proportion of variation",
        ylim=c(0,1),
        space=c(0,0,0)
        )

lower_errors <- c(0.13, 0.11, 0.04)
upper_errors <- c(0.13, 0.10, 0.06)

arrows(bp, ace_values - lower_errors, bp, ace_values + upper_errors, angle = 90, code = 3, length = 0.08, col = "black")