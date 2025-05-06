install.packages(c("tidyr", "tidyverse", "dplyr", "ggplot2", "vegan", "stringr",
                   "ggthemes", "ggrepel", "scales", "ggpattern", "lemon",
                   "RColorBrewer", "rstatix", "ggsci"))




#Load in Packages
library(grid)
library(cowplot)
library(gridExtra)
library(patchwork)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)
library(ggthemes)
library(ggrepel)
library(scales)
library(ggpattern)
library(lemon)
library(RColorBrewer)
library(rstatix)
library(ggsci)

#Load in dataset
mydata3 <- read.csv(file="Presence_absence.csv",header=T)

#Transpose dataset
microdata3 <- pivot_longer(mydata3, 2:226, names_to = "Taxonomy", values_to = "Frequency")

#Create Lysis method 'factor' column
microdata3$Lysis <- factor(ifelse(grepl("6|7|8|9|10", microdata3$Sample_ID), "Gram-positive Lysis", "Gram-negative Lysis"))

#Create Treatment Column
microdata3$Treatment <- factor(case_when(grepl("^L", microdata3$Sample_ID) ~ "Listeria Presence",
                                        grepl("^S", microdata3$Sample_ID) ~ "Salmonella Presence",
                                        grepl("^Z", microdata3$Sample_ID) ~ "Campylobacter Presence",
                                        TRUE ~ "Other"))
  
#Calculate the sum of reads, richness, evenness, shannon diversity
communitydata3 <- microdata3 %>%
  group_by(Sample_ID) %>%
  summarise(Reads = sum(Frequency), Richness = specnumber(Frequency), Shannon = diversity(Frequency, index = "shannon"), Evenness = diversity(Frequency) / log(specnumber(Frequency)))

#Join results to database
communitydata3 <- communitydata3 %>%
  left_join(microdata3 %>% select(Sample_ID, Lysis, Treatment), by = "Sample_ID") %>%
  distinct()

#reorder 
communitydata3 <- communitydata3 %>%
  select(Sample_ID, Treatment, Lysis, Reads, Richness, Shannon, Evenness)

#Remove Sample_ID column
communitydata3 <- communitydata3[, -1]

#Download data set as csv
write.csv(communitydata3, file = "communitydata3.csv", row.names = FALSE)

#mean and sd
summary_stats3 <- communitydata3 %>%
  group_by(Treatment) %>%
  summarise(
    Mean = mean(Shannon),
    SD = sd(Shannon)
  )

print(summary_stats3)

# A tibble: 3 × 3
#Treatment               Mean    SD
#<fct>                  <dbl> <dbl>
#1 Campylobacter Presence  21    4.60
#2 Listeria Presence       21.7  1.51
#3 Salmonella Presence     26.8  2.40

#mean and sd
summary_stats4 <- communitydata3 %>%
  group_by(Lysis) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness)
  )

print(summary_stats4)

# A tibble: 2 × 3
#Lysis                Mean    SD
#<fct>               <dbl> <dbl>
#1 Gram-negative Lysis  22.4  4.80
#2 Gram-positive Lysis  23.9  3.06

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(communitydata3, Richness, Evenness, Shannon, Reads)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 


# A tibble: 4 × 3
#variable statistic       p
#<chr>        <dbl>   <dbl>
#1 Evenness     0.820 0.00293 - no
#2 Reads        0.877 0.0236 - no
#3 Richness     0.931 0.204  - yes
#4 Shannon      0.814 0.00238 - no



#Create Histogram to show distribution
hist(communitydata3$Evenness, ylab="Frequency", xlab="Evenness", main = "Histogram of Evenness") 
hist(communitydata3$Reads, ylab="Frequency", xlab="Reads", main = "Histogram of Reads")
hist(communitydata3$Richness, ylab="Frequency", xlab="Richness", main = "Histogram of Richness")
hist(communitydata3$Shannon, ylab="Frequency", xlab="Shannon Diversity", main = "Histogram of Shannon Diversity")



# test for equal variance
levene_test(communitydata3, Evenness ~ Treatment)
# p = 0.0000437 - no
levene_test(communitydata3, Evenness ~ Lysis)
# p = 0.932 - yes
levene_test(communitydata3, Reads ~ Treatment)
# p = 0.383 - yes
levene_test(communitydata3, Reads ~ Lysis)
# p = 0.849 - yes
levene_test(communitydata3, Richness ~ Treatment)
# p = 0.385 - yes
levene_test(communitydata3, Richness ~ Lysis)
# p = 0.427 - yes
levene_test(communitydata3, Shannon ~ Treatment)
# p = 0.00174 - no
levene_test(communitydata3, Shannon ~ Lysis)
# p = 0.929 - yes

# if data are not both normal and equal variance - transform the data.
shapiro_test(rank(communitydata3$Evenness))
# p =  0.629 - yes
shapiro_test(rank(communitydata3$Reads))
# p =  0.629 - yes
shapiro_test(rank(communitydata3$Shannon))
# p = 0.629 - yes
levene_test(communitydata3, rank(Evenness) ~ Treatment)
# p = 1 - yes
levene_test(communitydata3, rank(Shannon) ~ Treatment)
# p = 1 - yes



lmeven1 <- lm(rank(Evenness) ~ Treatment * Lysis, data = communitydata3)
lmreads1 <- lm(rank(Reads) ~ Treatment * Lysis, data = communitydata3)
lmrich1 <- lm(Richness ~ Treatment * Lysis, data = communitydata3)
lmshan1 <- lm(rank(Shannon) ~ Treatment * Lysis, data = communitydata3)



anova(lmeven1)

#Analysis of Variance Table

#Response: rank(Evenness)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Treatment        2 432.00 216.000   77.76 1.351e-07 ***
#Lysis            1   2.72   2.722    0.98   0.34174    
#Treatment:Lysis  2  16.44   8.222    2.96   0.09017 .  
#Residuals       12  33.33   2.778                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(lmreads1)
#Analysis of Variance Table

#Response: rank(Reads)
#Df Sum Sq Mean Sq F value Pr(>F)
#Treatment        2  22.33  11.167  0.3317 0.7241
#Lysis            1   1.39   1.389  0.0413 0.8425
#Treatment:Lysis  2  56.78  28.389  0.8432 0.4543
#Residuals       12 404.00  33.667 

anova(lmrich1)
#Analysis of Variance Table

#Response: Richness
#Df  Sum Sq Mean Sq F value  Pr(>F)  
#Treatment        2 122.333  61.167  5.5327 0.01983 *
#Lysis            1   9.389   9.389  0.8492 0.37492  
#Treatment:Lysis  2   4.111   2.056  0.1859 0.83268  
#Residuals       12 132.667  11.056                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


anova(lmshan1)

#Analysis of Variance Table

#Response: rank(Shannon)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Treatment        2 432.00 216.000 69.4286 2.533e-07 ***
#Lysis            1   6.72   6.722  2.1607    0.1673    
#Treatment:Lysis  2   8.44   4.222  1.3571    0.2942    
#Residuals       12  37.33   3.111                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#If the p value reports there is a significant difference - determine 'size' of effect
eta_squared(lmeven1)
# effsize = 89% (treatment)
eta_squared(lmrich1)
# effsize = 46% (treatment)
eta_squared(lmshan1)
# effsize = 89% (treatment)


# If the Anova tests report a significant difference, and there are more then 2 groups, you can test which groups differ from which.
Eventukey1 <- tukey_hsd(communitydata3, rank(Evenness) ~ Treatment)
richtukey1 <- tukey_hsd(communitydata3, Richness ~ Treatment)
shantukey1 <- tukey_hsd(communitydata3, rank(Shannon) ~ Treatment)

#Download data set as csv
write.csv(Eventukey1, file = "Evennesstukey1.csv", row.names = FALSE)
write.csv(richtukey1, file = "Richnesstukey1.csv", row.names = FALSE)
write.csv(shantukey1, file = "Shannontukey1.csv", row.names = FALSE)


#Conduct the Kruskal Wallis test
kruskal_test(data = communitydata3, Evenness ~ Treatment)
# p = 0.000511 - yes
kruskal_test(data = communitydata3, Evenness ~ Lysis)
# p = 0.757 - no
kruskal_test(data = communitydata3, Reads ~ Treatment)
# 0.676 - no
kruskal_test(data = communitydata3, Reads ~ Lysis)
# p = 0.825 - no
kruskal_test(data = communitydata3, Richness ~ Treatment)
# p = 0.0134 - yes
kruskal_test(data = communitydata3, Richness ~ Lysis)
# p = 0.503 - no
kruskal_test(data = communitydata3, Shannon ~ Treatment)
# p = 0.000511 - yes
kruskal_test(data = communitydata3, Shannon ~ Lysis)
# p = 0.627 - no


#Test individually:

listeria <- "Listeria Presence"
listeria_data1 <- subset(communitydata3, Treatment == listeria)

campy <- "Campylobacter Presence"
campy_data1 <- subset(communitydata3, Treatment == campy)

salmonella <- "Salmonella Presence"
salmonella_data1 <- subset(communitydata3, Treatment == salmonella)


shapiro_test(listeria_data1, Richness)  
# A tibble: 1 × 3
#variable statistic      p
#<chr>        <dbl>  <dbl>
#1 Richness     0.767 0.0288 -no

shapiro_test(rank(listeria_data1$Richness))
# A tibble: 1 × 3
#variable                      statistic p.value
#<chr>                             <dbl>   <dbl>
#1 rank(listeria_data1$Richness)     0.775  0.0346

shapiro_test(campy_data1, Richness)  
# A tibble: 1 × 3
#variable statistic      p
#<chr>        <dbl>  <dbl>
#1 Richness     0.899 0.368 - yes

shapiro_test(salmonella_data1, Richness)  
# A tibble: 1 × 3
#variable statistic      p
#<chr>        <dbl>  <dbl>
#1 Richness     0.891 0.324 - yes


kruskal_test(data = listeria_data1, Richness ~ Lysis)
# p = 0.239 - no
kruskal_test(data = campy_data1, Richness ~ Lysis)
# p = 0.658 - no
kruskal_test(data = salmonella_data1, Richness ~ Lysis)
# p =  1  - no




#plot the data

Lysisrichplot3 <- ggplot(communitydata3, aes(x = Lysis, y = Richness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + theme(legend.position = "none") +
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Richness") + scale_color_npg(limits=c("Listeria Presence", "Salmonella Presence", "Campylobacter Presence" ), labels = c(expression(italic("Listeria") ~ "Presence"),
                                                                                                                                                        expression(italic("Salmonella") ~ "Presence"),
                                                                                                                                                        expression(italic("Campylobacter") ~ "Presence"))) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisrichplot3

Lysisshanplot3 <- ggplot(communitydata3, aes(x = Lysis, y = Shannon)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + theme(legend.position = "none")  +
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Shannon Diversity") + scale_color_npg(limits=c("Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisshanplot3

LegendkeyListeria <- expression(paste(italic("Listeria"), "Presence"))
LegendkeyCampy <- expression(paste(italic("Campylobacter"), "Presence"))
Legendkeysal <- expression(paste(italic("Salmonella"), "Presence"))

# Extract the legend
legend3 <- get_legend(Lysisrichplot3)

# Save the legend as a separate file
ggsave("legend_only_3.png", legend3, width = 3, height = 4, dpi = 300)

final_plot_3 <- plot_grid(
  Lysisrichplot3, 
  Lysisshanplot3,
  legend3,
  nrow = 1, 
  labels = c("A", "B", ""),
  rel_heights = c(1), rel_widths = c(2,2,1)
)
final_plot_3
ggsave("final_plot_3.png", plot = final_plot_3, width = 10, height = 6, dpi = 300)
