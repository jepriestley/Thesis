

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
microdata2 <- read.csv(file="R Programming results.csv",header=T)

#Transpose dataset
microdata2 <- pivot_longer(microdata2, 2:226, names_to = "Taxonomy", values_to = "Frequency")

#Create Lysis method 'factor' column
microdata2$Lysis <- factor(ifelse(grepl("6|7|8|9|10", microdata2$Sample_ID), "Gram-positive Lysis", "Gram-negative Lysis"))

#Create Treatment Column
microdata2$Treatment <- factor(case_when(grepl("^A", microdata2$Sample_ID) ~ "No Spike Control",
                                        grepl("^B", microdata2$Sample_ID) ~ "Enumeration Spike Control",
                                        grepl("^C", microdata2$Sample_ID) ~ "0.1cfu",
                                        grepl("^D", microdata2$Sample_ID) ~ "1cfu",
                                        grepl("^E", microdata2$Sample_ID) ~ "10cfu",
                                        grepl("^F", microdata2$Sample_ID) ~ "100cfu",
                                        grepl("^G", microdata2$Sample_ID) ~ "1000cfu",
                                        TRUE ~ "Other"))


#Download data set as csv
write.csv(microdata2, file = "microdata2.csv", row.names = FALSE)

#Calculate the sum of reads, richness, evenness, shannon diversity
communitydata2 <- microdata2 %>%
  group_by(Sample_ID) %>%
  summarise(Reads = sum(Frequency), Richness = specnumber(Frequency), Shannon = diversity(Frequency, index = "shannon"), Evenness = diversity(Frequency) / log(specnumber(Frequency)))

#Join results to database
communitydata2 <- communitydata2 %>%
  left_join(microdata2 %>% select(Sample_ID, Lysis, Treatment), by = "Sample_ID") %>%
  distinct()

#reorder 
communitydata2 <- communitydata2 %>%
  select(Sample_ID, Treatment, Lysis, Reads, Richness, Shannon, Evenness)

#Remove Sample_ID column
communitydata2 <- communitydata2[, -1]

#Download data set as csv
write.csv(communitydata2, file = "communitydata2.csv", row.names = FALSE)

#mean and sd
summary_stats2 <- communitydata2 %>%
  group_by(Treatment) %>%
  summarise(
    Mean = mean(Richness),
    SD = sd(Richness)
  )

print(summary_stats2)

# A tibble: 7 × 3
#Treatment                  Mean    SD
#<fct>                     <dbl> <dbl>
#1 0.1cfu                     32.5  4.81
#2 1000cfu                    33.6  2.50
#3 100cfu                     34    2.87
#4 10cfu                      34.7  6.43
#5 1cfu                       31    3.06
#6 Enumeration Spike Control  31.4  5.23
#7 No Spike Control           26.1  5.00

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(communitydata2, Richness, Evenness, Shannon, Reads)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 


# A tibble: 4 × 3
#variable statistic         p
#<chr>        <dbl>     <dbl>
#1 Evenness     0.958 0.0190 - no  
#2 Reads        0.902 0.0000462 - no
#3 Richness     0.983 0.453 - yes 
#4 Shannon      0.978 0.265 - yes



#Create Histogram to show distribution
hist(communitydata2$Evenness, ylab="Frequency", xlab="Evenness", main = "Histogram of Evenness") 
hist(communitydata2$Reads, ylab="Frequency", xlab="Reads", main = "Histogram of Reads")
hist(communitydata2$Richness, ylab="Frequency", xlab="Richness", main = "Histogram of Richness")
hist(communitydata2$Shannon, ylab="Frequency", xlab="Shannon Diversity", main = "Histogram of Shannon Diversity")



# test for equal variance
levene_test(communitydata2, Evenness ~ Treatment)
# p = 0.00562 - no
levene_test(communitydata2, Evenness ~ Lysis)
# p = 0.393 - yes
levene_test(communitydata2, Reads ~ Treatment)
# p = 0.185 - yes
levene_test(communitydata2, Reads ~ Lysis)
# p = 0.347 - yes
levene_test(communitydata2, Richness ~ Treatment)
# p = 0.0317 - no
levene_test(communitydata2, Richness ~ Lysis)
# p = 0.535 - yes
levene_test(communitydata2, Shannon ~ Treatment)
# p = 0.132 - yes
levene_test(communitydata2, Shannon ~ Lysis)
# p = 0.0000820 - no

# if data are not both normal and equal variance - transform the data.
shapiro_test(rank(communitydata2$Evenness))
# p =  0.0134 - no
shapiro_test(rank(communitydata2$Reads))
# p = 0.0134 - no
levene_test(communitydata2, rank(Evenness) ~ Treatment)
# p = 0.00954 - no
levene_test(communitydata2, rank(Richness) ~ Treatment)
# p = 0.230 - yes
levene_test(communitydata2, rank(Shannon) ~ Lysis)
# p = 0.685 - yes

lmrich2 <- lm(rank(Richness) ~ Treatment * Lysis, data = communitydata2)
lmshan2 <- lm(rank(Shannon) ~ Treatment * Lysis, data = communitydata2)

anova(lmrich2)
#Analysis of Variance Table

#Response: rank(Richness)
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Treatment        6  6956.8  1159.5  8.6405 1.163e-06 ***
#Lysis            1 11934.2 11934.2  3.663e-13 ***
#Treatment:Lysis  6  2020.2   336.7  2.5091   0.03201 *  
#Residuals       56  7514.7   134.2                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Compare individual treatment richness based on lysis method:
nospike <- "No Spike Control"
nospike_data1 <- subset(communitydata2, Treatment == nospike)

spike <- "Enumeration Spike Control"
spike_data1 <- subset(communitydata2, Treatment == spike)

enu0.1 <- "0.1cfu"
enu0.1_data1 <- subset(communitydata2, Treatment == enu0.1)

enu1 <- "1cfu"
enu1_data1 <- subset(communitydata2, Treatment == enu1)

enu10 <- "10cfu"
enu10_data1 <- subset(communitydata2, Treatment == enu10)

enu100 <- "100cfu"
enu100_data1 <- subset(communitydata2, Treatment == enu100)

enu1000 <- "1000cfu"
enu1000_data1 <- subset(communitydata2, Treatment == enu1000)


# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(nospike_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.880 0.130 - yes

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(spike_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.923 0.379 - yes

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(enu0.1_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.944 0.597 - yes

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(enu1_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.979 0.960 - yes


# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(enu10_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.946 0.621 - yes

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(enu100_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.915 0.314 - yes

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(enu1000_data1, Richness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 1 × 3
#variable statistic     p
#<chr>        <dbl> <dbl>
#1 Richness     0.924 0.384 - yes


# test for equal variance
levene_test(nospike_data1, Richness ~ Lysis)
# p = 0.242 - yes
levene_test(spike_data1, Richness ~ Lysis)
# p = 0.103 - yes
levene_test(enu0.1_data1, Richness ~ Lysis)
# p = 0.697 - yes
levene_test(enu1_data1, Richness ~ Lysis)
# p = 1  - yes
levene_test(enu10_data1, Richness ~ Lysis)
# p = 0.524 - yes
levene_test(enu100_data1, Richness ~ Lysis)
# p = 0.867 - yes
levene_test(enu1000_data1, Richness ~ Lysis)
# p = 0.481 - yes

#one-way anova
lmrichnospike <- lm(Richness ~ Lysis, data = nospike_data1)

anova(lmrichnospike)

#Analysis of Variance Table

#Response: Richness
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Lysis      1  202.5   202.5  72.321 2.805e-05 ***
#Residuals  8   22.4     2.8                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#one-way anova
lmrichspike <- lm(Richness ~ Lysis, data = spike_data1)

anova(lmrichspike)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Lysis      1   62.5  62.500  25.641 0.0009726 ***
#Residuals  8   19.5   2.438                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#one-way anova
lmrichenum0.1 <- lm(Richness ~ Lysis, data = enu0.1_data1)

anova(lmrichenum0.1)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Lysis      1   62.5  62.500  27.027 0.0008237 ***
#Residuals  8   18.5   2.312                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#one-way anova
lmrichenum1 <- lm(Richness ~ Lysis, data = enu1_data1)

anova(lmrichenum1)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value Pr(>F)
#Lysis      1   12.1  12.100  1.3948 0.2715
#Residuals  8   69.4   8.675    

#one-way anova
lmrichenum10 <- lm(Richness ~ Lysis, data = enu10_data1)

anova(lmrichenum10)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value   Pr(>F)   
#Lysis      1   57.6  57.600  19.692 0.002175 **
#Residuals  8   23.4   2.925                    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#one-way anova
lmrichenum100 <- lm(Richness ~ Lysis, data = enu100_data1)

anova(lmrichenum100)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value    Pr(>F)    
#Lysis      1   62.5  62.500  34.483 0.0003733 ***
#Residuals  8   14.5   1.813                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#one-way anova
lmrichenum1000 <- lm(Richness ~ Lysis, data = enu1000_data1)

anova(lmrichenum1000)

#Analysis of Variance Table

#Response: rank(Richness)
#Df Sum Sq Mean Sq F value Pr(>F)
#Lysis      1   19.6 19.6000  2.5747 0.1473
#Residuals  8   60.9  7.6125 

anova(lmshan2)

#Analysis of Variance Table

#Response: rank(Shannon)
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Treatment        6  9833.6  1638.9  7.6179 5.294e-06 ***
#Lysis            1  4657.7  4657.7 21.6495 2.046e-05 ***
#Treatment:Lysis  6  2038.2   339.7  1.5789    0.1704    
#Residuals       56 12048.0   215.1                      
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#If the p value reports there is a significant difference - determine 'size' of effect
eta_squared(lmrich2)
# effsize = 24% (treatment), 42% (lysis), 7% (treatment:lysis)
eta_squared(lmshan2)
# effsize = 34% (treatment), 16% (lysis), 7% (treatment:lysis)


# If the Anova tests report a significant difference, and there are more then 2 groups, you can test which groups differ from which.
richtukey2 <- tukey_hsd(communitydata2, rank(Richness) ~ Treatment)
shantukey2 <- tukey_hsd(communitydata2, rank(Shannon) ~ Treatment)

#Download dataset as csv
write.csv(richtukey2, file = "Richnesstukey2.csv", row.names = FALSE)
write.csv(shantukey2, file = "Shannontukey2.csv", row.names = FALSE)


#Conduct the Kruskal Wallis test
kruskal_test(data = communitydata2, Evenness ~ Treatment)
# p = 0.00000602 - yes
kruskal_test(data = communitydata2, Evenness ~ Lysis)
# p = 0.106 - no
kruskal_test(data = communitydata2, Reads ~ Treatment)
# p = 0.00221 - yes
kruskal_test(data = communitydata2, Reads ~ Lysis)
# p = 0.846 - no
kruskal_test(data = communitydata2, Richness ~ Treatment)
# p = 0.00971 - yes
kruskal_test(data = communitydata2, Richness ~ Lysis)
# p = 0.0000000736 - yes
kruskal_test(data = communitydata2, Shannon ~ Treatment)
# p = 0.000582  - yes
kruskal_test(data = communitydata2, Shannon ~ Lysis)
# p = 0.000798 - yes


#If the p value reports there is a significant difference - determine 'size' of effect
kruskal_effsize(data = communitydata2, Evenness ~ Treatment)
# effsize = 45% - large
kruskal_effsize(data = communitydata2, Reads ~ Treatment)
# effsize = 21% - large
kruskal_effsize(data = communitydata2, Richness ~ Treatment)
# effsize = 17% - large
kruskal_effsize(data = communitydata2, Richness ~ Lysis)
# effsize = 41% - large
kruskal_effsize(data = communitydata2, Shannon ~ Treatment)
# effsize = 28% - large
kruskal_effsize(data = communitydata2, Shannon ~ Lysis)
# effsize = 15% - large


# If the K-W tests report a significant difference, and there are more then 2 groups, you can test which groups differ from which with a post-hoc Dunn test
EvennessDunn2 <- dunn_test(data = communitydata2, Evenness ~ Treatment)
ReadsDunn2 <- dunn_test(data = communitydata2, Reads ~ Treatment)
RichnessDunn2 <- dunn_test(data = communitydata2, Richness ~ Treatment)
ShannonDunn2 <- dunn_test(data = communitydata2, Shannon ~ Treatment)


#Download data set as csv
write.csv(EvennessDunn2, file = "EvennessDunn2.csv", row.names = FALSE)
write.csv(ReadsDunn2, file = "ReadsDunn2.csv", row.names = FALSE)
write.csv(RichnessDunn2, file = "RichnessDunn2.csv", row.names = FALSE)
write.csv(ShannonDunn2, file = "ShannonDunn2.csv", row.names = FALSE)

#plot data

theme_set(theme_classic())

Lysisrichplot2 <- ggplot(communitydata2, aes(x = Lysis, y = Richness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + theme(legend.position = "none")  +
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Richness") + scale_color_npg(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisrichplot2

Lysisshanplot2 <- ggplot(communitydata2, aes(x = Lysis, y = Shannon)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  theme(legend.position = "none")  +
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Shannon Diversity") + scale_color_npg(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisshanplot2

# Extract the legend
legend2 <- get_legend(Lysisrichplot2)

# Save the legend as a separate file
ggsave("legend_only_2.png", legend2, width = 3, height = 4, dpi = 300)


final_plot_2 <- plot_grid(
  Lysisrichplot2, 
  Lysisshanplot2,
  legend2,
  nrow = 1, 
  labels = c("A", "B", ""),
  rel_heights = c(1), rel_widths = c(2,2,1)
)

final_plot_2
ggsave("final_plot_2.png", plot = final_plot_2, width = 10, height = 6, dpi = 300)
