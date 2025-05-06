

#Load in Packages
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
mydata1 <- read.csv(file="R Programming results copy.csv",header=T)

#Transpose dataset
microdata1 <- pivot_longer(mydata3, 2:226, names_to = "Taxonomy", values_to = "Frequency")

#Create Lysis method 'factor' column
microdata1$Lysis <- factor(ifelse(grepl("6|7|8|9|10", microdata1$Sample_ID), "Gram-positive Lysis", "Gram-negative Lysis"))

#Create Treatment Column
microdata1$Treatment <- factor(case_when(grepl("^A", microdata1$Sample_ID) ~ "No Spike Control",
                                        grepl("^B", microdata1$Sample_ID) ~ "Enumeration Spike Control",
                                        grepl("^C", microdata1$Sample_ID) ~ "0.1cfu",
                                        grepl("^D", microdata1$Sample_ID) ~ "1cfu",
                                        grepl("^E", microdata1$Sample_ID) ~ "10cfu",
                                        grepl("^F", microdata1$Sample_ID) ~ "100cfu",
                                        grepl("^G", microdata1$Sample_ID) ~ "1000cfu",
                                        grepl("^L", microdata1$Sample_ID) ~ "Listeria Presence",
                                        grepl("^S", microdata1$Sample_ID) ~ "Salmonella Presence",
                                        grepl("^Z", microdata1$Sample_ID) ~ "Campylobacter Presence",
                                        TRUE ~ "Other"))

#Download data set as csv
write.csv(microdata1, file = "microdata1.csv", row.names = FALSE)

#Calculate the sum of reads, richness, evenness, shannon diversity
communitydata1 <- microdata1 %>%
  group_by(Sample_ID) %>%
  summarise(Reads = sum(Frequency), Richness = specnumber(Frequency), Shannon = diversity(Frequency, index = "shannon"), Evenness = diversity(Frequency) / log(specnumber(Frequency)))

#Join results to database
communitydata1 <- communitydata1 %>%
  left_join(microdata1 %>% select(Sample_ID, Lysis, Treatment), by = "Sample_ID") %>%
  distinct()

#reorder 
communitydata1 <- communitydata1 %>%
  select(Sample_ID, Treatment, Lysis, Reads, Richness, Shannon, Evenness)

#Remove Sample_ID column
communitydata1 <- communitydata1[, -1]

#Download data set as csv
write.csv(communitydata1, file = "communitydata1.csv", row.names = FALSE)

# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(communitydata1, Richness, Reads, Shannon, Evenness)  

#inspect the p value in the output
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 4 × 3
#variable statistic             p
#<chr>        <dbl>         <dbl>
#1 Evenness     0.823 0.00000000687 - no
#2 Reads        0.897 0.00000380  - no 
#3 Richness     0.984 0.366 - yes     
#4 Shannon      0.840 0.0000000250 - no


#Create Histogram to show distribution
hist(communitydata1$Evenness, ylab="Frequency", xlab="Evenness", main = "Histogram of Evenness") 
hist(communitydata1$Reads, ylab="Frequency", xlab="Reads", main = "Histogram of Reads")
hist(communitydata1$Richness, ylab="Frequency", xlab="Richness", main = "Histogram of Richness")
hist(communitydata1$Shannon, ylab="Frequency", xlab="Shannon Diversity", main = "Histogram of Shannon Diversity")



# test for equal variance
levene_test(communitydata1, Evenness ~ Treatment)
# p = 4.49e-12 - no
levene_test(communitydata1, Evenness ~ Lysis)
# p = 0.538 - yes
levene_test(communitydata1, Reads ~ Treatment)
# p = 0.175 - yes
levene_test(communitydata1, Reads ~ Lysis)
# p = 0.594 - yes
levene_test(communitydata1, Richness ~ Treatment)
# p = 0.0171 - no
levene_test(communitydata1, Richness ~ Lysis)
# p = 0.437 - yes
levene_test(communitydata1, Shannon ~ Treatment)
# p = 0.00000361 - no
levene_test(communitydata1, Shannon ~ Lysis)
# p = 0.740 - yes

# if data are not both normal and equal variance - transform the data.
shapiro_test(rank(communitydata1$Evenness))
# p =  0.00383 - no
shapiro_test(rank(communitydata1$Reads))
# p = 0.00383 - no
shapiro_test(rank(communitydata1$Shannon))
# p = 0.00383 - no
levene_test(communitydata1, rank(Evenness) ~ Treatment)
# p = 0.000559 - no
levene_test(communitydata1, rank(Richness) ~ Treatment)
# p = 0.00966 - no
levene_test(communitydata1, rank(Shannon) ~ Treatment)
# p = 0.0186 - no

#Conduct the Kruskal Wallis test
kruskal_test(data = communitydata1, Evenness ~ Treatment)
# p = 0.00000000194 - yes
kruskal_test(data = communitydata1, Evenness ~ Lysis)
# p = 0.105 - no
kruskal_test(data = communitydata1, Reads ~ Treatment)
# 0.000506 - yes
kruskal_test(data = communitydata1, Reads ~ Lysis)
# p = 0.9 - no
kruskal_test(data = communitydata1, Richness ~ Treatment)
# p = 0.000000585 - yes
kruskal_test(data = communitydata1, Richness ~ Lysis)
# p = 0.000016 - yes
kruskal_test(data = communitydata1, Shannon ~ Treatment)
# p = 0.0000000106 - yes
kruskal_test(data = communitydata1, Shannon ~ Lysis)
# p = 0.00796 - yes

#If the p value reports there is a significant difference - determine 'size' of effect
kruskal_effsize(data = communitydata1, Evenness ~ Treatment)
# effsize = 64% - large
kruskal_effsize(data = communitydata1, Reads ~ Treatment)
# effsize = 27% - large
kruskal_effsize(data = communitydata1, Richness ~ Treatment)
# effsize = 48% - large
kruskal_effsize(data = communitydata1, Richness ~ Lysis)
# effsize = 21% - large
kruskal_effsize(data = communitydata1, Shannon ~ Treatment)
# effsize = 59% - large
kruskal_effsize(data = communitydata1, Shannon ~ Lysis)
# effsize = 7% - moderate



# If the K-W tests report a significant difference, and there are more then 2 groups, you can test which groups differ from which with a post-hoc Dunn test
EvennessDunn1 <- dunn_test(data = communitydata1, Evenness ~ Treatment)
ReadsDunn1 <- dunn_test(data = communitydata1, Reads ~ Treatment)
RichnessDunn1 <- dunn_test(data = communitydata1, Richness ~ Treatment)
ShannonDunn1 <- dunn_test(data = communitydata1, Shannon ~ Treatment)

#Download data set as csv
write.csv(EvennessDunn1, file = "EvennessDunn1.csv", row.names = FALSE)
write.csv(ReadsDunn1, file = "ReadsDunn1.csv", row.names = FALSE)
write.csv(RichnessDunn1, file = "RichnessDunn1.csv", row.names = FALSE)
write.csv(ShannonDunn1, file = "ShannonDunn1.csv", row.names = FALSE)


#Plot Lysis data

Lysisevenplot1 <- ggplot(communitydata1, aes(x = Lysis, y = Evenness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
 labs(x = "Lysis Method", y = "Evenness") + scale_color_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisevenplot1

Lysisreadsplot1 <- ggplot(communitydata1, aes(x = Lysis, y = Reads)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Sum of Reads") + scale_color_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
  Lysisreadsplot1

Lysisrichplot1 <- ggplot(communitydata1, aes(x = Lysis, y = Richness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Richness") + scale_color_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisrichplot1

Lysisshanplot1 <- ggplot(communitydata1, aes(x = Lysis, y = Shannon)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Treatment), width = 0.2, size = 1.5) +
  labs(x = "Lysis Method",y = "Shannon Diversity") + scale_color_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) + guides(color = guide_legend(title = NULL)) + scale_x_discrete(labels = c("Gram-negative Lysis" = "-", "Gram-positive Lysis" = "+"))
Lysisshanplot1

#save plots
ggsave("Lysisreadsplot1.png", plot = Lysisreadsplot1, width = 8, height = 6, dpi = 300)
ggsave("Lysisrichplot1.png", plot = Lysisrichplot1, width = 8, height = 6, dpi = 300)
ggsave("Lysisshanplot1.png", plot = Lysisshanplot1, width = 8, height = 6, dpi = 300)


#Plot Treatment data

treatevenplot1 <- ggplot(communitydata1, aes(x = Treatment, y = Evenness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Lysis), width = 0.2, size = 1.5) + guides(color = guide_legend(title = "Lysis Method")) + scale_color_discrete(labels= c("-", "+")) +
  labs(x = "Treatment", y = "Evenness") + scale_x_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" ))
treatevenplot1

treatreadsplot1 <- ggplot(communitydata1, aes(x = Treatment, y = Reads)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Lysis), width = 0.2, size = 1.5) +
  labs(x = "Treatment",y = "Sum of Reads") + guides(color = guide_legend(title = "Lysis Method")) + scale_color_discrete(labels = c("-", "+")) + scale_x_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" ))
treatreadsplot1

treatrichplot1 <- ggplot(communitydata1, aes(x = Treatment, y = Richness)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Lysis), width = 0.2, size = 1.5) + scale_x_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) +
  labs(x = "Treatment",y = "Richness")  + guides(color = guide_legend(title = "Lysis Method")) + scale_color_discrete(labels= c("-", "+"))
treatrichplot1

treatshanplot1 <- ggplot(communitydata1, aes(x = Treatment, y = Shannon)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color = Lysis), width = 0.2, size = 1.5) + scale_x_discrete(limits=c("No Spike Control", "Enumeration Spike Control", "0.1cfu", "1cfu", "10cfu", "100cfu", "1000cfu", "Listeria Presence", "Salmonella Presence", "Campylobacter Presence" )) +
  labs(x = "Treatment",y = "Shannon Diversity") + guides(color = guide_legend(title = "Lysis Method")) + scale_color_discrete(labels= c("-", "+"))
treatshanplot1

#save plots
ggsave("treatreadsplot1.png", plot = treatreadsplot1, width = 8, height = 6, dpi = 300)
ggsave("treatrichplot1.png", plot = treatrichplot1, width = 8, height = 6, dpi = 300)
ggsave("treatshanplot1.png", plot = treatshanplot1, width = 8, height = 6, dpi = 300)
