
#Load in dataset
mydata3 <- read.csv(file="R Programming results copy.csv",header=T)

#Transpose dataset
lysisdata3 <- pivot_longer(mydata3, 2:226, names_to = "Taxonomy", values_to = "Frequency")

#Create Lysis method 'factor' column
lysisdata3$Lysis <- factor(ifelse(grepl("6|7|8|9|10", lysisdata3$Sample_ID), "Gram-positive Lysis", "Gram-negative Lysis"))

#Create Treatment Column
lysisdata3$Treatment <- factor(case_when(grepl("^A", lysisdata3$Sample_ID) ~ "No Spike Control",
                                        grepl("^B", lysisdata3$Sample_ID) ~ "Enumeration Spike Control",
                                        grepl("^C", lysisdata3$Sample_ID) ~ "0.1cfu",
                                        grepl("^D", lysisdata3$Sample_ID) ~ "1cfu",
                                        grepl("^E", lysisdata3$Sample_ID) ~ "10cfu",
                                        grepl("^F", lysisdata3$Sample_ID) ~ "100cfu",
                                        grepl("^G", lysisdata3$Sample_ID) ~ "1000cfu",
                                        grepl("^L", lysisdata3$Sample_ID) ~ "Listeria Presence",
                                        grepl("^S", lysisdata3$Sample_ID) ~ "Salmonella Presence",
                                        grepl("^Z", lysisdata3$Sample_ID) ~ "Campylobacter Presence",
                                        TRUE ~ "Other"))


#Calculate All - Alpha Diversity 
Alpha_diversity3 <- lysisdata3 %>%
  group_by(Sample_ID) %>%
  summarise(Reads = sum(Frequency), Richness = specnumber(Frequency), Shannon = diversity(Frequency, index = "shannon"), Evenness = diversity(Frequency) / log(specnumber(Frequency)))


Alpha_diversity3 <- Alpha_diversity3 %>%
  left_join(lysisdata3 %>% select(Sample_ID, Lysis, Treatment), by = "Sample_ID") %>%
  distinct()

Alpha_diversity3 <- Alpha_diversity3 %>%
  select(Sample_ID, Treatment, Lysis, Reads, Richness, Shannon, Evenness)
Alpha_diversity3 <- Alpha_diversity3[, -1]


# test to see if the data are normally distributed with a shaprio-wilk test for each response
shapiro_test(Alpha_diversity3, Richness, Reads, Shannon, Evenness)  

#inspect the p value in the output - this is the info that will need to go into your dissertation
# a p value less then 0.05 means we can reject the hypothesis the data are normally distributed 

# A tibble: 2 × 3
variable statistic          p
<chr>        <dbl>      <dbl>
  1 Reads        0.897 0.00000380
2 Richness     0.984 0.366  



hist(Alpha_diversity3$Evenness)
hist(Alpha_diversity3$Richness)
hist(Alpha_diversity3$Shannon)
hist(Alpha_diversity3$Reads)

# test for equal variance
levene_test(Alpha_diversity3, rank(Richness) ~ Treatment)
# p = 0.230 - yes
levene_test(Alpha_diversity3, Richness ~ Lysis)
# p = 0.535 - yes
levene_test(Alpha_diversity3, rank(Shannon) ~ Treatment)
# p = 0.439 - yes
levene_test(Alpha_diversity3, Shannon ~ Lysis)
# p = 0.993 - yes
levene_test(Alpha_diversity3, Evenness ~ Treatment)
# p = 0.00562 - no
levene_test(Alpha_diversity3, Evenness ~ Lysis)
# p = 0.393 - yes
levene_test(Alpha_diversity3, Reads ~ Treatment)
# p = 0.185 - yes
levene_test(Alpha_diversity3, Reads ~ Lysis)
# p = 0.347 - yes

# if data are not both normal and equal variance then you can try transforming the data in an attempt to get them to the requirements for parametric tests.
# there are many, but ranking is one of many procedures used to transform data that do not meet the assumptions of normality. A rank-based procedure has been recommended as being robust to non-normal errors, resistant to outliers, and highly efficient for many distributions.
shapiro_test(rank(Alpha_diversity3$Evenness))
# p =  0.0134 - no
shapiro_test(rank(Alpha_diversity3$Reads))
# p = 0.0134 - no
shapiro_test(rank(Alpha_diversity3$Shannon))

lmrich3 <- lm(rank(Richness) ~ Treatment * Lysis, data = Alpha_diversity3)
anova(lmrich)

eta_squared(lmrich3)
#Analysis of Variance Table

#Response: rank(Richness)
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Treatment        6  6956.8  1159.5  8.6405 1.163e-06 ***
#Lysis            1 11934.2 11934.2 88.9346 3.663e-13 ***
#Treatment:Lysis  6  2020.2   336.7  2.5091   0.03201 *  
# Residuals       56  7514.7   134.2                      
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


kruskal_test(data = Alpha_diversity3, Reads ~ Treatment)
kruskal_test(data = Alpha_diversity3, Reads ~ Lysis)
kruskal_test(data = Alpha_diversity3, Richness ~ Treatment)
kruskal_test(data = Alpha_diversity3, Richness ~ Lysis)
kruskal_test(data = Alpha_diversity3, Evenness ~ Treatment)
kruskal_test(data = Alpha_diversity3, Evenness ~ Lysis)
kruskal_test(data = Alpha_diversity3, Shannon ~ Treatment)
kruskal_test(data = Alpha_diversity3, Shannon ~ Lysis)

kruskal_effsize(data = Alpha_diversity3, Reads ~ Treatment)
kruskal_effsize(data = Alpha_diversity3, Reads ~ Lysis)
kruskal_effsize(data = Alpha_diversity3, Richness ~ Treatment)
kruskal_effsize(data = Alpha_diversity3, Richness ~ Lysis)
kruskal_effsize(data = Alpha_diversity3, Evenness ~ Treatment)
kruskal_effsize(data = Alpha_diversity3, Evenness ~ Lysis)
kruskal_effsize(data = Alpha_diversity3, Shannon ~ Treatment)
kruskal_effsize(data = Alpha_diversity3, Shannon ~ Lysis)

wilcox.test(Richness ~ Lysis, data = Alpha_diversity3)

p1 <- ggplot(Alpha_diversity3, aes(x = Lysis, y = Richness, fill = Lysis)) +
  geom_boxplot() + #optional - delete the code in brackets to leave as default or replace with the colour of your choice
  labs(x = "Lysis Method", y = "Richness") + scale_fill_npg()
p1
ggsave("LysisRichness", p1, dpi = 1200, width = 15, height = 6)

p2 <- ggplot(Alpha_diversity3, aes(x = Treatment, y = Reads, fill = Treatment)) +
  geom_boxplot() + #optional - delete the code in brackets to leave as default or replace with the colour of your choice
  labs(x = "Lysis Method", y = "Reads") + scale_fill_npg()
p2