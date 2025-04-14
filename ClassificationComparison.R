#Load in Packages
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
library(ggsci)

#Load in data set's

silvablast <- read_csv("Trial1 silva BLAST taxonomy.csv")
silvabay <- read_csv("Trail1 silva Bay taxonomy.csv")
spikeblast <- read_csv("Trail1 spike BLAST taxonomy.csv")
spikebay <- read_csv("Trial1 Spike Bay taxonomy.csv")


#Add method and classifier column to each data set

silvablast$method <- "blast"
silvabay$method <- "bay"
spikeblast$method <- "blast"
spikebay$method <- "bay"

silvablast$classifier <- "silva"
silvabay$classifier <- "silva"
spikeblast$classifier <- "spike"
spikebay$classifier <- "spike"

#Combine into one data set

data <- bind_rows(silvablast, silvabay, spikeblast, spikebay)

#Do jaccard dissimilarity for presence/absence of taxa at genus level

silva_data <- data %>%
  filter(classifier == "silva")

taxa_Domain_silva_same <- silva_data %>%
  select(`Feature ID`, Domain, method) %>%
  spread(key = method, value = Domain) %>%
  mutate(same_domain = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_domain)

taxa_phylum_silva_same <- silva_data %>%
  select(`Feature ID`, Phylum, method) %>%
  spread(key = method, value = Phylum) %>%
  mutate(same_phylum = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_phylum)

taxa_class_silva_same <- silva_data %>%
  select(`Feature ID`, Class, method) %>%
  spread(key = method, value = Class) %>%
  mutate(same_class = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_class)

taxa_order_silva_same <- silva_data %>%
  select(`Feature ID`, Order, method) %>%
  spread(key = method, value = Order) %>%
  mutate(same_order = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_order)

taxa_family_silva_same <- silva_data %>%
  select(`Feature ID`, Family, method) %>%
  spread(key = method, value = Family) %>%
  mutate(same_family = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_family)

taxa_genus_silva_same <- silva_data %>%
  select(`Feature ID`, Genus, method) %>%
  spread(key = method, value = Genus) %>%
  mutate(same_genus = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_genus)


taxa_species_silva_same <- silva_data %>%
  select(`Feature ID`, Species, method) %>%
  spread(key = method, value = Species) %>%
  mutate(same_species = if_else(is.na(bay) & is.na(blast), NA,
                                if_else(is.na(bay) | is.na(blast), 0,
                                        if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_species)

combined_taxa_data <- taxa_Domain_silva_same %>%
  left_join(taxa_phylum_silva_same, by = 'Feature ID') %>%
  left_join(taxa_class_silva_same, by = 'Feature ID') %>%
  left_join(taxa_order_silva_same, by = 'Feature ID') %>%
  left_join(taxa_family_silva_same, by = 'Feature ID') %>%
  left_join(taxa_genus_silva_same, by = 'Feature ID') %>%
  left_join(taxa_species_silva_same, by = 'Feature ID')



jaccard_sim <- colMeans(combined_taxa_data[, -1], na.rm = TRUE)

jaccard_table <- data.frame(Mean_value = jaccard_sim)

rownames(jaccard_table) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

jaccard_table <- data.frame("Taxanomic_rank" = rownames(jaccard_table), jaccard_table)

rownames(jaccard_table) <- NULL

silva_only <- ggplot(jaccard_table, aes(x = Taxanomic_rank, y = mean_value, group = 1)) +
  geom_line() +
    geom_point() +
  labs(title = "Similarity between Blast and Bay classification at each taxanomic rank", x = "Taxanomic Rank", y = "Jaccard Similarity Score")




spike_data <- data %>%
  filter(classifier == "spike")

taxa_Domain_spike_same <- spike_data %>%
  select(`Feature ID`, Domain, method) %>%
  spread(key = method, value = Domain) %>%
  mutate(same_domain = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_domain)

taxa_phylum_spike_same <- spike_data %>%
  select(`Feature ID`, Phylum, method) %>%
  spread(key = method, value = Phylum) %>%
  mutate(same_phylum = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_phylum)

taxa_class_spike_same <- spike_data %>%
  select(`Feature ID`, Class, method) %>%
  spread(key = method, value = Class) %>%
  mutate(same_class = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_class)

taxa_order_spike_same <- spike_data %>%
  select(`Feature ID`, Order, method) %>%
  spread(key = method, value = Order) %>%
  mutate(same_order = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_order)

taxa_family_spike_same <- spike_data %>%
  select(`Feature ID`, Family, method) %>%
  spread(key = method, value = Family) %>%
  mutate(same_family = if_else(is.na(bay) & is.na(blast), NA,
                               if_else(is.na(bay) | is.na(blast), 0,
                                       if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_family)

taxa_genus_spike_same <- spike_data %>%
  select(`Feature ID`, Genus, method) %>%
  spread(key = method, value = Genus) %>%
  mutate(same_genus = if_else(is.na(bay) & is.na(blast), NA,
                              if_else(is.na(bay) | is.na(blast), 0,
                                      if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_genus)


taxa_species_spike_same <- spike_data %>%
  select(`Feature ID`, Species, method) %>%
  spread(key = method, value = Species) %>%
  mutate(same_species = if_else(is.na(bay) & is.na(blast), NA,
                                if_else(is.na(bay) | is.na(blast), 0,
                                        if_else(bay == blast, 1, 0)))) %>%
  select('Feature ID', same_species)

combined_taxa_data_spike <- taxa_Domain_spike_same %>%
  left_join(taxa_phylum_spike_same, by = 'Feature ID') %>%
  left_join(taxa_class_spike_same, by = 'Feature ID') %>%
  left_join(taxa_order_spike_same, by = 'Feature ID') %>%
  left_join(taxa_family_spike_same, by = 'Feature ID') %>%
  left_join(taxa_genus_spike_same, by = 'Feature ID') %>%
  left_join(taxa_species_spike_same, by = 'Feature ID')

jaccard_sim_spike <- colMeans(combined_taxa_data_spike[, -1], na.rm = TRUE)


jaccard_table_spike <- data.frame(mean_value = jaccard_sim_spike)

rownames(jaccard_table_spike) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

jaccard_table_spike <- data.frame("Taxanomic_rank" = rownames(jaccard_table_spike), jaccard_table_spike)

rownames(jaccard_table_spike) <- NULL

total_data_spike <- colSums(combined_taxa_data_spike[, -1], na.rm = TRUE)
plot_data <- bind_cols(jaccard_table, jaccard_table_spike)

colnames(plot_data) <- c("Taxanomic_rank", "Silva", "Spike")

plot_longer <- pivot_longer(plot_data, cols = 2:3, names_to = "Classifier", values_to = "Jaccard_Similarity_Score")

plot_longer$Taxanomic_rank <- factor(plot_longer$Taxanomic_rank, levels = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
                                     
write.csv(plot_data, file = "plot_data_all.csv", row.names = FALSE)

JaccardPlot <- ggplot(plot_longer, aes(x = Taxanomic_rank, y = Jaccard_Similarity_Score, group = Classifier, color = Classifier)) +
  geom_line() +
  geom_point() +
  scale_colour_npg() +
  labs(x = "Taxanomic Rank", y = "Jaccard Similarity Score") +
  theme(legend.background = element_blank(), panel.background = element_blank(), legend.position = "inside")
JaccardPlot
JaccardPlot <- reposition_legend(JaccardPlot, 'top right') 
ggsave("similarityplot.png", JaccardPlot, width = 8, height = 5, unit = "in", dpi = 600)

theme_set(theme_classic())

write.csv(jaccard_table, file = "jaccard_silva.csv", row.names = FALSE)
write.csv(jaccard_table_spike, file = "jaccard_spike.csv", row.names = FALSE)

jaccard_chi_silva <- read.csv("jaccard_silva.csv", row.names = 1)
jaccard_chi_spike <- read.csv("jaccard_spike.csv", row.names = 1)

chi_result_silva <-chisq.test(jaccard_chi_silva)
print(chi_result_silva)

#Pearson's Chi-squared test

#data:  jaccard_chi_silva
#X-squared = 1082.3, df = 6, p-value < 2.2e-16

chi_result_spike <-chisq.test(jaccard_chi_spike)
print(chi_result_spike)

#Pearson's Chi-squared test

#data:  jaccard_chi_spike
#X-squared = 2355.8, df = 6, p-value < 2.2e-16

salmonella_chi <- read.csv("salmonalladata.csv", row.names = 1)

chi_result_salmonella <-chisq.test(salmonella_chi)
print(chi_result_salmonella)

#Pearson's Chi-squared test

#data:  salmonella_chi
#X-squared = 241564, df = 89, p-value < 2.2e-16

salmonellab_chi <- read.csv("salmonellablastdata.csv", row.names = 1)

salmonellab_chi <- salmonellab_chi + 0.5

chi_result_salmonellab <-chisq.test(salmonellab_chi)
print(chi_result_salmonellab)

#Pearson's Chi-squared test

#data:  salmonellab_chi
##X-squared = 51.949, df = 89, p-value = 0.9994

listeria_chi <- read.csv("Listeriaspike.csv", row.names = 1)

chi_result_listeria <-chisq.test(listeria_chi)
print(chi_result_listeria)

#Pearson's Chi-squared test

#data:  listeria_chi
#X-squared = 2337.2, df = 89, p-value < 2.2e-16

listeriab_chi <- read.csv("Listeriablast.csv", row.names = 1)

listeriab_chi <- listeriab_chi + 0.5

chi_result_listeriab <-chisq.test(listeriab_chi)
print(chi_result_listeriab)


#Pearson's Chi-squared test

#data:  listeriab_chi
#X-squared = 0, df = 89, p-value = 1

campyb_chi <- read.csv("campyblast.csv", row.names = 1)

campyb_chi <- campyb_chi + 0.5

chi_result_campyb<-chisq.test(campyb_chi)
print(chi_result_campyb)

#Pearson's Chi-squared test

#data:  campyb_chi
#X-squared = 0, df = 89, p-value = 1

campyspike_chi <- read.csv("campyspike.csv", row.names = 1)

campyspike_chi <- campyspike_chi + 0.5

chi_result_campyspike <-chisq.test(campyspike_chi)
print(chi_result_campyspike)

#Pearson's Chi-squared test

#data:  campyspike_chi
#X-squared = 35.855, df = 89, p-value = 1

campy_chi <- read.csv("campycomparison.csv", row.names = 1)

campy_chi <- campy_chi + 0.5

chi_result_campy <-chisq.test(campy_chi)
print(chi_result_campy)

#Pearson's Chi-squared test

#data:  campy_chi
#X-squared = 75.546, df = 267, p-value = 1
