#E.coli-data from Elias vai-4th May 2023
setwd("E:/Others_project/Elias_vai_project")
ec1<-read.csv("E_coli.csv")
#making it tidy
library(dplyr)
library(tidyr)
ec2 <- ec1 %>%
  pivot_longer(cols = 12:31, 
               names_to = "Antibiotics", 
               values_to = "Susceptibility")


#MARI of pathotypes
library(dplyr)
ecM1<-ec1%>%
  select("MARI","pathotype")%>%
  filter(!(pathotype %in% c("Commensal","ETEC")))%>%
  #filter(pathotype !="Commensal" )
  mutate(pathotype = ifelse(pathotype == "aEPEC", "EPEC", pathotype))
   # Test to find out mean difference
      #T-test & wilcox.test() for mean diiference if the group number is two
      # T- test for normal distribution
      #wilcox.test() for non-normal distribution
  #Finding out whether my data is normally distributed
     #p-value from the Shapiro-Wilk test, whether your data is normally distributed or not
     # Shapiro-Wilk test assesses the null hypothesis that your data follows a normal distribution.
     #If the p-value from the test is greater than a chosen significance level (e.g., 0.05), it indicates that the data does not significantly deviate from a normal distribution
shapiro.test(ecM1$MARI)
  # p value =0.594>0.5 ,thus nomal distribution
library(ggplot2)
library(ggpubr)
# Perform t-test and extract p-value
resultM1 <- t.test(MARI ~ pathotype, data = ecM1)
pvalueM1 <- resultM1$p.value

p1 <- ggplot(ecM1, aes(x = pathotype, y = MARI, fill = pathotype)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    comparisons = list(c("EHEC", "EPEC")),
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(a)",
    x = "Pathotype",
    y = "MAR index"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )
p1

               #0.127, no significant difference.

#MARI of phylogroup
library(dplyr)
library(ggplot2)
library(ggpubr)

ecM2 <-ec1%>%
  select("MARI","phylogroup")%>%
  filter(phylogroup %in% c("B1","C"))
resultM2 <- t.test(MARI ~ phylogroup, data = ecM2)
pvalueM2 <- resultM2$p.value

p2 <-ggplot(ecM2, aes(x = phylogroup, y = MARI, fill = phylogroup)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    comparisons = list(c("B1", "C")),
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(b)",
    x = "Phylogroup"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )
p2
                #0.5, no significant difference

#MARI in different season
ecM3 <-ec1%>%
  select("MARI","season")%>%
  filter(season %in% c("Summer","Rainy"))
resultM3 <- t.test(MARI ~ season, data = ecM3)
pvalueM3 <- resultM3$p.value

p3 <-ggplot(ecM3, aes(x = season, y = MARI, fill = season)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    comparisons = list(c("Summer", "Rainy")),
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(C)",
    x = "Season"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p3

# Comparison of MARI Between Pathogenic and Commensal Isolates
ecM5 <- ec1 %>%
  select(MARI, diarrhea) %>%
  mutate(diarrhea = recode(diarrhea,
                           "non-diarrheagenic" = "Commensal",
                           "diarrheagenic" = "Pathogenic"))

# Plot
p4 <-ggplot(ecM5, aes(x = diarrhea, y = MARI, fill = diarrhea)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    method = "t.test",
    comparisons = list(c("Pathogenic", "Commensal")),
    label = "p.format",  # Use this if you decide to show the p-value later
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(d)",
    x = "Isolate Type",
    y = "MARI"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p4



#Comparison of Antibiotic Resistance Genes by Pathotype
library(ggplot2)
library(ggpubr)
ecA1<-ec1%>%
  select("ARG","pathotype")%>%
  filter(!(pathotype %in% c("Commensal","ETEC")))%>%
  #filter(pathotype !="Commensal" )
  mutate(pathotype = ifelse(pathotype == "aEPEC", "EPEC", pathotype))
resultA1 <- wilcox.test(ecA1$ARG ~ ecA1$pathotype,exact= FALSE)
pvalueA1 <- resultA1$p.value

p5 <-ggplot(ecA1, aes(x = pathotype, y = ARG, fill = pathotype)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("EHEC", "EPEC")),
    label = "p.format",
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(e)",
    x = "Pathotype",
    y = "ARG Count"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p5
#0.4,no significant difference



#Frequency of Antibiotic Resistance Genes by Phylogroup
ecA2 <-ec1%>%
  select("ARG","phylogroup")%>%
  filter(phylogroup %in% c("B1","C"))
shapiro.test(ecA2$ARG)
resultA2 <- wilcox.test(ecA2$ARG ~ ecA2$phylogroup,exact= FALSE)
pvalueA2 <- resultA2$p.value

p6 <-ggplot(ecA2, aes(x = phylogroup, y = ARG, fill = phylogroup)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    method = "wilcox.test",
    exact = FALSE,
    comparisons = list(c("B1", "C")),
    label = "p.format",
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(f)",
    x = "Phylogroup",
    y = "ARG Count"
  ) +
  scale_fill_brewer(palette = "Pastel2") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p6

#0.00073, Frequency of ARG is significantly higher than B1



#Frequency of Antibiotic Resistance Genes by Source
ecA4 <-ec1%>%
  select("ARG","source")%>%
  filter(source %in% c("poultry","cow"))
resultA4 <- wilcox.test(ARG ~ source, data = ecA4)
pvalueA4 <- resultA4$p.value

p7 <- ggplot(ecA4, aes(x = source, y = ARG, fill = source)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("cow", "poultry")),
    label = "p.format",
    vjust = 1.5,
    size = 4,
    scientific = FALSE
  ) +
  labs(
    title = "(g)",
    x = "Source",
    y = "ARG Count"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p7
#Frequency of ARG is sigificantly higher in poultry isolates.



#Comparison of ARG Frequency Between Pathogenic and Commensal Isolates
ecA5<-ec1%>%
  select("ARG","diarrhea")%>%
  mutate(diarrhea = recode(diarrhea, "non-diarrheagenic" = "Commensal", "diarrheagenic" = "Pathogenic"))
resultA5 <-  wilcox.test(ecA5$ARG ~ ecA5$diarrhea,exact= FALSE)
pvalueA5 <- resultA5$p.value

p8 <- ggplot(ecA5, aes(x = diarrhea, y = ARG, fill = diarrhea)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Pathogenic", "Commensal")),
    label = "p.format",
    vjust = 1.5,
    size = 4
  ) +
  labs(
    title = "(h)",
    x = "Isolate Type",
    y = "ARG Count"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p8
#0.27, no significant difference.


library(patchwork)

f1 <- (p1 + p2 + p3 + p4 + 
         plot_layout(ncol = 4, widths = rep(1, 4))) /
  (p5 + p6 + p7 + p8 + 
     plot_layout(ncol = 4, widths = rep(1, 4))) +
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(10, 20, 10, 20)) 
f1
# Save as high-resolution TIFF
ggsave("f1.tiff", plot = f1, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure1############################



# ARG in diferent biofilm forming groups
      normality_testA6 <- ec1 %>%
        group_by(biofilm) %>%
        summarize(p_value = shapiro.test(ARG)$p.value)
      bartlett_resultA6 <- bartlett.test(MARI ~ biofilm, data = ec1)
      modelM6<-TukeyHSD(aov(ARG ~ biofilm, data = ec1))






##Summaries based on sample types


##### Percentage of Resistant Isolates by Sample Type

ec3_antibiotic <- ec2 %>%
  # Group and count resistant samples
  group_by(sample_type, Antibiotics, Susceptibility) %>%  # Use correct column name
  summarise(count = n(), .groups = 'drop') %>%
  # Calculate percentage resistance
  group_by(sample_type, Antibiotics) %>%
  mutate(percentage = count / sum(count, na.rm = TRUE) * 100) %>%
  ungroup() %>% 
  # Keep only resistant samples
  filter(Susceptibility == "R") 

# Heatmap (resistance %)
p9 <-ggplot(ec3_antibiotic, aes(x = Antibiotics, y = sample_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add subtle grid lines
  scale_fill_gradient(
    low = "white", 
    high = "red", 
    na.value = "gray90",  # Color for NA values
    limits = c(0, 100)   # Fix scale from 0% to 100%
  ) +
  labs(
    x = "Antibiotic", 
    y = "Sample Type", 
    fill = "Resistance (%)",
    title = "(a)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()  # Remove default grid lines
  ) +
  # Optional: Add text labels if space allows
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")), 
    color = "black", 
    size = 3
  )

p9


#####Mean Multiple Antibiotic Resistance Index (MARI) by Sample Type
ec3_MARI <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and biofilm status
  group_by(sample_type) %>%
  summarise(
    mean_MARI = mean(MARI, na.rm = TRUE))
  
p10 <-ggplot(ec3_MARI, aes(x = 1, y = sample_type, fill = mean_MARI)) +  # x=1 creates a single column
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Mean MARI",
    limits = c(0, max(ec3_MARI$mean_MARI, na.rm = TRUE))  # Auto-scale upper limit
  ) +
  labs(
    x = NULL,  # Remove x-axis label
    y = "Sample Type",
    title = "(b)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = round(mean_MARI, 2)),  # Show values with 2 decimal places
    color = "black",
    size = 4
  )

p10


##### mean Antibiotic Resistance Gene frequency
ec3_ARG <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and biofilm status
  group_by(sample_type) %>%
  summarise(
    mean_ARG = mean(ARG, na.rm = TRUE))

p11 <-ggplot(ec3_ARG, aes(x = 1, y = sample_type, fill = mean_ARG)) +  # x=1 creates a single column
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Mean ARG",
    limits = c(0, max(ec3_ARG$mean_ARG, na.rm = TRUE))  # Auto-scale upper limit
  ) +
  labs(
    x = NULL,  # Remove x-axis label
    y = "Sample Type",
    title = "(c)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = round(mean_ARG, 2)),  # Show values with 2 decimal places
    color = "black",
    size = 4
  )

p11

library(patchwork)

f2 <- (p9 + 
         plot_layout(ncol = 1)) /
  (p10 + p11 + 
     plot_layout(ncol = 2, widths = c(1, 1))) +
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(15, 20, 15, 20))  # top, right, bottom, left
f2
# Save as high-resolution TIFF
ggsave("f2.tiff", plot = f2, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure2############################




library(dplyr)
library(ggplot2)

##### Biofilm Distribution by Sample Type
ec3_biofilm <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and biofilm status
  group_by(sample_type, biofilm) %>%
  # Count occurrences
  count() %>%
  # Calculate percentages within each sample_type
  group_by(sample_type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() 

p12 <- ggplot(ec3_biofilm, aes(x = biofilm, y = sample_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 100),
    na.value = "gray90"
  ) +
  labs(
    x = "Biofilm Strength",
    y = "Sample Type",
    title = "(a)",
    fill = "Percentage (%)"
  ) +
  theme_test(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")),
    color = "black",
    size = 3
  )

p12



##### Pathotypes Distribution by Sample Type
ec3_pathotypes <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and pathotypes status
  group_by(sample_type, pathotype) %>%
  # Count occurrences
  count() %>%
  # Calculate percentages within each sample_type
  group_by(sample_type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() 

p13 <- ggplot(ec3_pathotypes, aes(x = pathotype, y = sample_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 100),
    na.value = "gray90"
  ) +
  labs(
    x = "Pathotypes",
    y = "Sample Type",
    title = "(b)",
    fill = "Percentage (%)"
  ) +
  theme_test(base_size = 9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")),
    color = "black",
    size = 3
  )

p13



##### Phylogroup Distribution by Sample Type
ec3_phylogroup <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and phylogroup status
  group_by(sample_type, phylogroup) %>%
  # Count occurrences
  count() %>%
  # Calculate percentages within each sample_type
  group_by(sample_type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() 

p14 <- ggplot(ec3_phylogroup, aes(x = phylogroup, y = sample_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 100),
    na.value = "gray90"
  ) +
  labs(
    x = "Phylogroup",
    y = "Sample Type",
    title = "(c)",
    fill = "Percentage (%)"
  ) +
  theme_test(base_size = 9) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")),
    color = "black",
    size = 3
  )

p14



##### pathogenicity Distribution by Sample Type
ec3_diarrhea <- ec2 %>%
  # First filter for unique sample_ID (keeps first occurrence if duplicates exist)
  distinct(sample_ID, .keep_all = TRUE) %>%
  # Group by sample_type and diarrhea status
  group_by(sample_type, diarrhea) %>%
  # Count occurrences
  count() %>%
  # Calculate percentages within each sample_type
  group_by(sample_type) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup() 

p15 <- ggplot(ec3_diarrhea, aes(x = diarrhea, y = sample_type, fill = percentage)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 100),
    na.value = "gray90"
  ) +
  labs(
    x = "Isolate types",
    y = "Sample Type",
    title = "(d)",
    fill = "Percentage (%)"
  ) +
  theme_test(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  geom_text(
    aes(label = paste0(round(percentage, 1), "%")),
    color = "black",
    size = 3
  )

p15

library(patchwork)

f3 <- (p12 + p13 + 
         plot_layout(ncol = 2, widths = c(1, 1))) /
  (p14 + p15 + 
     plot_layout(ncol = 2, widths = c(1, 1))) +
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(15, 20, 15, 20))  # top, right, bottom, left

f3

# Save as high-resolution TIFF
ggsave("f3.tiff", plot = f3, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure3############################





# Function to create a labeled bar chart

make_bar_plot <- function(data, column_name, title) {
  df <- data %>%
    count(!!sym(column_name)) %>%
    mutate(
      perc = round(n / sum(n) * 100, 1),
      label = paste0(perc, "%")
    )
  
  max_perc <- max(df$perc)
  
  ggplot(df, aes(x = reorder(!!sym(column_name), -perc), y = perc, fill = !!sym(column_name))) +
    geom_col(width = 0.7) +
    geom_text(aes(label = label), vjust = -0.5, size = 4) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(
      limits = c(0, max_perc + 10),  # Add buffer to top
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
  labs(title = title, y = "Percentage")
}


### Pathotypes Distribution
p16 <- make_bar_plot(ec1, "pathotype", "(a)") +
  labs(x = "Pathotypes") +
  theme(axis.title.x = element_text(size = 14, face = "bold"))

### Phylogroup Distribution  
p17 <- make_bar_plot(ec1, "phylogroup", "(b)") +
  theme(axis.title.y = element_blank()) +
  labs(x = "Phylogroups") +
  theme(axis.title.x = element_text(size = 14, face = "bold"))

### Biofilm Formation
p18 <- make_bar_plot(ec1, "biofilm", "(c)") +
  theme(axis.title.y = element_blank()) +
  labs(x = "Biofilm classification") +
  theme(axis.title.x = element_text(size = 14, face = "bold"))

### Diarrheagenic potential
p19 <- make_bar_plot(ec1, "diarrhea", "(d)") +
  labs(x = "Diarrheagenic potential") +
  theme(axis.title.x = element_text(size = 14, face = "bold"))


p16
p17
p18
p19



### Biofilm gene percentage

file2 <-read.csv("file2.csv")

file2.1 <- file2 %>%
  filter(categories == "Biofilm") %>%
  mutate(label = paste0(round(percentage, 1), "%"))

max_perc_2.1 <- max(file2.1$percentage)

p20 <- ggplot(file2.1, aes(x = gene, y = percentage, fill = gene)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    limits = c(0, max_perc_2.1 + 10),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold.italic", size = 12),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(title = "(e)", x = "Biofilm genes")

p20





### Antibiotic Resistance Gene percentage
library(colorspace)

file2.2 <- file2 %>%
  filter(categories == "Antibiotic Resistance") %>%
  mutate(label = paste0(round(percentage, 1), "%"))

max_perc_2.2 <- max(file2.2$percentage)
my_colors <- qualitative_hcl(14, palette = "Set2")

p21 <- ggplot(file2.2, aes(x = gene, y = percentage, fill = gene)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), vjust = -0.5, size = 4) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(
    limits = c(0, max_perc_2.2 + 10),
    expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold.italic", size = 12),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  labs(title = "(f)", x = "Antibiotic Resistance Genes")

p21


library(patchwork)

# Combine plots with patchwork# Combine plots with pageneric.skeleton()tchwork
f4 <- (p16 + p17 + p18 + 
         plot_layout(ncol = 3, widths = c(1, 1, 1))) /
  (p19 + p20 + p21 + 
     plot_layout(ncol = 3, widths = c(1, 1, 1))) +
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(15, 20, 15, 20))  # top, right, bottom, left

f4


# Save as high-resolution TIFF
ggsave("f4.tiff", plot = f4, width = 16, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure4############################




###### Antibiotic Resistance Percentage ##### 
ecB6<-ec2%>%
  group_by(Antibiotics)%>%
  count(Susceptibility)%>%
  mutate(prob = n/sum(n)*100)
library(ggplot2)
library(ggpubr)

p22 <-ggplot(ecB6, aes(x = reorder(Antibiotics, -prob), y = prob, fill = Susceptibility)) +
  geom_col(position = "stack", color = "black", width = 0.75) +
  scale_fill_manual(
    values = c("S" = "#a6cee3", "I" = "#1f78b4", "R" = "#b2df8a")
  ) +
  labs(
    title = "Antibiotic Resistance Profile",
    x = "Antibiotics",
    y = "Isolates",
    fill = "Susceptibility"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))

# Save as high-resolution TIFF
ggsave("f5.tiff", plot = p22, width = 12, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure5############################





##### meanOD bassed analysis######################

#Comparison of meanOD Between EHEC and EPEC Pathotypes
library(ggplot2)
library(ggpubr)
ecB1<-ec1%>%
  select("meanOD","pathotype")%>%
  filter(!(pathotype %in% c("Commensal","ETEC")))%>%
  #filter(pathotype !="Commensal" )
  mutate(pathotype = ifelse(pathotype == "aEPEC", "EPEC", pathotype))
resultb1 <- wilcox.test(meanOD ~ pathotype, data = ecB1)
pvalueb1 <- resultb1$p.value
ecB1$pathotype <- factor(ecB1$pathotype, levels = c("EHEC", "EPEC"))

p23 <-ggplot(ecB1, aes(x = pathotype, y = meanOD, fill = pathotype)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2, alpha = 0.9) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "black") +
  stat_compare_means(
    label = "p.format",
    method = "wilcox.test",
    comparisons = list(c("EHEC", "EPEC")),
    size = 4,
    vjust = 1.5
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "(a)",
    x = "Pathotype",
    y = "Mean OD (Biofilm Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p23
#0.98, no significant difference





#Biofilm-Forming Ability by Phylogroup
library(dplyr)
library(ggplot2)
library(ggpubr)
ecB2 <-ec1%>%
  select("meanOD","phylogroup")%>%
  filter(phylogroup %in% c("B1","C"))
shapiro.test(ecB2$meanOD)
resultB2 <- t.test(ecB2$meanOD ~ ecB2$phylogroup,exact= FALSE)
pvalueB2 <- resultB2$p.value
ecB2$phylogroup <- factor(ecB2$phylogroup, levels = c("B1", "C"))

p24 <-ggplot(ecB2, aes(x = phylogroup, y = meanOD, fill = phylogroup)) +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.15), size = 1.8, shape = 21, stroke = 0.2, alpha = 0.6) +
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    comparisons = list(c("B1", "C")),
    size = 4,
    vjust = 1.4
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "(b)",
    x = "Phylogroup",
    y = "Mean OD (Biofilm Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p24 
#0.049, meanOD is higher in B1 phylogroup.





#Biofilm-Forming Ability by Phylogroup
ecB4 <-ec1%>%
  select("meanOD","source")%>%
  filter(source %in% c("poultry","cow"))
resultM4 <- wilcox.test(meanOD ~ source, data = ecB4)
pvalueM4 <- resultM4$p.value
#no significant meanOD difference

#p_label <- paste0("p = ", signif(pvalueM4, 3))

# Create the plot
p25 <-ggplot(ecB4, aes(x = source, y = meanOD, fill = source)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  stat_compare_means(
    label = "p.format",
    method = "t.test",
    comparisons = list(c("poultry", "cow")),
    size = 4,
    vjust = 1.4
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "(c)",
    x = "Source",
    y = "Mean OD (Biofilm Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p25




#Biofilm-Forming Ability: Pathogens vs Commensals
ecB5<-ec1%>%
  select("meanOD","diarrhea")%>%
  mutate(diarrhea = recode(diarrhea, "non-diarrheagenic" = "Commensal", "diarrheagenic" = "Pathogenic"))

resultB5 <- wilcox.test(meanOD ~ diarrhea, data = ecB5)
pvalueB5 <- resultB5$p.value
#0.23, no significant difference

ecB5$diarrhea <- factor(ecB5$diarrhea, levels = c("Pathogenic", "Commensal"))

p26 <-ggplot(ecB5, aes(x = diarrhea, y = meanOD, fill = diarrhea)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.6, shape = 21, stroke = 0.3, alpha = 0.6) +
  stat_compare_means(
    label = "p.format",   # You can also use "p.signif" for significance stars
    method = "wilcox.test",
    exact = FALSE,
    comparisons = list(c("Pathogenic", "Commensal")),
    size = 4,
    vjust = 1.5
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "(d)",
    x = "Isolate Type",
    y = "Mean OD (Biofilm Index)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

p26

# Combine plots with patchwork# Combine plots with pageneric.skeleton()tchwork
f6 <- (p23 + p24 + 
         plot_layout(ncol = 2, widths = c(1, 1))) /
  (p25 + p26 + 
     plot_layout(ncol = 2, widths = c(1, 1))) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(15, 20, 15, 20))  # top, right, bottom, left

f6

# Save as high-resolution TIFF
ggsave("f6.tiff", plot = f6, width = 12, height = 8, dpi = 300, units = "in", device = "tiff")
##################figure4############################





           #since three groups, need to perform ANOVA
           #Need to check two assumption for one way ANOVa
           #first assumption: Normality
           normality_testM3 <- ec1 %>%
               group_by(season) %>%
               summarize(p_value = shapiro.test(MARI)$p.value)
           #p values greater thsn 0.05, do not deviate significantly from a normal distribution.
            #2nd assumption:Equal variance
           bartlett_resultM3 <- bartlett.test(MARI ~ season, data = ec1)
           #p value greater than 0.05, variances are equal across the groups. 
       #performing ANOVA
modelM3<-TukeyHSD(aov(MARI ~ season, data = ec1))
#summer vs Rainy, significant difference in MARI

ggplot(ec1, aes(MARI,season))+
geom_boxplot() +
  geom_signif(comparisons = list(c("Summer", "Rainy"), c("Summer", "Late autumn"), c("Rainy", "Late autumn")),
              test = "aov", map_signif_level = TRUE, textsize = 3.5, vjust = -0.5, hjust = 0.5) +
  geom_text(aes(x = 1.5, y = max(ec1$MARI) + 0.2, 
                label = paste0("p adj = ", format(modelM3$p_adj, scientific = FALSE))),
            hjust = 0.5, size = 4, fontface = "bold") +
  labs(title = "Comparison of MARI by Season", x = "Season", y = "MARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
library(ggsignif)

bartlett_result <- bartlett.test(MARI ~ season, data = ec1)

# Print the test result
print(bartlett_result)
# Perform ANOVA
modelM2 <- aov(MARI ~ season, data = ec1)


# Perform Tukey HSD test
tukey_result <- TukeyHSD(modelM2)

# Extract p-values
p_adj_values <- tukey_result$p.adj[, "p adj"]

# Create the plot
ggplot(ec1, aes(season, MARI)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("Summer", "Rainy"), c("Summer", "Late autumn"), c("Rainy", "Late autumn")),
              test = "aov", map_signif_level = TRUE, textsize = 3.5, vjust = -0.5, hjust = 0.5) +
  geom_text(aes(x = 1.5, y = max(ec1$MARI) + 0.2, 
                label = paste0("p adj = ", format(p_adj_values, scientific = FALSE))),
            hjust = 0.5, size = 4, fontface = "bold") +
  labs(title = "Comparison of MARI by Season", x = "Season", y = "MARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


    # Create the plot
ggplot(ec1, aes(season, MARI)) +
  geom_boxplot() +
  stat_compare_means(label = "", 
                     method = "anova", label.y = 1.5) +
  annotate("text", x = 2, y = max(ecM1$MARI), 
           label = paste0("ANOVA, p-value = ", format(p_adj_values, scientific = FALSE, digits = 3)), 
           hjust = 0, vjust = 1.5) +
  labs(title = "Comparison of MARI by Season", x = "Season", y = "MARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

library(ggplot2)
library(ggsignif)
library(dplyr)

  # Perform ANOVA test
resultM3 <- aov(MARI ~ season, data = ec1)
pvalueM3 <- summary(resultM3)$"Pr(>F)"[1]
summary(resultM3)
TukeyHSD(resultM3)
  # Create the plot

ggplot(ec1, aes(season, MARI)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("Summer", "Rainy"), c("Summer", "Late autumn"),c("Rainy","Late autumn")), 
              
              step_increase = 0.05,
              vjust = -0.75,
              annotations = paste0("p = ", format(pvalueM3, scientific = FALSE, digits = 3)),
              y_position = max(ecM1$MARI) * 1.1) +
  labs(title = "Comparison of MARI by Season", x = "Season", y = "MARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))



#MARI in different source
library(dplyr)
library(ggplot2)
library(ggpubr)
ecM4 <-ec1%>%
  select("MARI","source")%>%
  filter(source %in% c("poultry","cow"))
resultM4 <- t.test(MARI ~ source, data = ecM4)
pvalueM4 <- resultM4$p.value
                 #0.072--pvalue>0.05, no significant difference

#MARI commensal vs pathogens
ecM5<-ec1%>%
  select("MARI","diarrhea")%>%
  mutate(diarrhea = recode(diarrhea, "non-diarrheagenic" = "Commensal", "diarrheagenic" = "Pathogenic"))
resultM5 <- t.test(MARI ~ diarrhea, data = ecM5)
pvalueM5 <- resultM5$p.value
ggplot(ecM5, aes(diarrhea, MARI)) +
  geom_boxplot() +
  stat_compare_means(label = "", method = "t.test", comparisons = list(c("Pathogenic", "Commensal"))) +
  labs(title = "Comparison of MARI between Pathogens and Commensals", x = "Pathogenic", y = "MARI") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
                    #0.36, no significant difference

#MARI in different biofilm forming groups
normality_testM6 <- ec1 %>%
  group_by(biofilm) %>%
  summarize(p_value = shapiro.test(MARI)$p.value)
bartlett_resultM6 <- bartlett.test(MARI ~ biofilm, data = ec1)
modelM6<-TukeyHSD(aov(MARI ~ biofilm, data = ec1))
                   #no significant difference.





#meanOD in different season
normality_testB3 <- ec1 %>%
  group_by(season) %>%
  summarize(p_value = shapiro.test(meanOD)$p.value)
        #each group is non-normally distributed, perform Kruskal-Wallis Test
install.packages("conover.test")
library(conover.test)
modelB3<-conover.test(kruskal.test(meanOD ~ season, data = ec1))
      

#meanOD in diferent antibiotic resistance groups
library(dplyr)


  

my_anova <- aov(meanOD ~ Susceptibiliy,data = ec2)
anova_summary <- summary(my_anova)
tukey_results <- TukeyHSD(my_anova)
pvalues <- anova_summary$p.value







