### let's set our working directory, where my datasets are kept
setwd("/home/anik/genomics/Projects/Hand_rinse")

### You can check the current directory
getwd()

### Uploading the dataset 
A1 <- read.csv("Arabi_Sampling.csv", check.names = FALSE)

### Cleaning the dataset (undetermined = undetected, any value above 35 transformed to undetected)
A2 <- A1
A2[A2 == "Undetermined" | A2 == "undetermined"] <- "undetected"


### Transforming the dataset for categorical comparison
### transforming CT value higher than 35 as undetected using dplyr
library(dplyr)   

A2 <- A2 %>%
  mutate(across(10:18, ~ if_else(. > 35, "undetected", as.character(.)))) %>%
  mutate(across(10:18, ~ {
    ifelse(. == "undetected", "undetected",
           format(round(as.numeric(.), 2), nsmall = 2))
  }))

summary(A2$mTEC)





##### P1: First Plot: Visualization of detected and undetected (Sample Type)
      
      # We will need to transform the CT values to compare detected and undetected
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)

A3 <- A2 %>%
  mutate(across(10:18, ~ if_else(. <= 35, "detected", as.character(.))))

A3_tidy <- A3 %>%
  pivot_longer(
    cols = 10:18,
    names_to = "pathogen",
    values_to = "status"
  )

A3_group <- A3_tidy %>% 
  group_by(Types, Time, pathogen, status) %>%
  summarise(
    count = n(),
    .groups = "drop_last"  # Important for next step
  ) %>%
  mutate(
    percentage = count / sum(count) * 100
  ) %>%
  ungroup() %>% 
  mutate(Time = factor(Time, levels = c("Before", "After")))  # Or your exact time labels


p1 <-ggplot(A3_group %>% filter(status == "detected"),  # Keep only "detected" status
       aes(x = pathogen, y = percentage, fill = Time)) +  # Fill by Time, not status
  geom_col(position = "dodge") +  # Dodge by Time
  facet_grid(Types ~ .) +         # Vertical facets for Types
  scale_fill_viridis_d(option = 'magma') +
  labs(title = "(A)",
       x = "Gene",
       y = "Detected (%)",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold.italic"),
        axis.text.y = element_text(vjust = 1, size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))


p1  ##Plot1





##### P2: Second Plot: Visualization of detected and undetected (Household ID)


          ## We can start with A3_tidy from P1

A4_group <- A3_tidy %>% 
  group_by(HID, Time, pathogen, status) %>%
  summarise(
    count = n(),
    .groups = "drop_last"  # Important for next step
  ) %>%
  mutate(
    percentage = count / sum(count) * 100
  ) %>%
  ungroup() %>% 
  mutate(Time = factor(Time, levels = c("Before", "After"))) %>%  # Or your exact time labels
  mutate(HID = factor(HID, levels = c('H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10')))


p2 <-ggplot(A4_group %>% filter(status == "detected"),  # Keep only "detected" status
       aes(x = pathogen, y = percentage, fill = Time)) +  # Fill by Time, not status
  geom_col(position = "dodge") +  # Dodge by Time
  facet_grid(HID ~ .) +         # Vertical facets for HID
  scale_fill_viridis_d(option = 'magma') +
  labs(title = "(B)",
       x = "Gene",
       y = "Detected (%)",
       fill = "Time Point") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold.italic"),
        axis.text.y = element_text(vjust = 1, size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))


p2  ##Plot2

figure_1 <- p1 + p2 + 
  plot_layout(guides = "collect")

# Save as high-resolution TIFF
ggsave("Figure_3.tiff", plot = figure_1, width = 14, height = 10, dpi = 300, units = "in", device = "tiff")

figure_1 ##########################################################################################




###############################################################--Figure2--###############
##### P3: third Plot: CT values comparison

          ## WIll have to start from A2
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(rstatix)

A5 <- A2 %>% 
  mutate(across(10:18, ~ {
    ifelse(. == "undetected", 40, as.numeric(.))  # Convert "undetected" to 40, others to numeric
  }))

write.csv(A5, file = "A5_data.csv")
 
          ## Need to convert to paired data, before that conversion to longer format (tidy) is necessary

A5_tidy <- A5 %>%
  pivot_longer(
    cols = 10:18,
    names_to = "Targets",
    values_to = "CT_value"
  ) %>% 
  select(-"Number") 

A5_paired <- A5_tidy %>%
  # Clean Time column just in case (e.g., trim spaces, fix capitalization)
  mutate(Time = str_trim(str_to_title(Time))) %>%
  
  # Pivot to wide format using multiple value columns
  pivot_wider(
    id_cols = c(HID, Types, Date, Targets),  # keep these constant
    names_from = Time,
    values_from = c(CT_value, mTEC),
    names_glue = "{.value}_{Time}"                 # final column names: CT_value_Before, etc.
  ) 

  
A5_paired_filtered <- A5_paired %>% ##### Removing events without pair (Food only before)
  # Keep only rows with complete CT_value data
  filter(!is.na(CT_value_Before) & !is.na(CT_value_After))

A5_paired_filtered_new_columns <- A5_paired_filtered %>% 
  mutate(CT_diff = CT_value_Before - CT_value_After) %>% 
  mutate(mTEC_diff = mTEC_Before - mTEC_After)



A5_CT_within_gene <- A5_paired_filtered %>%
  group_by(Targets) %>%
  summarise(
    p_value = wilcox.test(CT_value_Before, CT_value_After, paired = TRUE)$p.value,
    n = n(),
    .groups = "drop"
  )

A5_paired_filtered_new_columns_paired <- A5_paired_filtered_new_columns %>%
  mutate(Pair_ID = paste(HID, Types, sep = "_"))

a5_CT_friedman_test <- A5_paired_filtered_new_columns_paired %>%
  friedman_test(CT_diff ~ Targets | Pair_ID)

print(a5_CT_friedman_test)

a5_CT_wilcox_test <- A5_paired_filtered_new_columns %>% 
  wilcox_test(CT_diff~ Targets, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Targets")


# Calculate mean and median for each target
a5_summary_stats <- A5_paired_filtered_new_columns %>%
  group_by(Targets) %>%
  summarise(
    Mean = mean(CT_diff, na.rm = TRUE),
    Median = median(CT_diff, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -10.5, 10.5),
    Median_vjust = ifelse(Median > Mean, -10.5, 10.5)
  )

# Plot
p3 <- ggplot(A5_paired_filtered_new_columns, aes(x = Targets, y = CT_diff)) +
  geom_violin(aes(fill = Targets)) +
  scale_fill_viridis_d() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_test() +
  stat_pvalue_manual(a5_CT_wilcox_test, label = "p.signif", hide.ns = TRUE) + 
  labs(title = "(A)",
       y = "Difference in CT values (Before - After)", x = "Targets") +
  theme(
  axis.text.x = element_text(
    face = "bold.italic", hjust = 1, 
    margin = ggplot2::margin(b = 5), 
    size = 12, angle = 90,
  ),
  plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5),
  axis.text.y = element_text(size = 12, face = "bold"), 
  axis.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 14, face = "bold"),
  panel.border = element_rect(linewidth = 1.2),
  legend.position = "none"
) +
  
  # Add mean and median points
  geom_point(data = a5_summary_stats, aes(x = Targets, y = Mean), 
             color = "blue", size = 4, shape = 18) +
  geom_point(data = a5_summary_stats, aes(x = Targets, y = Median), 
             color = "red", size = 4, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = a5_summary_stats, aes(x = Targets, y = Mean, 
                                      label = paste0("Mean: ", round(Mean, 2)), 
                                      vjust = Mean_vjust), 
            color = "blue", size = 4) +
  geom_text(data = a5_summary_stats, aes(x = Targets, y = Median, 
                                      label = paste0("Median: ", round(Median, 2)), 
                                      vjust = Median_vjust), 
            color = "red", size = 4) +
  
  # Annotate p-value and number of samples
  geom_text(
    data = A5_CT_within_gene,
    aes(x = Targets, y = max(A5_paired_filtered_new_columns$CT_diff, na.rm = TRUE) + 6,
        label = paste0("p = ", signif(p_value, 3), "\n", "n = ", n)),
    size = 5
  )

p3    #### Plot 3






##### P4: fourth Plot: Colony difference among households

# Select relevant columns
A5_MTEC <- A5 %>%
  filter(Types != "Food") %>%  # Note the quotes around "Food"
  select(Time, mTEC) 

# 2. Create the plot
A5_MTEC_wilcox <- A5_MTEC %>% 
  wilcox_test(mTEC~ Time, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Time")

A5_MTEC_summary_stats <- A5_MTEC%>%
  group_by(Time) %>%
  summarise(
    Mean = mean(mTEC, na.rm = TRUE),
    Median = median(mTEC, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -6.5, 6.5),
    Median_vjust = ifelse(Median > Mean, -6.5, 6.5)
  )

# 2. Plot
p4 <- ggplot(A5_MTEC, aes(x = Time, y = mTEC)) +
  geom_violin(aes(fill = Time)) +
  scale_fill_viridis_d() +
  theme_test() +
  stat_pvalue_manual(A5_MTEC_wilcox, label = "p.signif", hide.ns = TRUE) + 
  labs(
    title = "(B)",
    x = "Time",
    y = "Colony Counts in mTEC"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold.italic"),
        axis.text.y = element_text(vjust = 1, size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  
  # Add mean and median points
  geom_point(data = A5_MTEC_summary_stats, aes(x = Time, y = Mean), 
             color = "blue", size = 3, shape = 18) +
  geom_point(data = A5_MTEC_summary_stats, aes(x = Time, y = Median), 
             color = "red", size = 3, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = A5_MTEC_summary_stats, aes(x = Time, y = Mean, 
                                              label = paste0("Mean: ", round(Mean, 2)), 
                                              vjust = Mean_vjust), 
            color = "blue", size = 4) +
  geom_text(data = A5_MTEC_summary_stats, aes(x = Time, y = Median, 
                                              label = paste0("Median: ", round(Median, 2)), 
                                              vjust = Median_vjust), 
            color = "red", size = 4) 

p4    #### Plot 4

figure_2 <- p3 + p4 + plot_layout(widths = c(4, 1), guides = "collect")

# Save as high-resolution TIFF
ggsave("Figure_2.tiff", plot = figure_2, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")

figure_2 ##########################################################################################









########################################################--Figure3--################
##### P6: sixth Plot: Colony counts difference by HID
A5_MTEC_HID <- A5 %>%
  select(HID, EID, mTEC) %>% 
  filter(EID != "E3" & EID != "E4")

A5_friedman_test_result <- A5_MTEC_HID %>%
  friedman_test(mTEC ~ HID | EID)

print(A5_friedman_test_result)

#### post hoc, wilcoxon ranked signed
library(rstatix)
a5_MTEC_HID_wilcox <- A5_MTEC_HID %>%           
  pairwise_wilcox_test(
    mTEC ~ HID,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_xy_position(x = "HID") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE ~ "NS"
  ))

A5_MTEC_HID_effect <- A5_MTEC_HID %>%
  wilcox_effsize(mTEC ~ HID, paired = TRUE)  # For Wilcoxon signed-rank effect sizes

p6a <-ggplot(A5_MTEC_HID_effect, aes(x=paste(group1, "vs", group2), y=effsize, fill=magnitude)) +
  geom_col() +
  geom_hline(yintercept = c(-0.3, 0.3, -0.5, 0.5), linetype = "dashed") +
  labs(title = "(C)", y = "Rank-Biserial r", x = "Comparison")+
  theme_minimal()+
  scale_fill_brewer(palette = "RdYlBu")+
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90),
    axis.title = element_text(size = 14, face = "bold")
  ) 


p6a ### Plot 6a

#### Plot 6

A5_MTEC_HID_summary_stats <- A5_MTEC_HID%>%
  group_by(HID) %>%
  summarise(
    Mean = mean(mTEC, na.rm = TRUE),
    Median = median(mTEC, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -4.5, 4.5),
    Median_vjust = ifelse(Median > Mean, -4.5, 4.5)
  )
p6b <- ggplot(A5_MTEC_HID, aes(x = HID, y = mTEC)) +
  geom_violin(aes(fill = HID), trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_brewer(palette = "Set3") +
  stat_pvalue_manual(a5_MTEC_HID_wilcox, label = "p", hide.ns = TRUE,
                     size = 5,
                     bracket.size = 0.7,
                     step.increase = 0.1,
                     y_position = max(A5_MTEC_HID$mTEC) * 1.05
  ) +
  geom_point(data = A5_MTEC_HID_summary_stats, aes(y = Mean), 
             color = "blue", size = 3, shape = 18) +
  geom_point(data = A5_MTEC_HID_summary_stats, aes(y = Median), 
             color = "red", size = 3, shape = 17) +
  geom_text(data = A5_MTEC_HID_summary_stats, aes(y = Mean, 
                                                    label = paste0("Mean: ", round(Mean, 2)), 
                                                    vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = A5_MTEC_HID_summary_stats, aes(y = Median, 
                                                    label = paste0("Median: ", round(Median, 2)), 
                                                    vjust = Median_vjust), 
            color = "red", size = 3) +
  labs(
    title = "(A)",
    x = "HID",
    y = "Colony Counts in mTEC"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold")
  ) 

p6b   #### Plot 6a





##### P7: Seventh Plot: Colony counts difference by EID
A5_MTEC_EID <- A5 %>%
  select(HID, EID, mTEC) %>% 
  filter(EID != "E3" & EID != "E4")

A5_MTEC_EID_friedman_test <- A5_MTEC_EID %>%
  friedman_test(mTEC ~ EID | HID)
print(A5_MTEC_EID_friedman_test)

#### post hoc, wilcoxon ranked signed
library(rstatix)
a5_MTEC_EID_wilcox <- A5_MTEC_EID %>%           
  pairwise_wilcox_test(
    mTEC ~ EID,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_xy_position(x = "EID") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE ~ "NS"
  ))

### mwasuring effect size For Wilcoxon signed-rank test
A5_MTEC_EID_effect <- A5_MTEC_EID %>%
  wilcox_effsize(mTEC ~ EID, paired = TRUE)  

p7a <- ggplot(A5_MTEC_EID_effect, aes(x=paste(group1, "vs", group2), y=effsize, fill=magnitude)) +
  geom_col() +
  geom_hline(yintercept = c(-0.3, 0.3, -0.5, 0.5), linetype = "dashed") +
  labs(title = "(D)", y = "Rank-Biserial r", x = "Comparison")+
  theme_minimal()+
  scale_fill_brewer(palette = "RdYlBu")+
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank()
  ) 
  

p7a  #### plot 7A

A5_MTEC_EID_summary_stats <- A5_MTEC_EID%>%
  group_by(EID) %>%
  summarise(
    Mean = mean(mTEC, na.rm = TRUE),
    Median = median(mTEC, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -4.5, 4.5),
    Median_vjust = ifelse(Median > Mean, -4.5, 4.5)
  )


p7b <- ggplot(A5_MTEC_EID, aes(x = EID, y = mTEC)) +
  geom_violin(aes(fill = EID), trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel1") +
  stat_pvalue_manual(a5_MTEC_EID_wilcox, label = "p.adj.signif", hide.ns = TRUE,
                     size = 5,
                     bracket.size = 0.7,
                     step.increase = 0.1,
                     y_position = max(A5_MTEC_EID$mTEC) * 1.05
  ) +
  geom_point(data = A5_MTEC_EID_summary_stats, aes(y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = A5_MTEC_EID_summary_stats, aes(y = Median), 
             color = "red", size = 2, shape = 17) +
  geom_text(data = A5_MTEC_EID_summary_stats, aes(y = Mean, 
                                                  label = paste0("Mean: ", round(Mean, 2)), 
                                                  vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = A5_MTEC_EID_summary_stats, aes(y = Median, 
                                                  label = paste0("Median: ", round(Median, 2)), 
                                                  vjust = Median_vjust), 
            color = "red", size = 3) +
  labs(
    title = "(B)",
    x = "EID",
    y = "Colony Counts in MTEC"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text= element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank()
  ) 

p7b  ### Plot 7b

top_row <- p6b + p7b + plot_layout(widths = c(2, 1))
bottom_row <- p6a + p7a + plot_layout(widths = c(2, 1))

figure_3 <- top_row / bottom_row
figure_3

# Save as high-resolution TIFF
ggsave("Figure_1.tiff", plot = figure_3, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")

figure_3##########################################################################################








##### PX:   
a5_CT_diff_types_friedman_test <- A5_paired_filtered_new_columns %>%
  group_by(Types) %>% 
  friedman_test(CT_diff ~ Targets | HID)
print(a5_CT_diff_types_friedman_test)
######No need to do the following few lines as results were insignificant

a5_CT_diff_types_wilcox_test <- A5_paired_filtered_new_columns %>% 
  group_by(Types) %>% 
  wilcox_test(CT_diff~ Targets, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Targets")

# Calculate mean and median for each target
a5_CT_diff_types_summary_stats <- A5_paired_filtered_new_columns %>%
  group_by(Types, Targets) %>%
  summarise(
    Mean = mean(CT_diff, na.rm = TRUE),
    Median = median(CT_diff, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -1.5, 1.5),
    Median_vjust = ifelse(Median > Mean, -1.5, 1.5)
  )

# Plot
ggplot(A5_paired_filtered_new_columns, aes(x = Targets, y = CT_diff)) +
  geom_violin(aes(fill = Targets)) +
  scale_fill_brewer() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_test() +
  stat_pvalue_manual(a5_CT_diff_types_wilcox_test, label = "p.signif", hide.ns = TRUE) + 
  labs(title = "(a)",
       y = "Difference in CT values (Before - After)", x = "Targets") +
  theme(
    axis.text.x = element_text(
      face = "italic", hjust = 1, 
      margin = ggplot2::margin(b = 5), 
      size = 12
    ),
    plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 1.2),
    legend.position = "none"
  ) +
  
  # Add mean and median points
  geom_point(data = a5_CT_diff_types_summary_stats, aes(x = Targets, y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = a5_CT_diff_types_summary_stats, aes(x = Targets, y = Median), 
             color = "orange", size = 2, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = a5_CT_diff_types_summary_stats, aes(x = Targets, y = Mean, 
                                                       label = paste0("Mean: ", round(Mean, 2)), 
                                                       vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = a5_CT_diff_types_summary_stats, aes(x = Targets, y = Median, 
                                                       label = paste0("Median: ", round(Median, 2)), 
                                                       vjust = Median_vjust), 
            color = "orange", size = 3)+
  facet_wrap(~Types, nrow = 3)




##### PY:   
a5_CT_diff_HID_friedman_test <- A5_paired_filtered_new_columns %>%
  group_by(HID) %>% 
  friedman_test(CT_diff ~ Targets | Types)
print(a5_CT_diff_HID_friedman_test)

a5_CT_diff_HID_wilcox_test <- A5_paired_filtered_new_columns %>% 
  group_by(HID) %>% 
  wilcox_test(CT_diff~ Targets, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Targets")

# Calculate mean and median for each target
a5_CT_diff_HID_summary_stats <- A5_paired_filtered_new_columns %>%
  group_by(HID, Targets) %>%
  summarise(
    Mean = mean(CT_diff, na.rm = TRUE),
    Median = median(CT_diff, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -1.5, 1.5),
    Median_vjust = ifelse(Median > Mean, -1.5, 1.5)
  )

# Plot
ggplot(A5_paired_filtered_new_columns, aes(x = Targets, y = CT_diff)) +
  geom_violin(aes(fill = Targets)) +
  scale_fill_brewer() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_test() +
  stat_pvalue_manual(a5_CT_diff_HID_wilcox_test, label = "p.signif", hide.ns = TRUE) + 
  labs(title = "(a)",
       y = "Difference in CT values (Before - After)", x = "Targets") +
  theme(
    axis.text.x = element_text(
      face = "italic", hjust = 1, 
      margin = ggplot2::margin(b = 5), 
      size = 12
    ),
    plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 1.2),
    legend.position = "none"
  ) +
  
  # Add mean and median points
  geom_point(data = a5_CT_diff_HID_summary_stats, aes(x = Targets, y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = a5_CT_diff_HID_summary_stats, aes(x = Targets, y = Median), 
             color = "orange", size = 2, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = a5_CT_diff_HID_summary_stats, aes(x = Targets, y = Mean, 
                                                     label = paste0("Mean: ", round(Mean, 2)), 
                                                     vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = a5_CT_diff_HID_summary_stats, aes(x = Targets, y = Median, 
                                                     label = paste0("Median: ", round(Median, 2)), 
                                                     vjust = Median_vjust), 
            color = "orange", size = 3)+
  facet_wrap(~HID, nrow = 5)




##### P5: fifth Plot: Colony counts difference by types
A5_MTEC_types <- A5 %>%
  filter(Types != "Food") %>%  # Note the quotes around "Food"
  select(Types, Time, mTEC) 

a5_MTEC_wilcox_test <- A5_MTEC_types%>% 
  group_by(Types) %>% 
  wilcox_test(mTEC~ Time, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Time")

A5_MTEC_types_summary_stats <- A5_MTEC_types%>%
  group_by(Types, Time) %>%
  summarise(
    Mean = mean(mTEC, na.rm = TRUE),
    Median = median(mTEC, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -1.5, 1.5),
    Median_vjust = ifelse(Median > Mean, -1.5, 1.5)
  )

p5 <- ggplot(A5_MTEC_types, aes(x = Time, y = mTEC)) +
  geom_violin(aes(fill = Time), trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel1") +
  stat_pvalue_manual(a5_MTEC_wilcox_test, label = "p", hide.ns = TRUE,
                     size = 5,
                     bracket.size = 0.7,
                     step.increase = 0.1,
                     y_position = max(A5_MTEC_types$mTEC) * 1.05
  ) +
  geom_point(data = A5_MTEC_types_summary_stats, aes(y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = A5_MTEC_types_summary_stats, aes(y = Median), 
             color = "orange", size = 2, shape = 17) +
  geom_text(data = A5_MTEC_types_summary_stats, aes(y = Mean, 
                                                    label = paste0("Mean: ", round(Mean, 2)), 
                                                    vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = A5_MTEC_types_summary_stats, aes(y = Median, 
                                                    label = paste0("Median: ", round(Median, 2)), 
                                                    vjust = Median_vjust), 
            color = "orange", size = 3) +
  labs(
    title = "(b)",
    x = "Time",
    y = "Colony Counts in MTEC"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = -0.1, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13)
  ) +
  facet_wrap(~ Types)

p5  #### Plot 5




####figure_5##########################################################################################
##### P8: Eight Plot: Comparing CT value of each target across EID
A5_CT_EID <- A5_tidy %>%
  select(Targets, HID, EID, CT_value) %>% 
  filter(EID != "E3" & EID != "E4")

A5_CT_EID_friedman_test <- A5_CT_EID %>%
  group_by(Targets) %>% 
  friedman_test(CT_value ~ EID | HID)
print(A5_CT_EID_friedman_test)

#### post hoc, wilcoxon ranked signed
library(rstatix)
a5_CT_EID_wilcox <- A5_CT_EID %>%  
  group_by(Targets) %>% 
  pairwise_wilcox_test(
    CT_value ~ EID,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_xy_position(x = "EID") %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE ~ "NS"
  ))

### mwasuring effect size For Wilcoxon signed-rank test
A5_CT_EID_effect <- A5_CT_EID %>%
  group_by(Targets) %>% 
  wilcox_effsize(CT_value ~ EID, paired = TRUE) 

p8.1 <-ggplot(A5_CT_EID_effect, aes(x=paste(group1, "vs", group2), y=effsize, fill=magnitude)) +
  geom_col() +
  geom_hline(yintercept = c(-0.3, 0.3, -0.5, 0.5), linetype = "dashed") +
  labs(y = "Rank-Biserial r", x = "Comparison")+
  theme_calc()+
  scale_fill_brewer(palette = "RdYlBu")+
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 14, face = "bold.italic"),
    strip.placement = "outside", strip.clip = "off",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90),
    axis.title = element_text(size = 14, face = "bold")
  )  +
  facet_wrap(~Targets)

p8.1

ggsave("Figure_5.tiff", plot = p8.1, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
###figure_5##########################################################################################





###figure_6##########################################################################################
library(dplyr)
A5_CT_EID_summary_stats <- A5_CT_EID %>%
  group_by(Targets, EID) %>%
  summarise(
    Mean = mean(CT_value, na.rm = TRUE),
    Median = median(CT_value, na.rm = TRUE)
  ) %>%
  mutate(
    Mean_vjust = case_when(
      Mean > Median ~ -5.5,
      Mean < Median ~ 5.5,
      TRUE ~ -3  # When Mean == Median
    ),
    Median_vjust = case_when(
      Median > Mean ~ -5.5,
      Median < Mean ~ 5.5,
      TRUE ~ 5  # When Mean == Median
    )
  )

p8.2 <- ggplot(A5_CT_EID, aes(x = EID, y = CT_value)) +
  geom_violin(aes(fill = EID), trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel1") +
  stat_pvalue_manual(a5_CT_EID_wilcox, label = "p", hide.ns = TRUE,
                     size = 5,
                     bracket.size = 0.7,
                     step.increase = 0.1,
                     y_position = max(A5_CT_EID$CT_value) * 1.05
  ) +
  geom_point(data = A5_CT_EID_summary_stats, aes(y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = A5_CT_EID_summary_stats, aes(y = Median), 
             color = "red", size = 2, shape = 17) +
  geom_text(data = A5_CT_EID_summary_stats, aes(y = Mean, 
                                                  label = paste0("Mean: ", round(Mean, 2)), 
                                                  vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = A5_CT_EID_summary_stats, aes(y = Median, 
                                                  label = paste0("Median: ", round(Median, 2)), 
                                                  vjust = Median_vjust), 
            color = "red", size = 3) +
  labs(x = "EID", y = "CT value") +
  theme_test() +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 14, face = "bold.italic"),
    strip.placement = "outside", strip.clip = "off",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 90),
    axis.title = element_text(size = 14, face = "bold")
  )  + scale_y_continuous(limits = c(10, 50), oob = scales::oob_keep) +
  facet_wrap(~Targets)
p8.2

ggsave("Figure_6.tiff", plot = p8.2, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
###figure_6##########################################################################################

p8 <- p8.2 + p8.1 + plot_layout(widths = c(3, 1)) ### Plot 8
p8





  ### P9: Ninth Plot Comparing CT value of each target across HID
A5_CT_HID <- A5_tidy %>%
  select(Targets, EID, HID, CT_value) %>% 
  filter(EID != "E3" & EID != "E4")

A5_CT_HID_friedman_test <- A5_CT_HID %>%
  group_by(Targets) %>% 
  friedman_test(CT_value ~ HID | EID)
print(A5_CT_HID_friedman_test)

#### post hoc, wilcoxon ranked signed
library(rstatix)
a5_CT_HID_wilcox <- A5_CT_HID %>%  
  group_by(Targets) %>% 
  pairwise_wilcox_test(
    CT_value ~ HID,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_xy_position(x = "HID") %>%
  mutate(p.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE ~ "NS"
  ))

### mwasuring effect size For Wilcoxon signed-rank test
A5_CT_HID_effect <- A5_CT_HID %>%
  group_by(Targets) %>% 
  wilcox_effsize(CT_value ~ HID, paired = TRUE) 

ggplot(A5_CT_HID_effect, aes(x=paste(group1, "vs", group2), y=effsize, fill=magnitude)) +
  geom_col() +
  geom_hline(yintercept = c(-0.3, 0.3, -0.5, 0.5), linetype = "dashed") +
  labs(y = "Rank-Biserial r", x = "Comparison")+
  theme(
    axis.text.x = element_text(angle = 90)) +
  facet_wrap(~Targets)

A5_CT_HID_summary_stats <- A5_CT_HID%>%
  group_by(Targets, HID) %>%
  summarise(
    Mean = mean(CT_value, na.rm = TRUE),
    Median = median(CT_value, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -1.5, 1.5),
    Median_vjust = ifelse(Median > Mean, -1.5, 1.5)
  )

p9 <- ggplot(A5_CT_HID, aes(x = HID, y = CT_value)) +
  geom_violin(aes(fill = HID), trim = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_brewer(palette = "Pastel1") +
  stat_pvalue_manual(a5_CT_HID_wilcox, label = "p.adj.signif", hide.ns = TRUE,
                     size = 5,
                     bracket.size = 0.7,
                     step.increase = 0.1,
                     y_position = max(A5_CT_HID$CT_value) * 1.05
  ) +
  geom_point(data = A5_CT_HID_summary_stats, aes(y = Mean), 
             color = "blue", size = 2, shape = 18) +
  geom_point(data = A5_CT_HID_summary_stats, aes(y = Median), 
             color = "orange", size = 2, shape = 17) +
  geom_text(data = A5_CT_HID_summary_stats, aes(y = Mean, 
                                                label = paste0("Mean: ", round(Mean, 2)), 
                                                vjust = Mean_vjust), 
            color = "blue", size = 3) +
  geom_text(data = A5_CT_HID_summary_stats, aes(y = Median, 
                                                label = paste0("Median: ", round(Median, 2)), 
                                                vjust = Median_vjust), 
            color = "orange", size = 3) +
  labs(
    title = "(b)",
    x = "HID",
    y = "CT value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = -0.1, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13)
  ) + facet_wrap(~Targets)

p9   ### Plot 9





###################################################################################
##### p10: Heatmap representing variance in HID and EID for each target
# Rename p columns before joining to clarify
eid_p <- A5_CT_EID_friedman_test %>%
  select(Targets, p_eid = p)

hid_p <- A5_CT_HID_friedman_test %>%
  select(Targets, p_hid = p)

# Combine the two
friedman_pvalues <- left_join(eid_p, hid_p, by = "Targets")
print(friedman_pvalues)

library(gt)

friedman_pvalues %>%
  gt() %>%
  fmt_number(columns = c(p_eid, p_hid), decimals = 4) %>%
  tab_header(title = "Friedman Test p-values by Target Pathogen") %>%
  cols_label(p_eid = "p (across Events)", p_hid = "p (across Households)")

library(ggplot2)
library(tidyr)
library(dplyr)

# Pivot to long format
friedman_long <- friedman_pvalues %>%
  pivot_longer(cols = starts_with("p_"), names_to = "Test", values_to = "p_value") %>%
  mutate(
    Test = recode(Test,
                  p_eid = "Across Events",
                  p_hid = "Across Households"),
    Sig = ifelse(p_value < 0.05, "*", "")
  )

# Custom color scale: dark red for low p-values
p10 <- ggplot(friedman_long, aes(x = Test, y = Targets, fill = p_value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(p_value, 3)), color = "black", size = 4) +
  scale_fill_gradientn(
    colors = c("darkred", "yellow", "white"),
    name = expression(italic(p)*"-value"),
    limits = c(0, 1)
  )  +
  theme_test() +
  labs(x = NULL, y = "Target Genes") +
  theme(axis.text.y = element_text(size = 14, face = "bold.italic"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"))+
  theme_wsj()

p10

# Save as high-resolution TIFF
ggsave("Figure_4.tiff", plot = p10, width = 8, height = 6, dpi = 300, units = "in", device = "tiff")
##########################################--figure4--####################################







##### P6: sixth Plot: CT values with events

library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)


# Run Friedman test per target
# Kruskal-Wallis test
A5_kw_results <- A5_tidy %>%
  group_by(Targets) %>%
  kruskal_test(CT_value ~ EID) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Pairwise post-hoc: Dunn test
A5_posthoc_dunn <- A5_tidy %>%
  group_by(Targets) %>%
  dunn_test(CT_value ~ EID, p.adjust.method = "bonferroni") %>%
  add_significance("p.adj") %>%
  mutate(y.position = 39)


ggplot(A5_tidy, aes(x = EID, y = CT_value, color = HID)) +
  geom_boxplot(aes(group = HID), color = "black", alpha = 0.3, width = 0.5, outlier.shape = NA) +  # All HIDs combined boxplot
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  stat_pvalue_manual(A5_posthoc_dunn, label = "p.adj.signif", tip.length = 0.01, size = 3, hide.ns = TRUE) +
  facet_wrap(~ Targets, nrow = 3, ncol = 3) +
  labs(
    title = "CT Value Comparison Across EIDs (Kruskal-Wallis + Dunn)",
    x = "EID",
    y = "CT Value",
    color = "HID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines"),
    legend.position = "bottom"
  ) +
  scale_y_continuous(limits = c(10, 40)) 


ggplot(A5_tidy, aes(x = EID, y = CT_value, color = HID)) +
  geom_boxplot(aes(group = HID), color = "black", alpha = 0.3, width = 0.5, outlier.shape = NA) +  # All HIDs combined boxplot
  geom_point(size = 2, alpha = 0.8, position = position_jitter(width = 0.2)) +  # Individual HID points
  facet_wrap(~ Targets, nrow = 3, ncol = 3) +
  labs(
    title = "CT Values by Sample and Pathogen Target",
    x = "Sample ID (EID)",
    y = "CT Value",
    color = "Household (HID)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_y_continuous(limits = c(10, 40))





library(lme4)
library(emmeans)
library(broom.mixed)

# Run model per target
A5_lme_results <- A5_tidy %>%
  group_by(Targets) %>%
  group_map(~ {
    model <- lmer(CT_value ~ EID + (1 | HID), data = .x)
    emmeans(model, pairwise ~ EID, adjust = "bonferroni")$contrasts %>%
      as.data.frame() %>%
      mutate(Targets = unique(.x$Targets))
  }) %>%
  bind_rows() %>%
  filter(p.value < 0.05) %>%
  mutate(
    y.position = 39,
    group1 = as.character(contrast) %>% stringr::str_extract("^[^ ]+"),
    group2 = as.character(contrast) %>% stringr::str_extract("[^ ]+$"),
    p.adj.signif = symnum(p.value, corr = FALSE, na = FALSE,
                          cutpoints = c(0, .001, .01, .05, .1, 1),
                          symbols = c("***", "**", "*", ".", "ns"))
  )
ggplot(A5_tidy, aes(x = EID, y = CT_value, color = HID)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  stat_pvalue_manual(A5_lme_results, label = "p.adj.signif", tip.length = 0.01, size = 3) +
  facet_wrap(~ Targets, nrow = 3, ncol = 3) +
  labs(
    title = "CT Value Comparison Across EIDs (Mixed Effects Model)",
    x = "EID",
    y = "CT Value",
    color = "HID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines"),
    legend.position = "bottom"
  ) +
  scale_y_continuous(limits = c(10, 40)) +
  scale_color_brewer(palette = "Set1")









##### Trying ML


###########################Based on untrasnformed mTEC############################
# Simple model with one predictor
simple_model <- glm(mTEC ~ `E. coli (uidA)`, 
                    family = gaussian(link = "identity"), 
                    data = A5)

# Full model with all candidate predictors
full_model <- glm(mTEC ~ `V. cholerae (ctxA)` + `Shigella (ipaH)` + 
                    `V. cholerae (toxR)` + `E. coli (uidA)` + 
                    `Salmonella (invA)` + Adenovirus + Rotavirus + 
                    Norovirus + `C. jejuni (cadF)`,
                  family = gaussian(link = "identity"), 
                  data = A5)

# Forward selection from simple model toward full model
sw_forward <- step(simple_model,
                   scope = list(lower = formula(simple_model), upper = formula(full_model)),
                   direction = "forward")
summary(sw_forward) 
plot(sw_forward)

# Backward elimination from full model
sw_backward <- step(full_model, direction = "backward")
summary(sw_backward) 
plot(sw_backward)

# Feature section based on boruta algorithm
set.seed(12345)
A5$rand_vals <- runif(n = 58, min = 1, max = 100)
library(Boruta)
A5_boruta <- Boruta(mTEC ~ `V. cholerae (ctxA)` + `Shigella (ipaH)` + 
                      `V. cholerae (toxR)` + `E. coli (uidA)` + 
                      `Salmonella (invA)` + Adenovirus + Rotavirus + 
                      Norovirus + `C. jejuni (cadF)` + rand_vals,
                           data = A5, doTrace = 1, maxRuns =1000)
summary(A5_boruta)
plot(A5_boruta)

## Checking correlation 
cor.test(A5$mTEC, A5$Norovirus, method = "spearman")
cor.test(A5$mTEC, A5$`V. cholerae (toxR)`, method = "spearman")

## Building randome forest model
library(randomForest)
rf_model <- randomForest(mTEC ~ `V. cholerae (ctxA)` + `Shigella (ipaH)` + 
                           `V. cholerae (toxR)` + `E. coli (uidA)` + 
                           `Salmonella (invA)` + Adenovirus + Rotavirus + 
                           Norovirus + `C. jejuni (cadF)` + rand_vals,
                         data = A5, importance = TRUE)

## checking quality of the model
varImpPlot(rf_model)





###########################Based on untrasnformed mTEC############################

A$log_mTEC <- log10(A5$mTEC + 1)  # +1 to avoid log(0)

simple_model <- glm(log_mTEC ~ `E. coli (uidA)`, data = A5)

full_model <- glm(log_mTEC ~ `V. cholerae (ctxA)` + `Shigella (ipaH)` + 
                    `V. cholerae (toxR)` + `E. coli (uidA)` + 
                    `Salmonella (invA)` + Adenovirus + Rotavirus + 
                    Norovirus + `C. jejuni (cadF)`, data = A5)

# Forward selection from simple model toward full model
sw_forward <- step(simple_model,
                   scope = list(lower = formula(simple_model), upper = formula(full_model)),
                   direction = "forward")
summary(sw_forward) 
plot(sw_forward)

# Backward elimination from full model
sw_backward <- step(full_model, direction = "backward")
summary(sw_backward) 
plot(sw_backward)

# Boruta
set.seed(12345)
A5$rand_vals <- runif(nrow(A5), 1, 100)
A5_boruta <- Boruta(log_mTEC ~ `V. cholerae (ctxA)` + `Shigella (ipaH)` + 
                      `V. cholerae (toxR)` + `E. coli (uidA)` + 
                      `Salmonella (invA)` + Adenovirus + Rotavirus + 
                      Norovirus + `C. jejuni (cadF)` + rand_vals,
                    data = A5, doTrace = 1, maxRuns =1000)
summary(A5_boruta)
plot(A5_boruta)

# Correlation
cor.test(A5$log_mTEC, A5$`V. cholerae (ctxA)`, method = "spearman")
####################################################################################



###multicolineartiy analysis######################################################
car::vif(full_model)
cor(A5[c("log_mTEC", "mTEC", "V. cholerae (ctxA)", "Shigella (ipaH)", 
         "V. cholerae (toxR)", "E. coli (uidA)",
         "Salmonella (invA)", "Adenovirus", "Rotavirus",
         "Norovirus", "C. jejuni (cadF)")], use = "pairwise.complete.obs")


pairs(A5[c("log_mTEC", "mTEC", "V. cholerae (ctxA)", "Shigella (ipaH)", 
           "V. cholerae (toxR)", "E. coli (uidA)",
           "Salmonella (invA)", "Adenovirus", "Rotavirus",
           "Norovirus", "C. jejuni (cadF)")])

library(psych)

pairs.panels(A5[c("log_mTEC", "mTEC", "V. cholerae (ctxA)", "Shigella (ipaH)", 
                  "V. cholerae (toxR)", "E. coli (uidA)",
                  "Salmonella (invA)", "Adenovirus", "Rotavirus",
                  "Norovirus", "C. jejuni (cadF)")], pch = ".")

library(corrplot)
library(Hmisc)

# Step 1: Select numeric columns from A5
cor_data <- A5[c("log_mTEC", "mTEC", "V. cholerae (ctxA)", "Shigella (ipaH)", 
                 "V. cholerae (toxR)", "E. coli (uidA)",
                 "Salmonella (invA)", "Adenovirus", "Rotavirus",
                 "Norovirus", "C. jejuni (cadF)")]

# Step 2: Compute correlation matrix, handle missing values
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

# Step 3: Visualize correlation matrix
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45)

# Calculate correlation matrix and p-values
cor_results <- rcorr(cor_matrix)

# extract correlation matrix and p-values
cor_matrix <- cor_results$r
p_values <- cor_results$p

#mask insignificant correlations 
significant_cor <- cor_matrix
significant_cor[p_values > 0.05] <- NA

pretty_labels <- c("log_mTEC", 
                   "mTEC", 
                   "V. cholerae\n(ctxA)", 
                   "Shigella\n(ipaH)", 
                   "V. cholerae\n(toxR)", 
                   "E. coli\n(uidA)", 
                   "Salmonella\n(invA)", 
                   "Adenovirus", 
                   "Rotavirus", 
                   "Norovirus", 
                   "C. jejuni\n(cadF)")
colnames(significant_cor) <- pretty_labels
rownames(significant_cor) <- pretty_labels


tiff("Figure_7.tiff", width = 12, height = 9, units = "in", res = 300)
corrplot(significant_cor, method = 'ellipse', # Use ellipses to show strength/direction
         na.label = "",                       # Hide NA cells
         addCoef.col = "black",               # Show correlation values inside cells
         number.cex = 0.7,                    # Size of text for coefficients
         tl.col = "black",                    # Text color for labels
         tl.pos = "d",                        # Labels on the diagonal
         diag = FALSE)                        # Hide diagonal                       
dev.off()

###figure_6##########################################################################################
## The oval-shaped object on each scatterplot is a correlation ellipse.
## It provides a simple visual indicator of correlation strength. In this dataset, there are no strong
## correlations, so the ovals are mostly flat; with stronger correlations, the ovals would be tilted
## upward or downward to indicate a positive or negative correlation. The dot at the center of the
## ellipse is a point reflecting the means of the x- and y-axis variables.
## performance::check_model(full_model)  # Diagnose assumptions

## The line superimposed across the scatterplot (blue) is called a loess curve. 
## It indicates the general relationship between the x-axis and y-axis variables.
## the loess curve can sometimes be quite dramatic with V- or U-shaped curves as well
## as stairstep patterns. Recognizing such patterns can assist later with developing
## a better-fitting regression model.





# Regression model based on log_MTEC
r1_model <- lm(log_mTEC ~ toxR + norovirus + ctx + uidA + invA +
                adeno + rotavirus + cadF + ipah,
              data = A5)
## options(scipen = 999) command turns off scientific notation to make the output easier to read:
options(scipen = 999)
r1_model
## Evaluating the model performance
summary(r1_model)
## residual is equal to the true value minus the predicted value

## the p-value, denoted by Pr(>|t|), provides an
## estimate of the probability that the true coefficient is zero given the value of the estimate.
## Small p-values suggest that the true coefficient is very unlikely to be zero, which means
## that the feature is extremely unlikely to have no relationship with the dependent variable.

## Multiple R-squared value (also called the coefficient of determination) provides a
## measure of how well our model as a whole explains the values of the dependent variable.
## It is similar to the correlation coefficient in that the closer the value is to 1.0, the better
## the model perfectly explains the data. Since the R-squared value is 0.01241, we know that
## the model explains about 1.2 percent of the variation in the dependent variable. Because
## models with more features always explain more variation, the Adjusted R-squared value
## corrects R-squared by penalizing models with a large number of independent variables.

## The direction of the relationship between the predictors and the target outcome 
## can also be understood simply by looking at the sign (positive or negative) 
## before the Estimate value.

# Regression model based on MTEC
r2_model <- lm(MTEC ~ toxR + norovirus + ctx + uidA + invA +
                adeno + rotavirus + cadF + ipah,
              data = A5)
## options(scipen = 999) command turns off scientific notation to make the output easier to read:
options(scipen = 999)
r2_model
## Evaluating the model performance
summary(r2_model)




# Random Forest
library(randomForest)
rf_model <- randomForest(log_MTEC ~ toxR + norovirus + ctx + uidA + invA +
                           adeno + rotavirus + cadF + ipah,
                         data = A5, importance = TRUE)
varImpPlot(rf_model)


A5 %>%
  group_by(HID) %>%
  pairwise_wilcox_test(log_MTEC ~ EID, p.adjust.method = "bonferroni")

