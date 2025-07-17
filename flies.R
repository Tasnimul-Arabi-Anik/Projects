### let's set our working directory, where my datasets are kept
setwd("/home/anik/genomics/Projects/Reham/")

### You can check the current directory
getwd()




###figure_1##########################################################################################
C1 <- read.csv("file3.csv", check.names = FALSE)

C1_tidy <- C1 %>% 
  select(House_ID, Pathotypes, Exposed, `Non-exposed`) %>%
  mutate(across(c(Exposed, `Non-exposed`), ~ as.numeric(.))) %>%
  pivot_longer(
    cols = c(Exposed, `Non-exposed`),
    names_to = "Groups",
    values_to = "CFU"
  ) %>%
  mutate(
    CFU = ifelse(is.na(CFU), 0, CFU)
  )

#C1 wilcox test
C1_wilcox_test <- C1_tidy %>% 
  wilcox_test(CFU ~ Groups, paired = TRUE) %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "NS"
  )) %>%
  add_xy_position(x = "Groups")

# Compute mean, median, and non-zero count
C1_summary_stats <- C1_tidy %>%
  group_by(Groups) %>%
  summarise(
    Mean = mean(CFU),
    Median = median(CFU),
    NonZeroCount = sum(CFU > 0)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -2.5, 2.5),
    Median_vjust = ifelse(Median > Mean, -2.5, -1.0)
  )

# Plot
p1 <-ggplot(C1_tidy, aes(x = Groups, y = CFU)) +
  geom_violin(aes(fill = Groups)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1, height = 0.1)) +
  scale_fill_brewer(palette = 'Set3') +
  theme_test() +
  stat_pvalue_manual(C1_wilcox_test, label = "p.signif", hide.ns = TRUE) + 
  labs(y = expression("Thermotolerant " ~ italic(E.~coli) ~ " in mTEC (CFU/g)"), x = "Groups") +
  theme(axis.text.x = element_text(face = "bold", hjust = 1, size = 12),
        axis.text.y = element_text(face = "bold", size = 12), 
        axis.title = element_text(size = 15, face = "bold"),
        panel.border = element_rect(linewidth = 1.2),
        legend.position = "none") +
  # Add Mean and Median lines
  geom_hline(data = C1_summary_stats, aes(yintercept = Mean), color = "blue", linetype = "dashed") +
  geom_hline(data = C1_summary_stats, aes(yintercept = Median), color = "red", linetype = "dotted") +
  # Annotate Mean, Median, and Non-zero Count
  # Add mean and median points
  geom_point(data = C1_summary_stats, aes(x = Groups, y = Mean), 
             color = "blue", size = 4, shape = 18) +
  geom_point(data = C1_summary_stats, aes(x = Groups, y = Median), 
             color = "red", size = 4, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = C1_summary_stats, aes(x = Groups, y = Mean, 
                                         label = paste0("Mean: ", round(Mean, 2)), 
                                         vjust = Mean_vjust, hjust = 1.4), 
            color = "blue", size = 7) +
  geom_text(data = C1_summary_stats, aes(x = Groups, y = Median, 
                                         label = paste0("Median: ", round(Median, 2)), 
                                         vjust = Median_vjust, hjust = 1.4), 
            color = "red", size = 7) +
  geom_text(data = C1_summary_stats, aes(x = Groups, y = max(C1_tidy$CFU, na.rm = TRUE) +2, 
                                         label = paste0("n>0: ", NonZeroCount)), vjust = -0.4,
            size = 7, color = "black")

p1

ggsave("Figure_1.tiff", plot = p1, width = 10, height = 8, dpi = 300, units = "in", device = "tiff")
###figure_1##########################################################################################





###figure_2##########################################################################################
#C1 summary

C1_summary <- C1 %>%
  count(Pathotypes, name = "Count") %>%
  mutate(Percentage = Count / sum(Count) * 100)

p2b <-ggplot(C1_summary, aes(x = reorder(Pathotypes, -Count), y = Count, fill = Pathotypes)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = paste0(Count, " (", round(Percentage, 1), "%)")), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_brewer() +
  theme_test(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold", hjust = 1, size = 12),
        axis.text.y = element_text(face = "bold", size = 12), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.position = "none") +
  labs(
    x = "Pathotypes",
    y = expression("Detection (%) of " ~ italic(E.~coli) ~ " pathotypes in exposed group"),
    title = "(B)",
  )
p2b
###figure_2b##########################################################################################





###figure_2b##########################################################################################

### Uploading the dataset 
A1 <- read.csv("file1.csv", check.names = FALSE)

### Cleaning the dataset (undetermined = undetected, any value above 35 transformed to undetected)
A2 <- A1
A2[A2 == "Undetermined" | A2 == "undetermined"] <- "undetected"

####### Comparative Analysis ###################################################

## Converting to numeric data
A1_quant <- A2 %>% 
  mutate(across(2:11, ~ {
    ifelse(. == "undetected", 40, as.numeric(.))  # Convert "undetected" to 40, others to numeric
  })) %>% 
  select(-exposed)

## Converting to tidy data
A1_tidy <- A1_quant %>%
  pivot_longer(
    cols = 2:11,
    names_to = "Targets",
    values_to = "CT_value")

# Separate Targets and Controls
A1_sep <- A1_tidy %>%
  mutate(Target_Type = ifelse(grepl("_C$", Targets), "Non-exposed", "Exposed"),
         Targets = gsub("_C$", "", Targets))

# Converting to categorical
A1_cat <- A1_sep %>% 
  mutate(Real_time_detection = ifelse(CT_value < 40, "Detected", "Not Detected")) %>% 
  group_by(Targets, Target_Type) %>% 
  count(Real_time_detection)

p2a <- ggplot(A1_cat, aes(x = Targets, y = n, fill = Real_time_detection)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.7),
    vjust = +0.1,
    size = 4,
    fontface = "bold"
  ) +
  facet_wrap(~Target_Type) +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "(A)", x = "Targets", y = "Frequency", fill = "Detection Status") + #
  theme_test() +
  theme(axis.text.x = element_text(size = 12, face = "bold.italic"),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        # This adds a black border to legend keys
        legend.key = element_rect(color = "black", size = 0.5))
p2a
p2 <- p2a + p2b + plot_layout(widths = c(2, 1)) 
p2

ggsave("Figure_2.tiff", plot = p2, width = 16, height = 8, dpi = 300, units = "in", device = "tiff")
###figure_2##########################################################################################


# Pivot Wider to create paired columns
A1_paired <- A1_sep %>%
  pivot_wider(names_from = Target_Type, values_from = CT_value) %>%
  filter(!is.na(Exposed) & !is.na(`Non-exposed`))  # Keep only complete pairs

## Create a diff colum
A1_paired_diff <- A1_paired %>%
  mutate(Diff = `Non-exposed` - `Exposed`)

library(broom)

## Normality check
A1_normality_check <- A1_paired_diff %>%
  group_by(Targets) %>%
  group_modify(~ tidy(shapiro.test(.x$Diff)))

## Stat test and visualization
# Calculate mean and median for each target
A1_summary_stats <- A1_paired_diff %>%
  group_by(Targets) %>%
  summarise(
    Mean = mean(Diff, na.rm = TRUE),
    Median = median(Diff, na.rm = TRUE)
  ) %>%
  mutate(
    # Dynamically adjust position to prevent overlap
    Mean_vjust = ifelse(Mean > Median, -1.5, 1.5),
    Median_vjust = ifelse(Median > Mean, -1.5, 1.5)
  )

# Compute differences and annotate p-values and n
A1_t_test_summary <- A1_paired_diff %>%
  group_by(Targets) %>%
  summarise(
    p_value = t.test(Exposed, `Non-exposed`, paired = TRUE)$p.value,
    n = n(),
    .groups = "drop"
  )

library(rstatix)

A1_wilcox_test <- A1_paired_diff %>% 
  pairwise_wilcox_test(
    Diff ~ Targets,
    paired = TRUE,
    p.adjust.method = "BH"
  ) %>%
  add_xy_position(x = "Targets") %>%
  mutate(p.signif = case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE ~ "NS"
  ))


# Plot
p3a <- ggplot(A1_paired_diff, aes(x = Targets, y = Diff)) +
  geom_violin(aes(fill = Targets)) +
  scale_fill_brewer()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_test() +
  stat_pvalue_manual(A1_wilcox_test, label = "p.adj.signif", hide.ns = TRUE) + 
  labs(title = "(A)",
       y = "Difference in CT values (Non-exposed - exposed)", x = "Targets") +
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
  geom_point(data = A1_summary_stats, aes(x = Targets, y = Mean), 
             color = "blue", size = 4, shape = 18) +
  geom_point(data = A1_summary_stats, aes(x = Targets, y = Median), 
             color = "red", size = 4, shape = 17) +
  
  # Add text labels for mean and median with dynamic adjustment
  geom_text(data = A1_summary_stats, aes(x = Targets, y = Mean, 
                                         label = paste0("Mean: ", round(Mean, 2)), 
                                         vjust = Mean_vjust), 
            color = "blue", size = 4) +
  geom_text(data = A1_summary_stats, aes(x = Targets, y = Median, 
                                         label = paste0("Median: ", round(Median, 2)), 
                                         vjust = Median_vjust), 
            color = "red", size = 4) +
  
  # Annotate p-value and number of samples
  geom_text(
    data = A1_t_test_summary,
    aes(x = Targets, y = max(A1_paired_diff$Diff, na.rm = TRUE) + 13,
        label = paste0("p = ", signif(p_value, 3), "\n", "n = ", n)),
    size = 5
  )
p3a

### mwasuring effect size For Wilcoxon signed-rank test
A1_effect <- A1_paired_diff %>%
  wilcox_effsize(Diff ~ Targets, paired = TRUE)  

p3b <-ggplot(A1_effect, aes(x=paste(group1, "vs", group2), y=effsize, fill=magnitude)) +
  geom_col() +
  geom_hline(yintercept = c(-0.3, 0.3, -0.5, 0.5), linetype = "dashed") +
  labs(title = "(B)", y = "Rank-Biserial r", x = "Comparison")+
  theme_minimal()+
  scale_fill_brewer()+
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold.italic", angle = 90),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    # This adds a black border to legend keys
    legend.key = element_rect(color = "black", size = 0.5))

p3b


p3 <- p3a + p3b + plot_layout(widths = c(2, 1)) 
p3

ggsave("Figure_3.tiff", plot = p3, width = 14, height = 8, dpi = 300, units = "in", device = "tiff")
###figure_3##########################################################################################






##### Trying ML
A3 <- A2 %>% 
  mutate(across(2:11, ~ {
    ifelse(. == "undetected", 40, as.numeric(.))  # Convert "undetected" to 40, others to numeric
  })) %>% 
  select(House_ID, ipaH, invA, uidA, toxR,ctxA, exposed)

write.csv(A3, file = "A3.csv")


###########################Based on untrasnformed mTEC############################
# Simple model with one predictor
simple_model <- glm(exposed ~ uidA, 
                    family = gaussian(link = "identity"), 
                    data = A3)
# Full model with all candidate predictors
full_model <- glm(exposed ~ ipaH+ invA+  uidA+ toxR+ ctxA,
                  family = gaussian(link = "identity"), 
                  data = A3)
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
A3$rand_vals <- runif(n = 100, min = 20, max = 40)
library(Boruta)
A3_boruta <- Boruta(exposed ~ ipaH+ invA+  uidA+ toxR +ctxA+ rand_vals,
                    data = A3, doTrace = 1, maxRuns =1000)
summary(A5_boruta)
plot(A5_boruta)


###########################Based on untrasnformed mTEC############################

A3$log_mTEC <- log10(A3$exposed + 1)  # +1 to avoid log(0)

# Simple model with one predictor
simple_model <- glm(log_mTEC ~ uidA,
                    data = A3)
# Full model with all candidate predictors
full_model <- glm(log_mTEC  ~ ipaH+ invA+  uidA+ toxR+ctxA,
                  data = A3)
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
A3$rand_vals <- runif(n = 100, min = 20, max = 40)
library(Boruta)
A3_boruta <- Boruta(log_mTEC  ~ ipaH+ invA+  uidA+ toxR + ctxA+rand_vals,
                    data = A3, doTrace = 1, maxRuns =1000)
summary(A5_boruta)
plot(A5_boruta)

###multicolineartiy analysis######################################################
car::vif(full_model)

cor(A3[c("ipaH", "invA", "uidA", "ctxA", "exposed", "log_mTEC")], use = "pairwise.complete.obs")


pairs(A3[c("ipaH", "invA", "uidA", "ctxA", "toxR", "exposed", "log_mTEC")])

library(psych)

pairs.panels(A3[c("ipaH", "invA", "ctxA","uidA", "toxR", "exposed", "log_mTEC")], pch = ".")

library(corrplot)
library(Hmisc)

# Step 1: Select numeric columns from A5
cor_data <- A3[c("ipaH", "invA", "uidA", "toxR", "ctxA")]

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

     # pretty_labels <- c("log_mTEC", "mTEC",  "V. cholerae\n(ctxA)", "Shigella\n(ipaH)", "V. cholerae\n(toxR)", 
                         "E. coli\n(uidA)", 
                         "Salmonella\n(invA)", 
                         "Adenovirus", 
                         "Rotavirus", 
                         "Norovirus", 
                         "C. jejuni\n(cadF)")
    #  colnames(significant_cor) <- pretty_labels
    #  rownames(significant_cor) <- pretty_labels


tiff("Figure_4.tiff", width = 12, height = 9, units = "in", res = 300)
corrplot(significant_cor, method = 'ellipse', # Use ellipses to show strength/direction
         na.label = "",                       # Hide NA cells
         addCoef.col = "black",               # Show correlation values inside cells
         number.cex = 0.9,                    # Size of text for coefficients
         tl.col = "black",                    # Text color for labels
         tl.pos = "d",                        # Labels on the diagonal
         tl.srt  = 90,
         diag = FALSE)                       # Hide diagonal
                                           
dev.off()









B1 <- read.csv("file2.csv", check.names = FALSE)
B1_chi_data <- B1 %>%
  group_by(Targets, Conventional_PCR) %>%
  summarise(
    Undetermined = sum(RT_Undetermined, na.rm = TRUE),
    Determined = sum(RT_Determined, na.rm = TRUE),
    .groups = "drop"
  )
library(dplyr)

# Remove rows where both counts are zero
B1_chi_data_clean <- B1_chi_data %>%
  filter(Undetermined > 0 | Determined > 0)

# Perform Chi-square or Fisher's test
B1_chi_results <- B1_chi_data_clean %>%
  group_by(Targets) %>%
  summarise(
    # Only proceed if there are exactly 2 rows (e.g., Positive & Negative)
    test_type = if (n() == 2) {
      # Check if all expected values are >= 5
      if (all(chisq.test(cbind(Undetermined, Determined))$expected >= 5)) {
        "Chi-square"
      } else {
        "Fisher"
      }
    } else {
      NA
    },
    test_result = list(
      if (!is.na(test_type) && test_type == "Chi-square") {
        chisq.test(cbind(Undetermined, Determined))
      } else if (!is.na(test_type) && test_type == "Fisher") {
        fisher.test(cbind(Undetermined, Determined))
      } else {
        NULL
      }
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(test_result, function(x) if (!is.null(x)) x$p.value else NA),
    statistic = sapply(test_result, function(x) if (!is.null(x)) x$statistic else NA),
    testable = ifelse(is.na(p_value), "Not Testable", "Tested")
  )

# Display results
print(B1_chi_results)

library(ggplot2)

# Plot visualization
ggplot(B1_chi_results %>% filter(!is.na(p_value)), aes(x = Targets, y = -log10(p_value), fill = test_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = signif(p_value, 3)), vjust = -0.5, size = 3) +
  labs(y = "-log10(p-value)", x = "Targets", fill = "Test Type") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = -log10(0.05) + 0.3, label = "Significance Threshold (0.05)", size = 3, color = "red")







