setwd("/home/anik/genomics/Projects/fahim")
a1 <-read.csv('fahimA.csv', check.names = FALSE) 

library(dplyr)
library(tidyr)

a2 <- a1 %>%
  pivot_longer(
    cols = 3:11,
    names_to = "Antibiotics",
    values_to = "Zone_Diameter"
  ) %>% 
  mutate(Zone_Diameter = as.numeric(Zone_Diameter))

paired_data <- a2 %>%
  pivot_wider(names_from = Mestal, values_from = Zone_Diameter) %>%
  filter(!is.na(Presence) & !is.na(Absence)) 

#Shapiro test
paired_data %>%
  group_by(Antibiotics) %>%
  summarise(
    p_value = shapiro.test(Presence - Absence)$p.value
  )

# Compute differences and annotate p-values and n
t_test_summary <- paired_data %>%
  group_by(Antibiotics) %>%
  summarise(
    p_value = t.test(Presence, Absence, paired = TRUE)$p.value,
    n = n(),
    .groups = "drop"
  )

# Add a column of differences for plotting
paired_data <- paired_data %>%
  mutate(Diff = Presence - Absence)


# Create the boxplot with annotations
ggplot(paired_data, aes(x = Antibiotics, y = Diff)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  ylab("Difference in Zone Diameter (Presence - Absence)") +
  xlab("Antibiotics") +
  geom_text(
    data = t_test_summary,
    aes(x = Antibiotics, y = max(paired_data$Diff, na.rm = TRUE) + 1,
        label = paste0("p = ", signif(p_value, 3), "\n", "n = ", n)),
    size = 3
  )
