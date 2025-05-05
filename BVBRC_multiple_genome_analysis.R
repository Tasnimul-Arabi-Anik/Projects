setwd("/home/anik/genomics/Projects/delftia/patric")
A1 <- read.delim("Delftia_output_pathways.tsv", header = TRUE, sep = "\t")
A1$genome_name <- gsub("^Delftia tsuruhatensis ", "", A1$genome_name)
B1 <- read.delim("Delftia_output_subsystems.tsv", header = TRUE, sep = "\t")
B1$genome_name <- gsub("^Delftia tsuruhatensis ", "", B1$genome_name)

library(dplyr)
library(ggplot2)


#######Plot1
A2 <- A1 %>%
  group_by(genome_name, pathway_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(genome_name = fct_reorder(genome_name, n, sum))  # Order by total count

# Create plot
pub_plot <- ggplot(A2, aes(x = genome_name, y = n, fill = pathway_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +  # Use geom_col for pre-aggregated data
  facet_wrap(
    ~ stringr::str_wrap(pathway_class, width = 30),  # Wrap text at 15 characters
    scales = "free_x", 
    nrow = 3
  ) +
  
  # Improved aesthetics
  scale_fill_viridis_d(option = "plasma", guide = "none") +  # Remove redundant legend
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Pad top only
  
  # Professional theme adjustments
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 11),
    strip.placement = "outside", strip.clip = "off",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(
    x = "Genomes",
    y = "Gene Count"
  )

# Save with publication quality
ggsave(
  filename = "figure2.tiff",
  plot = pub_plot,
  device = "tiff",
  width = 16,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)







#######Plot2
A3 <-A1%>%
  group_by(genome_name, pathway_class, pathway_name) %>%
  count()
ggplot(A3, aes(x = pathway_name, y = n, fill = pathway_class)) +  
  geom_violin(alpha = 0.6, trim = FALSE) + 
  facet_wrap(
    ~ stringr::str_wrap(pathway_class, width = 30),  # Wrap text at 15 characters
    scales = "free_x", 
    nrow = 2
  )+
# Improved aesthetics
  scale_fill_viridis_d(option = "plasma", guide = "none") +  # Remove redundant legend
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Pad top only
  
  # Professional theme adjustments
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 10),
    strip.placement = "outside", strip.clip = "off",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")
  ) +
  labs(
    x = "Pathway Class",
    y = "Gene Count",
  )

# Save with publication quality
ggsave(
  filename = "figure3.tiff",
  device = "tiff",
  width = 16,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)








#######Plot3
B2 <-B1%>%
  group_by(genome_name, superclass, class) %>%
  count()%>%
  mutate(superclass = tolower(superclass))
# 1. Complete missing combinations with zeros (only genome_name × class)
B2_complete <- B2 %>%
  ungroup() %>%
  tidyr::complete(genome_name, superclass, class, fill = list(n = 0))

superclass_info <- data.frame(
  superclass = c("cell envelope", "cellular processes", "dna processing", "energy",
                 "membrane transport", "metabolism", "miscellaneous", "protein processing",
                 "regulation and cell signaling", "rna processing", 
                 "stress response, defense, virulence", ""),
  superclass_label = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "Z"),
  description = c(
    "A: Cell wall/membrane biogenesis",
    "B: Cell division, motility, etc.",
    "C: DNA replication/repair",
    "D: Energy production/conversion",
    "E: Transport across membranes",
    "F: Metabolic pathways",
    "G: Miscellaneous functions",
    "H: Protein folding/modification",
    "I: Signal transduction",
    "J: RNA processing",
    "K: Pathogenicity/defense mechanisms",
    "Z: Unclassified"
  )
)

# Join the mapping to your data
B2 <- B2 %>%
  left_join(superclass_info, by = "superclass")

library(cowplot)
library(grid)
library(ggplot2)

# Your existing plot
plot3 <- ggplot(B2, aes(x = class, y = genome_name, fill = n)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colors = c("white", "#ffeda0", "#feb24c", "#850101")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_test() +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 12, margin = margin(1, 1, 1, 1)),
    strip.placement = "outside", 
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust =1, vjust=0.5, size=11, face="bold"),
    axis.text.y = element_text(size=9, face="bold"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
    axis.title = element_blank(),
    plot.margin = margin(t=10, r=10, l=10, b=10),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  ) +
  facet_grid(~ superclass_label, scales = "free_x", space = "free_x")
plot3 <-plot3 + labs(
  x = "Subsystem Classes",
  y = "Strains")

# Description title and content
desc_title <- ggdraw() +
  draw_text("Superclass", x = 0, hjust = 0, fontface = "bold", size = 14)

desc_text <- paste(superclass_info$description, collapse = "\n")
desc_panel <- ggdraw() +
  draw_text(desc_text, x = 0, y = 1, hjust = 0, vjust = 1, size = 12, lineheight = 1.2, fontface = "italic")

# Combine title + text vertically
desc_column <- plot_grid(desc_title, desc_panel, ncol = 1, rel_heights = c(0.1, 0.9))

# Combine plot and the labeled description
final_plot <- plot_grid(plot3, desc_column, ncol = 2, rel_widths = c(4, 1))
final_plot
# 3. Save with dynamic sizing
ggsave(
  filename = "figure4.tiff",
  device = "tiff",
  width = 18,  # Adjust width to better fit x-axis labels
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)










#######Plot4
library(tidyverse)
B3 <-B1%>%
  group_by(genome_name, subclass) %>%
  count()
# 1. Complete missing combinations with zeros (only genome_name × subclass)
B3_complete <- B3 %>%
  ungroup() %>%
  tidyr::complete(genome_name, subclass, fill = list(n = 0))

plot4 <-ggplot(B3_complete, aes(x = subclass, y = genome_name, fill = n)) +
  geom_tile(color = "white", linewidth = 0.3)+
  scale_fill_gradientn(
    colors = c("white", "#ffeda0", "#feb24c", "#850101")) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove extra padding
  scale_y_discrete(expand = c(0, 0)) +
  theme_test() +
  theme(
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 15)),
    plot.margin = margin(t= 10, r = 20, l = 0, b =9),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )+
  labs(
    x = "Subsystem Subclasses",
    y = "Strains")

# 3. Save with dynamic sizing
ggsave(
  filename = "figure5.tiff",
  device = "tiff",
  width = 18,  # Adjust width to better fit x-axis labels
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)


library(patchwork)
###Combined_plot
combined_plot <- plot3 / plot4 
# Save
ggsave("combined_vertical.tiff", combined_plot, 
       width = 18, height = 20, dpi = 300)








####Plot5
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggdendro)
library(reshape2)
library(ggalign)


C1 <- read.csv("antismash.csv", check.names = FALSE)  

C2 <- C1 %>%
  pivot_longer(
    cols = 4:35,  # Columns 4 through 30 (inclusive)
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    feature = case_when(
      between(match(variable, names(C1)), 4, 13) ~ "host (+/-)",
      between(match(variable, names(C1)), 14, 22) ~ "BGC (+/-)",
      between(match(variable, names(C1)), 23, 29) ~ "knowncluster (identity)",
      between(match(variable, names(C1)), 30, 34) ~ "secretion systems (+/-)",
      match(variable, names(C1)) == 35 ~ "P",
      TRUE ~ NA_character_
    )
  )

ggplot(C2, aes(x = variable, y = Strain, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3)+
  scale_fill_gradientn(
    colors = c("white", "#ffeda0", "#feb24c", "#850101")) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove extra padding
  scale_y_discrete(expand = c(0, 0)) +
  theme_test() +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 13),
    strip.placement = "outside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )+ 
  labs(y = "Strains")+
  facet_grid(~ feature, scales = "free_x", space = "free_x")


ani_data <- read.table("ANIclustermap_matrix.tsv", header = TRUE)

# Get genome IDs from column names
genome_ids <- colnames(ani_data)

# Create new data frame with genome IDs as first column
ani_data_with_ids <- data.frame(
  Genome1 = genome_ids,  # New column with genome IDs
  ani_data,
  row.names = NULL,    # Remove automatic row numbers
  check.names = FALSE  # Preserve original column names
)
C3 <- ani_data_with_ids %>%
  pivot_longer(
    cols = 2:27,  # Columns 4 through 30 (inclusive)
    names_to = "Genome2",
    values_to = "ANI"
  )

C3$Genome1 <- sub("^(GCA_\\d+\\.\\d+).*", "\\1", C3$Genome1)
C3$Genome2 <- sub("^(GCA_\\d+\\.\\d+).*", "\\1", C3$Genome2)

C4 <- left_join(
  C2, 
  C3, 
  by = c("Sequences" = "Genome1")  # Match C1$Sequences with C2$Genome1
)
# Create a lookup table from C4 (Sequences → Strain)
strain_lookup <- C4 %>%
  distinct(Sequences, Strain)  # Remove duplicate pairs if they exist

# Replace Genome2 in C4 with matching Strain values
C4_updated <- C4 %>%
  left_join(strain_lookup, by = c("Genome2" = "Sequences")) %>%
  mutate(
    Genome2 = ifelse(is.na(Strain.y), Genome2, Strain.y)  # Replace if match found
  ) %>%
  select(-Strain.y)  

heatmap <- ggplot(C4_updated, aes(x = variable, y = Strain.x, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3)+
  scale_fill_gradientn(
    colors = c("white", "#ffeda0", "#feb24c", "#850101")) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove extra padding
  scale_y_discrete(expand = c(0, 0)) +
  theme_test() +
  theme(
    strip.background = element_rect(fill = "lightblue"),
    strip.text = element_text(size = 13),
    strip.placement = "outside", strip.clip = "off",
    panel.spacing.y = unit(1.5, "lines"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"), 
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(t = 15)),
    axis.title.x = element_blank(),
    plot.margin = margin(t= 10, r = 20, l = 0, b =9),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )+ labs(y = "Strains") +
  facet_grid(~ feature, scales = "free_x", space = "free_x")


ani_matrix <- C4_updated %>%
  select(Strain.x, Genome2, ANI) %>%
  distinct(Strain.x, Genome2, .keep_all = TRUE) %>%  # Keep only the first occurrence
  pivot_wider(
    names_from = Genome2, 
    values_from = ANI
  ) %>%
  column_to_rownames(var = "Strain.x")

rownames(ani_matrix) <- ani_matrix[[1]]
ani_matrix <- ani_matrix[-1]
distance_matrix <- 100 - as.matrix(ani_matrix)

# Hierarchical clustering
hc <- hclust(as.dist(distance_matrix), method = "average")

# Convert hclust object to a dendrogram for ggplot
dendro_data <- dendro_data(hc)

# Plot
dendro_plot <- ggdendrogram(dendro_data, rotate = TRUE) +
  coord_flip() +
  scale_y_reverse() +  
  labs(title = "FastANI Comparison") +
  theme_test() +
  theme(
    axis.text.y = element_text(size = 8, hjust = 0),    # Align to right
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.y.right = element_text(size = 8, hjust = 0),  # Forces text to appear on right
    plot.margin = margin(t = 10, r = 0, b = 200, l = 20)  # ← Add space below (b = bottom)
  ) +
  scale_x_discrete(position = "right")

# Ensure correct order for alignment
strain_order <- hc$labels[hc$order]
C4_updated$Strain.x <- factor(C4_updated$Strain.x, levels = strain_order)

library(cowplot)
combined_plot <- plot_grid(
  dendro_plot,
  heatmap,
  align    = "h",   # Horizontal alignment
  axis     = "l",   # Align left Y-axis (shared labels)
  rel_widths = c(0.2, 0.8)  # Adjust widths
)

# Save with publication quality
ggsave(
  filename = "figure6.tiff",
  device = "tiff",
  width = 20,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)













####Filtering data
library(dplyr)
T4SS <- B1 %>%
  filter(subsystem_name %in% c('Type 4 secretion and conjugative transfer', 'Type 4 secretion system')) %>%
  group_by(genome_name, gene, subsystem_name) %>%
  count()

T6SS <- B1 %>%
  group_by(genome_name, gene, subsystem_name) %>%
  filter(subsystem_name == 'Type VI secretion system') %>%
  count()
T4P <- B1 %>%
  group_by(genome_name, gene, subsystem_name) %>%
  filter(subsystem_name == 'Type IV pilus') %>%
  count()
Flagellum <- B1 %>%
  group_by(genome_name, gene, subsystem_name) %>%
  filter(subsystem_name == 'Flagellum') %>%
  count()
T1SS <- B1 %>%
  group_by(genome_name, gene, subsystem_name) %>%
  filter(subsystem_name == 'Type I secretion systems disambiguation') %>%
  count()

T2SS <- B1 %>%
  group_by(genome_name, gene, subclass) %>%
  filter(subclass == 'Protein secretion system, Type II') %>%
  count()

Xenobiotics <- A1 %>% 
  group_by(pathway_class, pathway_name) %>%
  filter(pathway_class == "Xenobiotics Biodegradation and Metabolism") %>%
  count()



#####Not useful

B2 <-B1%>%
  group_by(genome_name, superclass, class) %>%
  count()%>%
  mutate(superclass = tolower(superclass))
ggplot(B2, aes(x = class, y = n, fill = class)) +  
  geom_violin(alpha = 0.6, trim = FALSE) +  
  geom_jitter(aes(color = class), width = 0.2, size = 2, alpha = 0.8) +  
  scale_fill_viridis_d(option = "C") +  
  scale_color_viridis_d(option = "C") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none"  # Remove legend
  ) +
  facet_wrap(~ superclass, scales = "free_x", nrow = 2)

ggplot(B2, aes(x = class, y = genome_name, fill = n)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "grey", high = "black", midpoint = 50) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


B3 <-B1%>%
  group_by(genome_name, superclass, subclass) %>%
  count()%>%
  mutate(superclass = tolower(superclass))

ggplot(B3, aes(x = subclass, y = n, fill = subclass)) +
  geom_violin(alpha = 0.6, trim = FALSE) +  
  geom_jitter(aes(color = subclass), width = 0.1, size = 2, alpha = 0.8) +  
  scale_fill_viridis_d(option = "C") +  
  scale_color_viridis_d(option = "C") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 9),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none"  # Remove legend
  ) +
  facet_grid(~ superclass, scales = "free_x", space = "free_x")
  facet_wrap(~ superclass, scales = "free_x", nrow = 2)
  
B3 <-B1%>%
  group_by(genome_name, superclass, subclass) %>%
  count()%>%
  mutate(superclass = tolower(superclass)) 
ggplot(B3, aes(x = subclass, y = genome_name, fill = n)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "Red", midpoint = 50) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


B4 <-B1%>%
  group_by(genome_name, subsystem_name) %>%
  count()
ggplot(B4, aes(x = subsystem_name, y = genome_name, fill = n)) +
  geom_tile() +
  scale_fill_gradient2(low = "cyan", mid = "white", high = "red", midpoint = 50)
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

ggplot(B3, aes(x = subsystem_name, y = n, fill = subclass)) + 
  geom_jitter(aes(color = subsystem_name), width = 0.2, size = 2, alpha = 0.8) +
  geom_violin() +
  scale_fill_viridis_d() +  # Automatically generates distinct colors
  scale_color_viridis_d() +
  theme_test() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5), legend.position = "none") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = function(x) stringr::str_wrap(x, width = 10))  # Wrap long labels

  geom_violin(alpha = 0.6, trim = FALSE) +  
  geom_jitter(aes(color = subsystem_name), width = 0.1, size = 2, alpha = 0.8) +  
  scale_fill_viridis_d(option = "C") +  
  scale_color_viridis_d(option = "C") +
  theme_test() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none"  # Remove legend
  ) +
  facet_wrap(~ superclass, scales = "free_x", nrow = 4) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))


  
  
 







library(tidyverse)
B4 <-B1%>%
  group_by(genome_name, subsystem_name) %>%
  count()
# 1. Complete missing combinations with zeros (only genome_name × subsystem_name)
B4_complete <- B4 %>%
  ungroup() %>%
  tidyr::complete(genome_name, subsystem_name, fill = list(n = 0))

ggplot(B4_complete, aes(x = subsystem_name, y = genome_name, fill = n)) +
  geom_tile(color = "white", linewidth = 0.3) +  # Add white borders
  scale_fill_gradientn(
    colors = c("white", "#ffeda0", "#feb24c", "#850101"),  # Yellow-orange-red
    na.value = "grey90",
    limits = c(0, max(B4_complete$n)),
    name = "Gene Count"
  ) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove extra padding
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Functional subsystem_name",
    y = "Genome"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  )

# 3. Save with dynamic sizing
ggsave(
  filename = "plot5.tiff",
  device = "tiff",
  width = 20,  # Adjust width to better fit x-axis labels
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)








