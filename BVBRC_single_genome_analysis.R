setwd("/home/anik/genomics/Projects/prova")

f1 <- read.csv("BVBRC_pathway.csv") %>%
  filter(genome_name == "Acinetobacter baumannii DUEMBL6")
f2 <- read.csv("BVBRC_subsystem.csv") %>% 
  filter(genome_name == "Acinetobacter baumannii DUEMBL6")

library(dplyr)
library(ggplot2)
f1.1 <-f1%>%
  group_by(pathway_class, pathway_name) %>%
  count()

ggplot(f1.1, aes(x = pathway_name, y = n, fill = pathway_class)) +  
  # Create the lollipop sticks
  geom_segment(aes(x = pathway_name, xend = pathway_name, 
                   y = 0, yend = n), 
               color = "gray") +
  # Add the points (lollipop heads)
  geom_point(aes(color = pathway_class), size = 2, show.legend = FALSE) + 
  
  facet_wrap(
    ~ stringr::str_wrap(pathway_class, width = 30),
    scales = "free_x", 
    nrow = 2
  ) +
  
  # Improved aesthetics
  scale_color_viridis_d(option = "plasma", guide = "none") +  # Remove legend
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  # Professional theme adjustments
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gainsboro"),
    strip.text = element_text(size = 8, face = "bold", color = "black"),
    strip.placement = "outside", 
    strip.clip = "off",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
    axis.title.x = element_text(face = "bold", color = "black"),  # Darker x-axis label
    axis.title.y = element_text(face = "bold", color = "black"),  # Darker y-axis label
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
    legend.position = "none"  # Ensure legend is removed
  ) +
  labs(
    x = "Pathway Class",
    y = "Gene Count"
  )

# Save with publication quality
ggsave(
  filename = "figure2.tiff",
  device = "tiff",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)
  
  
f1.2 <- f1 %>% 
  filter(pathway_class == "Lipid Metabolism") %>% 
  group_by(pathway_name, ec_description) %>%
  count()

ggplot(f1.2, aes(x = ec_description, y = n, fill = pathway_name)) +  
  geom_col(width = 0.7, show.legend = FALSE) +  # Barplot instead of lollipop
  
  facet_wrap(
    ~ stringr::str_wrap(pathway_name, width = 30),
    scales = "free_x", 
    nrow = 3
  ) +
  
  # Wrap x-axis labels to avoid overlap
  scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 30)  
  ) +
  
  # Nice color palette
  scale_fill_viridis_d(option = "plasma") +
  
  # Y-axis expansion for breathing room
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  # Theme polishing
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gainsboro"),
    strip.text = element_text(size = 9, face = "bold", color = "black"),
    strip.placement = "outside",
    strip.clip = "off",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
    axis.title.x = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
    legend.position = "none"
  ) +
  labs(
    x = "Product",
    y = "Gene Count"
  )

# Save high-res figure
ggsave(
  filename = "figure4.tiff",
  device = "tiff",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)





f1.3 <- f1 %>% 
  filter(pathway_class == "Xenobiotics Biodegradation and Metabolism") %>% 
  mutate(
    pathway_name = case_when(
      pathway_name == "1,1,1-Trichloro-2,2-bis(4-chlorophenyl)ethane (DDT) degradation" ~ "DDT degradation",
      TRUE ~ pathway_name  # Keep other descriptions unchanged
    )
  ) %>% 
  group_by(pathway_name, ec_description) %>%
  count()

library(readr)
write_csv(f1.3, "biodegradation.csv")

ggplot(f1.3, aes(x = ec_description, y = n, fill = pathway_name)) +  
  geom_col(width = 0.7, show.legend = FALSE) +  # Barplot instead of lollipop
  
  facet_wrap(
    ~ stringr::str_wrap(pathway_name, width = 25),
    scales = "free_x", 
    nrow = 4
  ) +
  
  # Wrap x-axis labels to avoid overlap
  scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 30)  
  ) +
  
  # Force y-axis to show only integers
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    breaks = scales::pretty_breaks()  # Ensures integer breaks
  ) +
  
  # Nice color palette
  scale_fill_viridis_d(option = "plasma") +
  
  # Theme polishing
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gainsboro"),
    strip.text = element_text(size = 9, face = "bold", color = "black"),
    strip.placement = "outside",
    strip.clip = "off",
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.title.x = element_text(face = "bold", color = "black", size = 11),
    axis.title.y = element_text(face = "bold", color = "black", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm"),
    legend.position = "none"
  ) +
  labs(
    x = "Xenobiotics degrading enzymes",
    y = "Gene Count"
  )

# Save with publication quality
ggsave(
  filename = "figure_5.tiff",
  device = "tiff",
  width = 24,
  height = 14,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)


f2.1 <-f2%>%
  group_by(superclass, class, subclass) %>%
  count()%>%
  mutate(superclass = tolower(superclass))

ggplot(f2.1, aes(x = subclass, y = n, fill = superclass)) + 
  geom_bar(stat = "identity", width = 0.7) +  # Use stat="identity" for raw 'n' values
  
  facet_wrap(
    ~ stringr::str_wrap(class, width = 25),  # Wrap facet titles
    scales = "free_x",  # Independent x-axes per facet
    nrow = 3
  ) +
  # Wrap x-axis labels (products) to avoid overlap
  scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 25)  # Adjust width as needed
  ) +
  
  # Use a palette with 12 distinct colors (e.g., "Paired", "Set3", or "viridis")
  scale_fill_manual(
    values = RColorBrewer::brewer.pal(12, "Paired")
  )+ # Alternative: scale_fill_viridis_d(option = "turbo")  # 12+ distinct colors
  
  # Adjust y-axis padding
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  # Theme improvements
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gainsboro"),
    strip.text = element_text(size = 8, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.title = element_text(face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),  # Cleaner look
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",  # Legend below plot
    legend.text = element_text(size = 8),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(0.4, "cm"),
    panel.border = element_rect(color = "gray80", fill = NA, size = 0.5),
    panel.grid.major.y = element_line(color = "gray90")
  ) +
  
  labs(
    x = "Subclass",
    y = "Count",
    fill = "Superclass"  # Legend title
  )
  
# Save with publication quality
ggsave(
  filename = "figure3.tiff",
  device = "tiff",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

Xenobiotics <- f1.1 %>% 
  group_by(pathway_class, pathway_name) %>%
  filter(pathway_class == "Xenobiotics Biodegradation and Metabolism") %>%
  count()
