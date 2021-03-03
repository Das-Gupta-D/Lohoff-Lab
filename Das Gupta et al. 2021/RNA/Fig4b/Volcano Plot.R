install.packages("ggrepel")
install.packages("svglite")
library(tidyverse)
library(readr)
library(ggrepel)
library(svglite)

volcano <- read_delim("data/volcano.csv", ";", 
    escape_double = FALSE, col_types = cols(log2FoldChange = col_double()), 
    locale = locale(decimal_mark = ",", grouping_mark = "."), 
    trim_ws = TRUE
    )

volcano_top <- volcano %>%
  slice_max(order_by = `adjusted P -log10`, prop = .005) 

volcano_interest <- volcano %>%
  filter(geneSymbol %in% c("Ccnd3", "Igll1", "Vpreb2", "Igkv5-43", "Irf4", "Rag1", "Ebf1", "Vpreb1", "Pecam1", "Arid3b", "Ikzf1", "Ikzf3"))

volcano <- volcano %>%
  mutate(
    relevant = ifelse(
      ((log2FoldChange > 1) | (log2FoldChange < -1)) & (`adjusted P -log10` >= 5), TRUE, FALSE)
  )

volcano %>%
  count(relevant)

## Plot function for the volcano plot
ggplot(mapping = aes(x = log2FoldChange, y = `adjusted P -log10`)) +
  
  # general scatter plot with all the genes
  geom_point(data = volcano, 
             aes(color = relevant), 
             size = 1,
             alpha = 1
             ) +
  
  # highlighting the most important genes
  geom_point(data = volcano_interest, 
             color = "#BA2B34", 
             alpha = 1, size = 1
             ) +
  
  # annotating the highlighted genes
  geom_text(data = volcano_interest, 
            aes(label = geneSymbol, x = log2FoldChange-0.1), 
            size = 4, hjust = 1, vjust = 0
            ) +
  
  # adjusting theme, removing grid, and color legend
  theme_bw() +
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      ) +
  guides(color = FALSE) +
  
  # adjusting color pallet to fit the final graphic
  scale_color_manual(values = c("#A0A09F", "#A7CFEE")) +
  
  # adding x- and y-titles
  labs(y = "BH-adjusted P (-log10)",
       x = "log2(fold change)"
       ) +

  # adding reference lines for P value and fold change cutoff
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed")
  
# outputing the file
ggsave("volcano.svg", device = "svg", width = 8, height = 8, units = "cm")
