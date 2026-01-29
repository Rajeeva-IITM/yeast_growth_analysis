# Set theming options for consistency across figures 

library(ggplot2)
library(purrr)

save_plot <- partial(ggsave, dpi=300, bg="transparent")
# Colorblind palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
paper_theme <- theme_set(theme_bw(base_size = 12))
theme_update(
  text = element_text(size = 12),
  plot.title = element_text(hjust = 0.5, size=14),
  axis.text = element_text(size=11),
  strip.background = element_rect(color = "grey0", fill = "grey"),
  strip.text = element_text(colour = "black", hjust = 0.5, size = 12),
  plot.tag = element_text(size = 15,)
)