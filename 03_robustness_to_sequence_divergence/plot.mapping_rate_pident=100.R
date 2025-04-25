#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)

infile <- "mapping_rate.tsv"
outfile <- "mapping_rate_pident100.svg"
w <- 8
h <- 6
dpi <- 300

df <- read.csv(infile, sep = "\t") %>% filter(ident == 100)

df$tool <- factor(df$tool, levels = unique(df$tool), ordered = TRUE)
df$qlen <- factor(df$qlen, levels = unique(df$qlen), ordered = TRUE)

df2 <- df %>% group_by(ident, qlen, tool) %>%
    summarise(recall_mean = mean(recall), recall_max = max(recall))

#  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors <- colorblind_pal()(8)
n <- length(unique(df$tool))

p <-
  ggplot(df, aes(
    x = tool,
    y = recall,
    color = tool,
    # shape = tool,
    # color = genome_size/1000000,
  )) +
  
  geom_quasirandom(dodge.width = 0.1, size = 0.8) +
  geom_text(data = df2, aes(x=tool, y=recall_max + 0.15, 
                            label = sprintf("%.4f", recall_mean)),
            angle=30, hjust=0.3, vjust=0, size = 3.5) + 
  
  scale_color_manual(values = colors[1:n]) +
  ylim(98.5, 100.4) +
  
  facet_wrap(qlen ~ .) +
  ylab("Alignment rate (%)") +
  xlab(NULL) +
  # guides(color = guide_legend(title="", override.aes = list(size = 2))) +
  labs(color = "Tools")


p <- p +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "grey20", size = 0.8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(size = 0.6),
    axis.ticks.x = element_line(size = 0.6),
    axis.text.x = element_text(size = 11, angle=30, hjust=1, vjust=1),
    
    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    
    legend.text = element_text(size = 11),
    # legend.position = "right",
    legend.position = "none", # c(0.91, 0.2),
    legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.6, "cm"),
    # legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    # space between key and text
    
    text = element_text(size = 11, family = "arial"),
  )


if (grepl("tiff?$", outfile, perl = TRUE, ignore.case = TRUE)) {
  ggsave(
    p,
    file = outfile,
    width = w,
    height = h,
    dpi = dpi,
    compress = "lzw"
  )
} else {
  ggsave(
    p,
    file = outfile,
    width = w,
    height = h,
    dpi = dpi
  )
}
