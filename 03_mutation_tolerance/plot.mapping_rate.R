#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)

infile <- "mapping_rate.tsv"
outfile <- "mapping_rate.svg"
w <- 8
h <- 6
dpi <- 300

df <- read.csv(infile, sep = "\t")

df$tool <- factor(df$tool, levels = unique(df$tool), ordered = TRUE)
df$qlen <- factor(df$qlen, levels = unique(df$qlen), ordered = TRUE)

p <-
  ggplot(df, aes(
    x = ident,
    y = recall,
    color = tool,
    # shape = tool,
    # color = genome_size/1000000,
  )) +
  
  geom_quasirandom(dodge.width = 0.6, size = 0.8) +
  geom_vline(xintercept = c(85,90,95), linewidth = 0.3, linetype=2, color = "grey70") + 
  #geom_hline(yintercept = c(90,95), linewidth = 0.3, linetype=2, color = "grey70") +
  # scale_color_stata() +
  scale_color_colorblind() +
  # scale_color_continuous_tableau("Red") + 
  # ylim(80, 100) +
  
  facet_wrap(qlen ~ .) +
  ylab("Alignment rate (%)") +
  xlab("Query identity (%)") +
  guides(color = guide_legend(title="", override.aes = list(size = 2))) +
  labs(color = "Tools")


for (a in unique(df$assembly)) {
  p <- p + geom_line(data = filter(df, assembly==a),
                     aes(x = ident, y = recall, color = tool),
                     linewidth = 0.3, alpha = 0.3)
}

p <- p +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "grey20", size = 0.8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(size = 0.6),
    axis.ticks.x = element_line(size = 0.6),
    
    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    
    legend.text = element_text(size = 11),
    # legend.position = "right",
    legend.position = c(0.9, 0.15),
    legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.6, "cm"),
    # legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    # legend.spacing.y = unit(0.1, "cm"),
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
