#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)

infile <- "aligned_fraction.tsv"
outfile <- "aligned_fraction.png"
w <- 5
h <- 3
dpi <- 300

df <- read.csv(infile, sep = "\t")

df$tool <- factor(df$tool, levels = unique(df$tool), ordered = TRUE)
df$qlen <- factor(df$qlen, levels = unique(df$qlen), ordered = TRUE)

p <-
  ggplot(df,
         aes(x = qlen,
             y = recall,
             color = tool,)) +
  
  geom_quasirandom(dodge.width = 0.7) +
  # scale_color_stata() +
  scale_color_colorblind() +
  ylim(80, 100) +
  
  # facet_grid(. ~ qlen) +
  ylab("Aligned fraction (%)") +
  xlab("Query length (bp)") +
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
    
    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.x = element_text(size = 11),
    strip.text.y = element_text(size = 11),
    
    legend.text = element_text(size = 11),
    # legend.position = c(args$lx,args$ly),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.6, "cm"),
    # legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    # legend.spacing.y = unit(0.1, "cm"),
    # space between key and text
    
    text = element_text(size = 11,
                        family = "arial"),
  )


if (grepl("tiff?$",
          outfile,
          perl = TRUE,
          ignore.case = TRUE)) {
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
