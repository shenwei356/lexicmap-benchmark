#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(ggthemes)

infile <- "pt2.mapping_rate.lexicmap.tsv"
outfile <- "pt2.mapping_rate.lexicmap.svg"
w <- 6
h <- 7
dpi <- 300

df <- read.csv(infile, sep = "\t")

df <- df %>% mutate(qlen = paste0("Query length: ", qlen, " bp")) 

df$wordsize <- factor(df$wordsize, levels = unique(df$wordsize), ordered = TRUE)
df$qlen <- factor(df$qlen, levels = unique(df$qlen), ordered = TRUE)

df <- df %>% group_by(ident,qlen,wordsize) %>% summarise(rate=mean(recall))

p <-
  ggplot(df, aes(
    x = ident,
    y = rate,
    color = wordsize,
  )) +
  
  geom_line() +
  geom_vline(xintercept = c(80, 85,90,95, 100), linewidth = 0.3, linetype=2, color = "grey70") + 
  geom_hline(yintercept = c(0, 25, 50, 75), linewidth = 0.3, linetype=2, color = "grey70") +
  # scale_color_stata() +
  scale_color_colorblind() +
  # scale_color_continuous_tableau("Red") + 
  # ylim(80, 100) +
  
  facet_grid(qlen ~ .) +
  ylab("Average alignment rate (%)") +
  xlab("Query identity (%)") +
  guides(color = guide_legend(title="Minimum lengths of\nmultiple and single anchors", override.aes = list(size = 2)))


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
    legend.position = c(0.8, 0.15),
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
