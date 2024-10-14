#!/usr/bin/env Rscript
library(ggplot2) 

infile <- "prefix.freq.tsv"
outfile <- "prefix.freq.svg"

w <- 3.5
h <- 3
dpi <- 300

df <- read.csv(infile, sep = "\t")

p <-
  ggplot(df, aes(
    x = prefix,
    y = freq,
  )) +
  
  geom_col( fill = "grey35") +
  
  ylab("Frequency") +
  xlab("Length of the shared prefix\nbetween a probe and captured k-mers") 


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
    legend.position = "none",
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
