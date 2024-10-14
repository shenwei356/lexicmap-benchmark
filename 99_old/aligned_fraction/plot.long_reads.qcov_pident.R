#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)

infile <- "qcov_pident.tsv"
outfile <- "qcov_pident.png"
w <- 9
h <- 4
dpi <- 300

df <- read.csv(infile, sep = "\t")


t <- theme_bw() +
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


p1 <-
  ggplot(df,
         aes(x = qcov_blastn,
             y = qcov_lexicmap,
             color = qlen, )) +
  geom_segment(
    x = 0,
    y = 0,
    xend = 200,
    yend = 200,
    color = "grey70",
    size = 0.3
  ) +
  geom_point(alpha = 0.7) +
  scale_color_gradient_tableau(palette = "Blue", trans = "log10") +
  xlim(70, 100) +
  ylim(70, 100) +
  ylab("Query coverage (LexicMap)") +
  xlab("Query coverage (Blastn)") +
  labs(color = "log10(query len)") +
  t

p2 <-
  ggplot(df,
         aes(x = pident_blastn,
             y = pident_lexicmap,
             color = qlen, )) +
  geom_segment(
    x = 0,
    y = 0,
    xend = 200,
    yend = 200,
    color = "grey70",
    size = 0.3
  ) +
  geom_point(alpha = 0.7) +
  scale_color_gradient_tableau(palette = "Blue", trans = "log10") +
  xlim(60, 100) +
  ylim(60, 100) +
  ylab("Percentage identity (LexicMap)") +
  xlab("Percentage identity (Blastn)") +
  labs(color = "log10(query len)") +
  t

legend <- get_legend(# create some space to the left of the legend
  p1)

p <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  legend,
  ncol = 3,
  labels = c("a", "b"),
  rel_widths = c(1, 1, 0.35)
) + theme (# fill the gap in sub figures
  panel.background = element_rect(fill = "white", colour = NA),)


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
