#!/usr/bin/env Rscript
library(tidyr)
library(ggplot2)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)

infile <- "marker_gene.tsv"
outfile <- "marker_gene.png"
w <- 7
h <- 3
dpi <- 300

df0 <- read.csv(infile, sep = "\t")

df <- df0 %>% select(-query) %>% gather(key="tool", value="recall")
df$tool <- factor(df$tool, levels = unique(df$tool), ordered = TRUE)


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
         aes(x = tool,
             y = recall,
             color = tool)) +
  
  geom_quasirandom() +
  scale_color_colorblind() +
  ylim(80, 100) +
  ylab("Recall (%)") +
  xlab(NULL) +
  labs(color = "Tools") +
  t + theme(
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
  

p2 <-
  ggplot(df0,
         aes(x = Blastn,
             y = LexicMap)) +
  geom_segment(
    x = 0,
    y = 0,
    xend = 200,
    yend = 200,
    color = "grey70",
    size = 0.3
  ) +
  geom_point(alpha = 0.7, color = "black") +
  xlim(98, 100) +
  ylim(98, 100) +
  ylab("Recall (LexicMap)") +
  xlab("Recall (Blastn)") +
  t

p <- plot_grid(
  p1,
  p2,
  ncol = 2,
  labels = c("a", "b"),
  rel_widths = c(1, 1)
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
