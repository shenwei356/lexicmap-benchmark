#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(ggthemes)
library(ggbeeswarm)
library(cowplot)

infile <- "pt.mapping_rate.tsv"
outfile <- "pt.mapping_rate_index_size.svg"
w <- 10
h <- 10
dpi <- 300

df <- read.csv(infile, sep = "\t")

df <- df %>% mutate(comb = paste0(window, "-", density)) %>%
  mutate(qlen = paste0("Query length: ", qlen, " bp")) %>%
  mutate(probes = paste0("Probes: ", probes))

df$comb <- factor(df$comb, levels = unique(df$comb), ordered = TRUE)
df$qlen <- factor(df$qlen, levels = unique(df$qlen), ordered = TRUE)
df$probes <- factor(df$probes, levels = unique(df$probes), ordered = TRUE)

df <- df %>% group_by(ident,qlen,comb,probes) %>% summarise(rate=mean(recall))

p1 <-
  ggplot(df, aes(
    x = ident,
    y = rate,
    color = comb,
  )) +

  geom_vline(xintercept = c(80, 85,90,95), linewidth = 0.3, linetype=2, color = "grey70") +
  geom_hline(yintercept = c(0, 25, 50, 75), linewidth = 0.3, linetype=2, color = "grey70") +
  geom_line() +
  scale_color_colorblind() +

  facet_grid(qlen ~ probes) +
  ylab("Average alignment rate (%)") +
  xlab("Query identity (%)") +
  labs(color = "Maximum seed distance\nand spacing of seeds\nadded to deserts")

t <- theme_bw() +
  theme(
    panel.border = element_rect(color = "grey20", size = 0.8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(size = 0.6),
    axis.ticks.x = element_line(size = 0.6),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),

    text = element_text(size = 11, family = "arial"),

    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),

    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    # legend.position = "right",
    legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.6, "cm"),
    # legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    legend.spacing.y = unit(0, "cm"),
    # space between key and text

  )

p1 <- p1 + t +
  theme(legend.title = element_text(size = 11),
        legend.position = c(0.9, 0.153),)

# -------------------------------------------------------------------
# index size

infile2 <- "pt.index_info.tsv"

df2 <- read.csv(infile2, sep = "\t")

df2 <- df2 %>% mutate(comb = paste0(window, "-", density)) %>%
  mutate(probes = paste0("Probes: ", probes))

df2$comb <- factor(df2$comb, levels = unique(df2$comb), ordered = TRUE)
df2$probes <- factor(df2$probes, levels = unique(df2$probes), ordered = TRUE)

df2 <- df2 %>%
    group_by(comb,probes) %>%
    summarise(asize=mean(size), amem=mean(mem), atime=mean(time))

p2 <-
  ggplot(df2, aes(
    x = atime,
    y = asize/1000000,
    color = comb,
    size = amem/1000,
  )) +

  geom_vline(xintercept = c(0, 2, 4, 6), linewidth = 0.3, linetype=2, color = "grey70") +
  geom_hline(yintercept = c(0, 3, 6, 9, 12), linewidth = 0.3, linetype=2, color = "grey70") +
  geom_point() +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  facet_grid(. ~ probes) +
  xlim(0, max(df2$atime)+0.3) +
  ylim(0, max(df2$asize/1000000)+0.6) +
  ylab("Average index size (MB)") +
  xlab("Average indexing time (second)") +
  t +
  labs(color = "Maximum seed distance\nand spacing of seeds\nadded to deserts",
       size = "Average indexing\nmemory(MB)") +
  theme(legend.position = "right",
        # axis.text.x = element_text(
        #   angle = 30, hjust = 1, vjust = 1
        # ),
        legend.title = element_text(size = 11),
        legend.box = "horizontal",
        legend.box.spacing = unit(0.2, "cm"),
        legend.margin = unit(0.1, "cm"),
  )

# -------------------------------------------------------------------

p <- plot_grid(
  p1,
  p2,
  nrow = 2,
  labels = c("a", "b"),
  rel_heights = c(3, 1.1)
) + theme ( # fill the gap in sub figures
  panel.background = element_rect(fill = "white", colour = NA),
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
