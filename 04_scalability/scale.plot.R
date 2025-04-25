#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(ggthemes)
library(cowplot)

t <- theme_bw() +
  theme(
    panel.border = element_rect(color = "grey20", size = 0.8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(size = 0.6),
    axis.ticks.y = element_line(size = 0.6),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    
    strip.background = element_rect(
      colour = "white",
      fill = "grey95",
      size = 0.2
    ),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    
    legend.text = element_text(size = 11),
    # legend.position = "right",
    legend.background = element_rect(fill = "transparent"),
    # legend.key.size = unit(0.6, "cm"),
    # legend.key = element_blank(),
    legend.text.align = 0,
    legend.box.just = "left",
    # space between key and text
    legend.spacing.y = unit(0, "mm"),
    
    text = element_text(size = 12, family = "arial"),
  )

# -----------------------------------------------------------------------------

w <- 8
h <- 6
dpi <- 300

outfile <- "scalability.svg"


infile <- "indexing.tsv"
df <- read.csv(infile, sep = "\t")

infile2 <- "searching.tsv"
df2 <- read.csv(infile2, sep = "\t")

# tmp <- df %>% filter(genomes == 100000) %>% 
#   arrange(desc(index_size))

df$tool <- factor(df$tool, levels = unique(df2$tool), ordered = TRUE)

colors <- colorblind_pal()(8)

tools1 <- unique(df$tool)
tools2 <- unique(df2$tool)
colors2 <- colors[1:length(tools2)]
dfc <- data.frame(tool = tools2, color = colors2 )
colors1 <- dfc %>% filter(tool %in% tools1) %>% select(color) %>% unlist() %>% unname()


df$task <- rep("Indexing", each=nrow(df))


p1 <-
  ggplot(df, aes(
    x = genomes,
    y = index_size/1000000000,
    color = tool,
  )) +
  geom_line(size=0.7) +
  geom_point(aes(size = indexing_mem/1000), alpha = 0.7) +
  # scale_color_colorblind() +
  scale_color_manual(values = colors1) + 
  xlab("Number of genomes") +
  ylab("Index size (GB)") +
  scale_x_log10(breaks=10^(0:6),
                labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=10^(-3:4),
                labels=trans_format("log10", math_format(10^.x))) +
  guides(color = guide_legend(title="Tool", order = 1),
         size = guide_legend(title="Indexing memory (GB)"), order = 2) +
  facet_grid(. ~ task) +
  t +
  theme(legend.position = c(0.27, 0.71))

# -----------------------------------------------------------------------------



# tmp <- df2 %>% filter(genomes == 100000) %>% 
#   arrange(desc(seconds))
# df2$tool <- factor(df2$tool, levels = append(tmp$tool,"Blastn(ws=15)"), ordered = TRUE)

df2 <- df2 %>% group_by(genomes, tool) %>%
    summarise(mean_seconds=mean(seconds),
             mean_kb=mean(kb))

df2$task <- rep("Searching", each=nrow(df2))

log60_trans <- function() {
  trans_new(
    name = "log60",
    transform = function(x) log(x, 60), 
    inverse = function(x) 60^x
  )
}

p2 <-
  ggplot(df2, aes(
    x = genomes,
    y = mean_seconds,
    color = tool,
  )) +
  geom_line(size=0.7) +
  geom_point(aes(size = mean_kb/1000000), alpha = 0.7) +
  # scale_color_colorblind() +
  scale_color_manual(values = colors2) + 
  xlab("Number of genomes") +
  ylab("Searching time") +
  scale_x_log10(breaks=10^(0:6),
                labels=trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    trans = log60_trans(),
    breaks = 60^c(0, 1, 2, 2.776206),
    labels = function(x) {
      case_when(
        as.integer(x) == 1 ~ "1 second",
        as.integer(x) == 60 ~ "1 minute",
        as.integer(x) == 3600 ~ "1 hour",
        as.integer(x) == 3600*24 ~ "1 day",
        TRUE ~ as.character(x)
      )
    }
  ) +
  guides(color = guide_legend(title="Tool", order = 1),
         size = guide_legend(title="Searching memory (GB)", order = 2)) +
  facet_grid(. ~ task) +
  t +
  theme(legend.position = c(0.31, 0.71))



# -----------------------------------------------------------------------------

p <- plot_grid(
  p1,
  p2,
  ncol = 2,
  labels = c("a", "b"),
  rel_widths = c(1, 1)
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
