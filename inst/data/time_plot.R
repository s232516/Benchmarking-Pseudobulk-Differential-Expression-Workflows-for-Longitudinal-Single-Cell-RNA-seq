library(ggplot2)

runtime_df <- data.frame(
  dataset = c(
    "splatPop","splatPop","splatPop",
    "rescueSim","rescueSim","rescueSim",
    "rescueSim (norm)","rescueSim (norm)","rescueSim (norm)"
  ),
  method = c(
    "MAST fixed-effects model","MAST mixed-effects model","DESeq2",
    "MAST fixed-effects model","MAST mixed-effects model","DESeq2",
    "MAST fixed-effects model","MAST mixed-effects model","DESeq2"
  ),
  time_s = c(
    3.61, 29.55, 3.08,
    46.35, 6184.15, 1.85,
    46.44, 7839.68, 2.01
  )
)

runtime_df$dataset <- factor(runtime_df$dataset,
                             levels = c("splatPop","rescueSim","rescueSim (norm)"))

pastel_cols <- c(
  "MAST fixed-effects model" = "#79CDCD",
  "MAST mixed-effects model" = "#FFAEB9",
  "DESeq2"                   = "#9BCD9B"
)

p <- ggplot(runtime_df, aes(x = dataset,
                            y = time_s,
                            group = method,
                            color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = pastel_cols) +
  scale_y_log10() +  # <- remove this for linear scale
  labs(
    title = "Runtime by method and simulation",
    x = "Dataset",
    y = "Runtime (seconds)",
    color = "Method"
  ) +
  theme_bw(base_size = 29) +
  theme(
    plot.title   = element_text(size = 32),
    axis.title.x = element_text(size = 29),
    axis.title.y = element_text(size = 29),
    axis.text.x  = element_text(size = 29, angle = 20, hjust = 1),
    axis.text.y  = element_text(size = 29),
    legend.title = element_text(size = 29),
    legend.text  = element_text(size = 29),
    strip.text   = element_text(size = 29),
    aspect.ratio = 1,
    legend.position = "right"
  )


p

# optional save
# ggsave("runtime_lineplot.png", p, width = 8, height = 6, dpi = 300)
