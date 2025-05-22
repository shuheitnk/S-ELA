library(dplyr)
library(reshape2)
library(ggplot2)

setwd("path")

files <- c("accuracy_score_medium_25d.csv",
           "accuracy_score_medium_50d.csv",
           "accuracy_score_large_25d.csv",
           "accuracy_score_large_50d.csv",
           "accuracy_score_s_ela.csv")

dfs <- lapply(files, function(file) {
  df <- read.csv(file)
  suffix <- tools::file_path_sans_ext(file)
  colnames(df)[2] <- suffix
  df
})

merged_df <- Reduce(function(x, y) full_join(x, y, by = "run"), dfs)

accuracy_long <- melt(merged_df, id.vars = "run", variable.name = "variable", value.name = "Value")

rename_map <- c(
  "accuracy_score_medium_25d" = "medium_25d",
  "accuracy_score_medium_50d" = "medium_50d",
  "accuracy_score_large_25d"  = "large_25d",
  "accuracy_score_large_50d"  = "large_50d",
  "accuracy_score_s_ela"      = "S-ELA"
)
accuracy_long$Category <- rename_map[accuracy_long$variable]

accuracy_long$Category <- factor(accuracy_long$Category, levels = c(
  "medium_25d",
  "medium_50d",
  "large_25d",
  "large_50d",
  "S-ELA"
))

ggplot(accuracy_long, aes(x = Category, y = Value, fill = Category)) +
  stat_boxplot(geom = "errorbar", width = 0.2, size = 1, color = "#31518C") +
  geom_boxplot(outlier.color = "#909089", color = "#31518C", linewidth = 1, alpha = 0.7) +
  labs(x = NULL, y = "Accuracy Score", title = NULL) +
  theme(
    text = element_text(family = "serif", margin = margin(r = 10)),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 21),
    axis.text.y = element_text(size = 21),
    axis.title.y = element_text(size = 25, margin = margin(r = 15)),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 1),
    panel.grid.major = element_line(color = "gray80", size = 0.5),
    panel.grid.minor = element_line(color = "gray90", size = 0.25)
  ) +
  scale_y_continuous(limits = c(0.73, 0.82), breaks = seq(0.73, 0.82, by = 0.02)) +
  scale_fill_manual(values = c(
    "medium_25d" = "grey",
    "medium_50d" = "grey",
    "large_25d"  = "grey",
    "large_50d"  = "grey",
    "S-ELA"      = "#F87B12"
  ))


# image size: width = 1100px, height = 400px



