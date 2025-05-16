generate_plot <- function(b1_file, b2_file, res_df, output_file = "venn_bar_combined.png") { 
  
  # Load data
  b1 <- read.csv(b1_file, stringsAsFactors = FALSE)
  b2 <- read.csv(b2_file, stringsAsFactors = FALSE)
  
  # Set label from file name
  bioset1_label <- tools::file_path_sans_ext(b1_file)
  bioset2_label <- tools::file_path_sans_ext(b2_file)
  
  # Format
  read_and_format <- function(df) {
    df <- df[!duplicated(df$Gene), c("Gene", "Fold.Change")]
    colnames(df) <- c("Gene", "FC")
    df$Gene <- tolower(trimws(df$Gene))
    df
  }
  b1 <- read_and_format(b1)
  b2 <- read_and_format(b2)
  
  
  df <- res_df[res_df$Pair != "Total", ]
  df$logP <- -log10(df$Pvalue)
  df$Pair <- factor(df$Pair, levels = c("upup", "downdown", "updown", "downup"))
  ymax <- max(df$logP) * 1.4
  bar_colors <- c("#E41A1C", "#4DAF4A", "#FFD700", "#CCCC00")
  
  p <- ggplot(df, aes(x = Pair, y = logP, fill = Pair)) +
    geom_col(color = "black", width = 0.9, show.legend = FALSE) +
    geom_text(aes(label = paste0(OverlapGenes, " genes\np = ", signif(Pvalue, 2))),
              vjust = -0.3, size = 5.2, family = "Helvetica") +
    scale_fill_manual(values = bar_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(y = "-log10(p-value)") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = "Helvetica") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 30, 20, 10)
    )
  
  venn_file <- tempfile(fileext = ".png")
  venn.diagram(
    x = setNames(list(unique(b1$Gene), unique(b2$Gene)),
                 c(bioset1_label, bioset2_label)),
    category.names = c(bioset1_label, bioset2_label),
    filename = venn_file,
    fill = c("#66CCFF", "#66FF66"),
    cex = 1.6,
    fontfamily = "Helvetica",
    cat.cex = 1.6,
    cat.fontfamily = "Helvetica",
    cat.pos = c(-20, 20),
    margin = 0.1,
    scaled = TRUE,
    imagetype = "png",
    resolution = 300,
    height = 2000,
    width = 2000
  )
  
  venn_img <- rasterGrob(png::readPNG(venn_file), interpolate = TRUE)
  pval_text <- paste0("Overlap p-value = ", format(res_df$Pvalue[res_df$Pair == "Total"], scientific = TRUE, digits = 2))
  pval_grob <- ggplot() +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    geom_text(aes(x = 0.5, y = 0.5, label = pval_text),
              size = 6, family = "Helvetica", hjust = 0.5)
  
  venn_with_label <- plot_grid(venn_img, pval_grob, ncol = 1, rel_heights = c(1, 0.15))
  final_plot <- plot_grid(venn_with_label, p, nrow = 1, rel_widths = c(1, 1.5))
  final_plot <- ggdraw(final_plot) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(output_file, final_plot, width = 10, height = 5, dpi = 300,
         units = "in", limitsize = FALSE)
}
