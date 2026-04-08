###############################################################################
# theme_nature.R — Nature Publishing Group–style ggplot2 theme
# Usage: source(here::here("scripts/utils/theme_nature.R"))
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
})

#' Nature-style ggplot2 theme
#'
#' 7 pt font, no gridlines, minimal axes. Follows Nature figure guidelines.
#' @param base_size Base font size in points (default 7).
#' @param base_family Font family (default "Helvetica" / "Arial").
theme_nature <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Text
      text              = element_text(size = base_size, colour = "black"),
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0, margin = margin(b = 4)),
      plot.subtitle     = element_text(size = base_size, hjust = 0,
                                       margin = margin(b = 3)),
      # Axes
      axis.title        = element_text(size = base_size, face = "bold"),
      axis.text         = element_text(size = base_size - 0.5, colour = "black"),
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "pt"),
      # Legend
      legend.title      = element_text(size = base_size, face = "bold"),
      legend.text       = element_text(size = base_size - 0.5),
      legend.key.size   = unit(3, "mm"),
      legend.background = element_blank(),
      legend.spacing    = unit(1, "mm"),
      # Strips (for facets)
      strip.text        = element_text(size = base_size, face = "bold"),
      strip.background  = element_blank(),
      # Panel
      panel.grid        = element_blank(),
      panel.border      = element_blank(),
      panel.background  = element_blank(),
      # Margins
      plot.margin       = margin(5, 5, 5, 5, "pt")
    )
}

#' Save a plot in Nature-compliant formats (PDF vector + PNG raster)
#' @param p A ggplot object.
#' @param filename Base filename without extension.
#' @param width Width in inches (default: single column = 3.5).
#' @param height Height in inches.
#' @param dir Output directory (default: results/figures).
save_nature_fig <- function(p, filename, width = 3.5, height = 3.5,
                            dir = here::here("results", "figures")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(filename, ".pdf")), plot = p,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(dir, paste0(filename, ".png")), plot = p,
         width = width, height = height, units = "in", dpi = 300)
  message("✔ Saved: ", filename, " (.pdf + .png)")
}

#' Add panel label (a, b, c…) to a ggplot — for multi-panel figures
#' @param p A ggplot object.
#' @param label Character label (e.g., "a").
add_panel_label <- function(p, label) {
  p + labs(tag = label) +
    theme(plot.tag = element_text(size = 10, face = "bold"))
}

message("✔ theme_nature.R loaded")
