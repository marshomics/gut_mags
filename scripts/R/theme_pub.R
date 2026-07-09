# theme_pub.R -- shared ggplot theme + dual-format save (PNG + editable-text SVG)
# Source this from every R figure script. svglite writes <text> elements, so all
# labels remain selectable/editable in vector editors.

suppressPackageStartupMessages({ library(ggplot2); library(svglite) })

theme_pub <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.key.size = unit(3, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

niche_colors <- function(cfg) unlist(cfg$figures$palette)

# Save a ggplot (or grid object) to png + svg. width/height in mm.
save_pub <- function(plot, path_no_ext, cfg, width = 120, height = 90) {
  dpi <- cfg$figures$dpi
  fmts <- cfg$figures$formats
  if ("png" %in% fmts)
    ggsave(paste0(path_no_ext, ".png"), plot, width = width, height = height,
           units = "mm", dpi = dpi)
  if ("svg" %in% fmts)
    ggsave(paste0(path_no_ext, ".svg"), plot, width = width, height = height,
           units = "mm", device = svglite::svglite)
}

# Save a BASE-graphics plotting expression to png + svg (for phytools simmap).
save_base <- function(expr, path_no_ext, cfg, width = 120, height = 120) {
  dpi <- cfg$figures$dpi
  inch_w <- width / 25.4; inch_h <- height / 25.4
  if ("png" %in% cfg$figures$formats) {
    png(paste0(path_no_ext, ".png"), width = inch_w, height = inch_h,
        units = "in", res = dpi); force(expr); dev.off()
  }
  if ("svg" %in% cfg$figures$formats) {
    svglite::svglite(paste0(path_no_ext, ".svg"), width = inch_w, height = inch_h)
    force(expr); dev.off()
  }
}
