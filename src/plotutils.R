get_default_theme <- function() {
  My_Theme = theme(
    text = element_text(size = 6),
    line = element_line(size = 0.5 * 5 / 14),
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.key.size = unit(0.1, "in"),
    axis.ticks = element_line(size = 0.5 * 5 / 14),
    axis.ticks.length = unit(1, 'pt'),
    axis.line  = element_line(size = 0.5 * 5 / 14)
    )
  return(My_Theme)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[n]
}
