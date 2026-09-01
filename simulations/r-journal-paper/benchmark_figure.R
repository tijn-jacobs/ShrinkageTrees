################################################################################
# Figure for the timing benchmark.
#
# Sourced by benchmark_timings.R. Kept separate so the figure can also be
# rebuilt from outputs/benchmark_timings.rds without refitting:
#   source("benchmark_figure.R")
#   raw <- readRDS("outputs/benchmark_timings.rds")
#   ggplot2::ggsave("outputs/benchmark_scaling.pdf",
#                   benchmark_figure(benchmark_summary(raw)),
#                   width = 6, height = 3.5)
################################################################################

fn_display <- c(SurvivalBART = "SurvivalBART()", SurvivalDART = "SurvivalDART()",
                HorseTrees = "HorseTrees()", `BART::mc.abart` = "abart()")

# raw: long-format timings, one row per fit (fn, n, seconds).
benchmark_summary <- function(raw) {
  ok <- raw[!is.na(raw$seconds), ]
  agg_mean <- aggregate(seconds ~ fn + n, data = ok, FUN = mean)
  agg_sd   <- aggregate(seconds ~ fn + n, data = ok, FUN = sd)
  names(agg_mean)[3] <- "mean"
  names(agg_sd)[3]   <- "sd"
  agg <- merge(agg_mean, agg_sd, by = c("fn", "n"))
  agg$sd[is.na(agg$sd)] <- 0
  agg <- agg[order(agg$fn, agg$n), ]
  agg$label <- factor(fn_display[as.character(agg$fn)],
                      levels = unname(fn_display))
  agg
}

benchmark_figure <- function(agg, breaks = seq(0, 15, by = 5)) {
  ggplot2::ggplot(
    agg, ggplot2::aes(x = n, y = mean, colour = label, fill = label)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_y_continuous(breaks = breaks, limits = c(0, NA)) +
    ggplot2::labs(x = "Sample size", y = "Wall-clock seconds",
                  colour = NULL, fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
