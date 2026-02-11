# ============================================================
# Linear vs Mixed effects — 2 plots (same “volcano module” style)
# - explicit calls: ggplot2::, lme4::, stats::, scales::
# - theme_bw(base_size=16), aspect.ratio=1, pretty breaks
# - pastel palette (points + group lines)
# ============================================================

# If not already loaded in your script:
# library(lme4)

# ----------------------------
# 1) Simulate grouped data
# ----------------------------
set.seed(1)

n_groups <- 12
n_per_group <- 40
group <- factor(rep(seq_len(n_groups), each = n_per_group))
N <- n_groups * n_per_group

x <- stats::runif(N, 0, 10)

# True model: common slope + group-specific intercept
beta0 <- 2
beta1 <- 0.8
sd_group_intercept <- 2.2
sd_error <- 1.2

b0 <- stats::rnorm(n_groups, mean = 0, sd = sd_group_intercept)
y  <- beta0 + beta1 * x + b0[group] + stats::rnorm(N, 0, sd_error)

dat <- data.frame(group = group, x = x, y = y)

# ----------------------------
# 2) Fit models
# ----------------------------
fit_lm   <- stats::lm(y ~ x, data = dat)
fit_lmer <- lme4::lmer(y ~ x + (1 | group), data = dat)

# ----------------------------
# 3) Prediction grids
# ----------------------------
grid <- data.frame(x = seq(min(dat$x), max(dat$x), length.out = 200))
grid$yhat_lm <- stats::predict(fit_lm, newdata = grid)

grid_by_group <- base::expand.grid(
  x = seq(min(dat$x), max(dat$x), length.out = 200),
  group = levels(dat$group)
)
# per-group fitted lines (includes random effects)
grid_by_group$yhat_mixed_group <- stats::predict(fit_lmer, newdata = grid_by_group, re.form = NULL)

# population line only (fixed effects)
grid_pop <- data.frame(x = seq(min(dat$x), max(dat$x), length.out = 200))
grid_pop$yhat_mixed_pop <- stats::predict(fit_lmer, newdata = grid_pop, re.form = NA)

# ----------------------------
# 4) Volcano-style theme bits
# ----------------------------
pastel_blue <- "#A6CEE3"  # your previous pastel blue

base_theme <- ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# ----------------------------
# Pastel palette for groups
# ----------------------------
# Uses RColorBrewer pastel palettes (very common + readable).
# If you don't have it installed: install.packages("RColorBrewer")
pastel_groups <- RColorBrewer::brewer.pal(n = n_groups, name = "Set3")

# ----------------------------
# PLOT 1: Linear model (one global line) — pastel
# ----------------------------
p_linear <- ggplot2::ggplot(
  dat,
  ggplot2::aes(x = x, y = y)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size  = 3,
    alpha = 0.7
  ) +
  ggplot2::geom_line(
    data = grid,
    ggplot2::aes(x = x, y = yhat_lm),
    color = pastel_blue,
    linewidth = 1.2
  ) +
  ggplot2::labs(
    title = "Fixed-effects model",
    x = "x",
    y = "y"
  ) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
  base_theme

# ----------------------------
# PLOT 2: Mixed model (group lines + population line) — pastel groups
# ----------------------------
p_mixed <- ggplot2::ggplot(
  dat,
  ggplot2::aes(x = x, y = y, color = group)
) +
  ggplot2::geom_point(
    size  = 3,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(
    data = grid_by_group,
    ggplot2::aes(x = x, y = yhat_mixed_group, group = group),
    linewidth = 0.9,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(
    data = grid_pop,
    ggplot2::aes(x = x, y = yhat_mixed_pop),
    color = pastel_blue,
    linewidth = 1.2,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(values = pastel_groups) +
  ggplot2::labs(
    title = "Mixed-effects model",
    x = "x",
    y = "y"
  ) +
  ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
  base_theme

# ----------------------------
# Print
# ----------------------------
p_linear
p_mixed
