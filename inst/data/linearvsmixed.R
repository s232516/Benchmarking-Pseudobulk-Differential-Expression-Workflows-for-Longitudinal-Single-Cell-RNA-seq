# ============================================================
# Longitudinal toy data: Linear vs Mixed-effects (same plot format)
# - explicit calls: ggplot2::, lme4::, stats::, scales::, RColorBrewer::
# - theme_bw(base_size=16), aspect.ratio=1, text size 18
# - pastel colors (NO colorspace)
# - fewer subjects + NOT parallel lines (random intercept + random slope)
# ============================================================

set.seed(1)

# ----------------------------
# 1) Simulate longitudinal data
# ----------------------------
n_subjects <- 8
timepoints <- c(0, 1, 2, 3, 4)
n_time <- length(timepoints)

dat_long <- data.frame(
  subjectID = factor(rep(seq_len(n_subjects), each = n_time)),
  time = rep(timepoints, times = n_subjects)
)

# Population effects
beta0 <- 5.0
beta_time <- 1.0

# Subject-specific deviations (random intercept + random slope)
sd_intercept <- 1.6
sd_slope     <- 0.45
cor_is       <- -0.25

Sigma <- matrix(
  c(sd_intercept^2, cor_is * sd_intercept * sd_slope,
    cor_is * sd_intercept * sd_slope, sd_slope^2),
  nrow = 2, byrow = TRUE
)

L <- chol(Sigma)
Z <- matrix(stats::rnorm(n_subjects * 2), ncol = 2) %*% L
b0 <- Z[, 1]  # random intercepts
b1 <- Z[, 2]  # random slopes

sd_error <- 0.9
dat_long$y <- (beta0 + b0[dat_long$subjectID]) +
  (beta_time + b1[dat_long$subjectID]) * dat_long$time +
  stats::rnorm(nrow(dat_long), mean = 0, sd = sd_error)

# ----------------------------
# 2) Fit models
# ----------------------------
fit_lm   <- stats::lm(y ~ time, data = dat_long)
fit_lmer <- lme4::lmer(y ~ time + (time | subjectID), data = dat_long)

# ----------------------------
# 3) Straight lines for plotting (abline)
# ----------------------------
# Linear model global line
b_lm <- stats::coef(fit_lm)  # (Intercept), time

# Mixed model population line
b_fix <- lme4::fixef(fit_lmer)  # (Intercept), time

# Mixed model subject-specific intercepts & slopes
re_subj <- lme4::ranef(fit_lmer)[["subjectID"]]
b_re0 <- re_subj[, "(Intercept)"]
b_re1 <- re_subj[, "time"]
subj_levels <- rownames(re_subj)

lines_df <- data.frame(
  subjectID = factor(subj_levels, levels = levels(dat_long$subjectID)),
  intercept = unname(b_fix["(Intercept)"] + b_re0),
  slope     = unname(b_fix["time"] + b_re1)
)

# ----------------------------
# 4) Plot style (same as your QC plots)
# ----------------------------
pastel_blue <- "#A6CEE3"

base_theme <- ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Pastel colors for subjects (no colorspace)
pal <- RColorBrewer::brewer.pal(n_subjects, "Pastel2")
names(pal) <- levels(dat_long$subjectID)

# ----------------------------
# PLOT 1: Linear model (one global trend)
# ----------------------------
p_lm <- ggplot2::ggplot(
  dat_long,
  ggplot2::aes(x = time, y = y)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size  = 3,
    alpha = 0.7
  ) +
  ggplot2::geom_abline(
    intercept = unname(b_lm["(Intercept)"]),
    slope     = unname(b_lm["time"]),
    linewidth = 1.2,
    color     = pastel_blue
  ) +
  ggplot2::labs(
    title = "Linear model",
    x = "Time",
    y = "Outcome"
  ) +
  ggplot2::scale_x_continuous(
    breaks = timepoints
  ) +
  base_theme

# ----------------------------
# PLOT 2: Mixed model (subject-specific trajectories, not parallel)
# ----------------------------
p_mixed <- ggplot2::ggplot(
  dat_long,
  ggplot2::aes(x = time, y = y, color = subjectID)
) +
  ggplot2::geom_point(
    size  = 3,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  # subject-specific straight lines
  ggplot2::geom_abline(
    data = lines_df,
    ggplot2::aes(intercept = intercept, slope = slope, color = subjectID),
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  # population straight line (fixed effects)
  ggplot2::geom_abline(
    intercept = unname(b_fix["(Intercept)"]),
    slope     = unname(b_fix["time"]),
    linewidth = 1.2,
    linetype  = 2,
    color     = "black"
  ) +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::labs(
    title = "Mixed-effects model",
    x = "Time",
    y = "Outcome"
  ) +
  ggplot2::scale_x_continuous(
    breaks = timepoints
  ) +
  base_theme

# Print
p_lm
p_mixed
