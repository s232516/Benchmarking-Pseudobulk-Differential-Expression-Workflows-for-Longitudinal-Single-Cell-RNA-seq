library(SingleCellExperiment)

board <- pins::board_local("nrgn-data")
simulation <- pins::pin_read(board, "simulations_final")
splatpop <- simulation$splatpop
rescueSim <- simulation$rescueSim
rescueSim_norm <- simulation[["rescueSim norm"]]

# ---------- splatPop ----------
sce <- splatpop

genes_splat <- nrow(sce)
cells_splat <- ncol(sce)

# subject/sample inferred from colnames: "sample_01:Cell1"
subject_splat <- sub(":.*", "", colnames(sce))
n_subjects_splat <- length(unique(subject_splat))
cells_per_subject_splat <- table(subject_splat)

# condition balance (cells)
cells_by_condition_splat <- table(colData(sce)$Condition)

# condition per subject (should be constant within subject)
cond_by_subject_splat <- tapply(as.character(colData(sce)$Condition),
                                subject_splat, function(x) unique(x))

list(
  genes = genes_splat,
  cells = cells_splat,
  subjects = n_subjects_splat,
  cells_per_subject = summary(as.numeric(cells_per_subject_splat)),
  cells_by_condition = cells_by_condition_splat
)

# ---------- rescueSim ----------
sce <- rescueSim

genes_rescue <- nrow(sce)
cells_rescue <- ncol(sce)

n_subjects_rescue <- length(unique(colData(sce)$subjectID))
n_samples_rescue  <- length(unique(colData(sce)$sampleID))
cells_per_sample_rescue <- table(colData(sce)$sampleID)
cells_per_subject_rescue <- table(colData(sce)$subjectID)

list(
  genes = genes_rescue,
  cells = cells_rescue,
  subjects = n_subjects_rescue,
  samples = n_samples_rescue,
  cells_per_sample = summary(as.numeric(cells_per_sample_rescue)),
  cells_per_subject = summary(as.numeric(cells_per_subject_rescue))
)



qc_cell <- function(sce, assay_name = "counts") {
  x <- assay(sce, assay_name)
  libsize <- colSums(x)
  detected_genes <- colSums(x > 0)
  pct_zeros_cell <- colMeans(x == 0) * 100
  data.frame(libsize, detected_genes, pct_zeros_cell)
}

qc_splat  <- qc_cell(splatpop)
qc_rescue <- qc_cell(rescueSim)


summary(qc_splat$libsize)
summary(qc_rescue$libsize)

# splatPop has only 50 genes, so use proportions too:
qc_splat$cdr <- qc_splat$detected_genes / nrow(splatpop)  # 0..1
summary(qc_splat$cdr)

cor(qc_rescue$libsize, qc_rescue$detected_genes, method="spearman")
cor(qc_splat$libsize, qc_splat$detected_genes, method="spearman")


gene_mean_var <- function(sce, assay_name = "counts") {
  x <- assay(sce, assay_name)
  data.frame(
    mean = rowMeans(x),
    var  = apply(x, 1, var),
    zero_frac = rowMeans(x == 0)
  )
}

gv_splat  <- gene_mean_var(splatpop)
gv_rescue <- gene_mean_var(rescueSim)


par(mfrow=c(2,2))
# ---------- rescueSim norm ----------
sce <- rescueSim_norm

genes_rescue_norm <- nrow(sce)
cells_rescue_norm <- ncol(sce)

n_subjects_rescue_norm <- length(unique(colData(sce)$subjectID))
n_samples_rescue_norm  <- length(unique(colData(sce)$sampleID))
cells_per_sample_rescue_norm <- table(colData(sce)$sampleID)
cells_per_subject_rescue_norm <- table(colData(sce)$subjectID)

list(
  genes = genes_rescue_norm,
  cells = cells_rescue_norm,
  subjects = n_subjects_rescue_norm,
  samples = n_samples_rescue_norm,
  cells_per_sample = summary(as.numeric(cells_per_sample_rescue_norm)),
  cells_per_subject = summary(as.numeric(cells_per_subject_rescue_norm))
)



qc_cell <- function(sce, assay_name = "counts") {
  x <- assay(sce, assay_name)
  libsize <- colSums(x)
  detected_genes <- colSums(x > 0)
  pct_zeros_cell <- colMeans(x == 0) * 100
  data.frame(libsize, detected_genes, pct_zeros_cell)
}

qc_splat  <- qc_cell(splatpop)
qc_rescue <- qc_cell(rescueSim)
qc_rescue_norm <- qc_cell(rescueSim_norm)

summary(qc_splat$libsize)
summary(qc_rescue$libsize)
summary(qc_rescue_norm$libsize)


qc_splat$cdr <- qc_splat$detected_genes / nrow(splatpop)
summary(qc_splat$cdr)

qc_rescue$cdr <- qc_rescue$detected_genes / nrow(rescueSim)
summary(qc_rescue$cdr)


qc_rescue_norm$cdr <- qc_rescue_norm$detected_genes / nrow(rescueSim_norm)
summary(qc_rescue_norm$cdr)


cor(qc_splat$libsize, qc_splat$detected_genes, method="spearman")
cor(qc_rescue$libsize, qc_rescue$detected_genes,
    method="spearman")
cor(qc_rescue_norm$libsize, qc_rescue_norm$detected_genes,
    method="spearman")


gene_mean_var <- function(sce, assay_name = "counts") {
  x <- assay(sce, assay_name)
  data.frame(
    mean = rowMeans(x),
    var  = apply(x, 1, var),
    zero_frac = rowMeans(x == 0)
  )
}

gv_splat  <- gene_mean_var(splatpop)
gv_rescue <- gene_mean_var(rescueSim)
gv_rescue_norm <- gene_mean_var(rescueSim_norm)


par(mfrow=c(2,2))
# splatPop -------------------------------------------------------------------
pastel_blue <- "gray"

# Histogram: library sizes (volcano-style)
p_lib_hist_splat <- ggplot2::ggplot(
  qc_splat,
  ggplot2::aes(x = libsize)
) +
  ggplot2::geom_histogram(
    color = "black",
    fill  = pastel_blue,
    alpha = 0.7,
    bins  = 40
  ) +
  ggplot2::labs(
    title = "Histogram of UMIs per cell",
    x = "UMIs/cell",
    y = "Count"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Detected genes vs UMIs
p_detect_vs_umi_splat <- ggplot2::ggplot(
  qc_splat,
  ggplot2::aes(x = libsize, y = detected_genes)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::labs(
    title = "Detected genes vs UMIs",
    x = "UMIs/cell",
    y = "Detected genes"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Mean–variance (log-log)
p_mean_var_splat <- ggplot2::ggplot(
  gv_splat,
  ggplot2::aes(x = mean, y = var)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Mean–variance trend",
    x = "Gene mean",
    y = "Gene variance"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Zeros vs mean (log-x)
p_zero_vs_mean_splat <- ggplot2::ggplot(
  gv_splat,
  ggplot2::aes(x = mean, y = zero_frac)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Zero fraction vs mean",
    x = "Gene mean",
    y = "Zero fraction"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )


p_lib_hist_splat
p_detect_vs_umi_splat
p_mean_var_splat
p_zero_vs_mean_splat
# Rescue----------------------------------------------------------------------

pastel_blue <- "gray"

#  Histogram: library sizes (volcano-style)
p_lib_hist <- ggplot2::ggplot(
  qc_rescue,
  ggplot2::aes(x = libsize)
) +
  ggplot2::geom_histogram(
    color = "black",
    fill  = pastel_blue,
    alpha = 0.7,
    bins  = 40
  ) +
  ggplot2::labs(
    title = "Histogram of UMIs per cell",
    x = "UMIs/cell",
    y = "Count"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Detected genes vs UMIs
p_detect_vs_umi <- ggplot2::ggplot(
  qc_rescue,
  ggplot2::aes(x = libsize, y = detected_genes)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::labs(
    title = "Detected genes vs UMIs",
    x = "UMIs/cell",
    y = "Detected genes"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Mean–variance (log-log)
p_mean_var <- ggplot2::ggplot(
  gv_rescue,
  ggplot2::aes(x = mean, y = var)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Mean–variance trend",
    x = "Gene mean",
    y = "Gene variance"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Zeros vs mean (log-x)
p_zero_vs_mean <- ggplot2::ggplot(
  gv_rescue,
  ggplot2::aes(x = mean, y = zero_frac)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Zero fraction vs mean",
    x = "Gene mean",
    y = "Zero fraction"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )


p_lib_hist
p_detect_vs_umi
p_mean_var
p_zero_vs_mean


# Rescue norm-------------------------------------------------------------------

pastel_blue <- "gray"

#  Histogram: library sizes (volcano-style)
p_lib_hist_norm <- ggplot2::ggplot(
  qc_rescue_norm,
  ggplot2::aes(x = libsize)
) +
  ggplot2::geom_histogram(
    color = "black",
    fill  = pastel_blue,
    alpha = 0.7,
    bins  = 40
  ) +
  ggplot2::labs(
    title = "Histogram of UMIs per cell",
    x = "UMIs/cell",
    y = "Count"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Detected genes vs UMIs
p_detect_vs_umi_norm <- ggplot2::ggplot(
  qc_rescue_norm,
  ggplot2::aes(x = libsize, y = detected_genes)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::labs(
    title = "Detected genes vs UMIs",
    x = "UMIs/cell",
    y = "Detected genes"
  ) +
  ggplot2::scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Mean–variance (log-log)
p_mean_var_norm <- ggplot2::ggplot(
  gv_rescue_norm,
  ggplot2::aes(x = mean, y = var)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::scale_y_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Mean–variance trend",
    x = "Gene mean",
    y = "Gene variance"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )

# Zeros vs mean (log-x)
p_zero_vs_mean_norm <- ggplot2::ggplot(
  gv_rescue_norm,
  ggplot2::aes(x = mean, y = zero_frac)
) +
  ggplot2::geom_point(
    color = pastel_blue,
    size = 3,
    alpha = 0.7
  ) +
  ggplot2::scale_x_log10(breaks = scales::breaks_log(n = 6)) +
  ggplot2::labs(
    title = "Zero fraction vs mean",
    x = "Gene mean",
    y = "Zero fraction"
  ) +
  ggplot2::theme_bw(base_size = 16) +
  ggplot2::theme(
    aspect.ratio = 1,
    text = ggplot2::element_text(size = 18)
  )


p_lib_hist_norm
p_detect_vs_umi_norm
p_mean_var_norm
p_zero_vs_mean_norm



subject <- sub(":.*", "", colnames(splatpop))
stopifnot(all(table(subject) == 100))   # expect 100 cells/subject

# One condition per subject?
table(sapply(split(as.character(colData(splatpop)$Condition), subject),
             function(x) length(unique(x))))

# Samples per subject (should be 2 if baseline+EoT)
table(table(colData(rescueSim)$sampleID[ , drop=TRUE]))  # (not always meaningful)
table(table(colData(rescueSim)$subjectID))               # cells per subject

# Better: number of samples per subject
samples_per_subject <- tapply(as.character(colData(rescueSim)$sampleID),
                              colData(rescueSim)$subjectID,
                              function(x) length(unique(x)))
table(samples_per_subject)


# Samples per subject (should be 2 if baseline+EoT)
table(table(colData(rescueSim_norm)$sampleID[ , drop=TRUE]))  # (not always meaningful)
table(table(colData(rescueSim_norm)$subjectID))               # cells per subject

# Better: number of samples per subject
samples_per_subject_norm <- tapply(as.character(colData(rescueSim_norm)$sampleID),
                              colData(rescueSim_norm)$subjectID,
                              function(x) length(unique(x)))
table(samples_per_subject)




