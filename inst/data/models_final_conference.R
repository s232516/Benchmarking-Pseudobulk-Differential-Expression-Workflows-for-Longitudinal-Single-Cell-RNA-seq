library(BiocParallel)
library(scater)
library(SummarizedExperiment)
library(MAST)
library(dplyr)
library(edgeR)
library(limma)

# Methods ----------------------------------------------------------------------
# pseudobulk -------------------------------------------------------------------


make_pseudobulk <- function(
    sce,
    ids,
    min_counts = 10,
    verbose = TRUE
) {
  if (is.null(ids)) {
    stop("Please provide ids for pseudobulk aggregation.")
  }

  if (verbose) message("Aggregating cells to pseudobulk...")

  pb <- scater::aggregateAcrossCells(
    sce,
    ids = ids,
    statistics = "sum"
  )

  counts_mat <- SummarizedExperiment::assay(pb, "counts")

  keep_genes <- rowSums(counts_mat) >= min_counts

  if (verbose) {
    message(
      "Keeping ", sum(keep_genes),
      " genes with total pseudobulk counts >= ", min_counts
    )
  }

  pb <- pb[keep_genes, ]

  keep_samples <- colSums(SummarizedExperiment::assay(pb, "counts")) > 0

  if (verbose && any(!keep_samples)) {
    message("Removing ", sum(!keep_samples), " pseudobulk samples with zero total counts.")
  }

  pb <- pb[, keep_samples]

  pb
}

## limma-voom pseudobulk --------------------------------------------------------

run_limma_from_pseudobulk <- function(
    pb,
    design_formula,
    coef_name,
    relevel = NULL,
    verbose = TRUE
) {
  t_start <- proc.time()[["elapsed"]]

  counts <- SummarizedExperiment::assay(pb, "counts")

  keep_samples <- colSums(counts) > 0

  if (any(!keep_samples)) {
    if (verbose) {
      message(
        "Removing ", sum(!keep_samples),
        " pseudobulk samples with zero total counts."
      )
    }

    pb <- pb[, keep_samples]
    counts <- SummarizedExperiment::assay(pb, "counts")
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(pb))
  metadata <- metadata[, colSums(!is.na(metadata)) > 0, drop = FALSE]

  metadata[] <- lapply(metadata, function(x) {
    if (is.character(x)) factor(x) else x
  })

  if (!is.null(relevel)) {
    for (var in names(relevel)) {
      ref <- relevel[[var]]

      if (var %in% colnames(metadata)) {
        metadata[[var]] <- stats::relevel(as.factor(metadata[[var]]), ref = ref)

        if (verbose) {
          message("Relevelled ", var, " with reference = ", ref)
        }
      }
    }
  }

  design <- stats::model.matrix(design_formula, data = metadata)

  if (verbose) {
    message("limma design matrix columns:")
    print(colnames(design))
  }

  coef_idx <- match(coef_name, colnames(design))

  if (is.na(coef_idx)) {
    stop(
      "Could not match coef_name '", coef_name,
      "' to any limma design coefficient. Available coefficients: ",
      paste(colnames(design), collapse = ", ")
    )
  }

  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y)

  v <- limma::voom(
    y,
    design = design,
    plot = FALSE
  )

  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)

  tab <- limma::topTable(
    fit,
    coef = coef_idx,
    number = Inf,
    sort.by = "P"
  )

  df <- as.data.frame(tab)
  df$gene <- rownames(df)

  df <- df |>
    dplyr::transmute(
      gene = gene,
      log2FoldChange = logFC,
      AveExpr = AveExpr,
      stat = t,
      pvalue = P.Value,
      padj = adj.P.Val,
      lfcSE = ifelse(
        is.finite(t) & t != 0,
        abs(logFC / t),
        NA_real_
      )
    )

  df$time <- proc.time()[["elapsed"]] - t_start

  if (verbose) {
    message("limma-voom pseudobulk completed: ", nrow(df), " genes tested.")
  }

  return(df)
}
## edgeR ----------------------------------------------------------------------
run_edger_from_pseudobulk <- function(
    pb,
    design_formula,
    coef_name = NULL,
    contrast = NULL,
    relevel = NULL,
    verbose = TRUE
) {
  t_start <- proc.time()[["elapsed"]]

  counts <- SummarizedExperiment::assay(pb, "counts")

  keep_samples <- colSums(counts) > 0

  if (any(!keep_samples)) {
    if (verbose) {
      message(
        "Removing ", sum(!keep_samples),
        " pseudobulk samples with zero total counts."
      )
    }

    pb <- pb[, keep_samples]
    counts <- SummarizedExperiment::assay(pb, "counts")
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(pb))
  metadata <- metadata[, colSums(!is.na(metadata)) > 0, drop = FALSE]

  metadata[] <- lapply(metadata, function(x) {
    if (is.character(x)) factor(x) else x
  })

  if (!is.null(relevel)) {
    for (var in names(relevel)) {
      ref <- relevel[[var]]

      if (var %in% colnames(metadata)) {
        metadata[[var]] <- stats::relevel(as.factor(metadata[[var]]), ref = ref)

        if (verbose) {
          message("Relevelled ", var, " with reference = ", ref)
        }
      }
    }
  }

  design <- stats::model.matrix(design_formula, data = metadata)

  if (verbose) {
    message("edgeR design matrix columns:")
    print(colnames(design))
  }

  y <- edgeR::DGEList(counts = counts)

  y <- edgeR::calcNormFactors(y)

  y <- edgeR::estimateDisp(y, design)

  fit <- edgeR::glmQLFit(y, design)

  if (!is.null(contrast)) {
    qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  } else {
    if (is.null(coef_name)) {
      stop("Please provide either coef_name or contrast for edgeR.")
    }

    coef_candidates <- if (length(coef_name) > 1) coef_name else c(coef_name)

    matched_coef <- coef_candidates[coef_candidates %in% colnames(design)]

    if (length(matched_coef) == 0) {
      stop(
        "Could not match coef_name to any edgeR design coefficient. Requested: ",
        paste(coef_candidates, collapse = ", "),
        ". Available coefficients: ",
        paste(colnames(design), collapse = ", ")
      )
    }

    if (length(matched_coef) > 1) {
      stop(
        "Multiple coef_name values matched edgeR coefficients: ",
        paste(matched_coef, collapse = ", ")
      )
    }

    coef_idx <- match(matched_coef, colnames(design))

    if (length(coef_idx) == 0) {
      stop(
        "Could not match coef_name '", coef_name,
        "' to any edgeR design coefficient. Available coefficients: ",
        paste(colnames(design), collapse = ", ")
      )
    }

    if (length(coef_idx) > 1) {
      stop(
        "coef_name '", coef_name,
        "' matches multiple edgeR coefficients: ",
        paste(colnames(design)[coef_idx], collapse = ", ")
      )
    }

    qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
  }

  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table

  df <- as.data.frame(tab)
  df$gene <- rownames(df)

  df <- df |>
    dplyr::transmute(
      gene = gene,
      log2FoldChange = logFC,
      logCPM = logCPM,
      stat = F,
      pvalue = PValue,
      padj = FDR,
      lfcSE = dplyr::if_else(
        is.finite(F) & F > 0,
        abs(logFC / sqrt(F)),
        NA_real_
      )
    )

  df$time <- proc.time()[["elapsed"]] - t_start

  if (verbose) {
    message("edgeR pseudobulk completed: ", nrow(df), " genes tested.")
  }

  return(df)
}
## MAST ------------------------------------------------------------------------
run_mast_pseudobulk <- function(
    pb,
    formula = ~ group * time + cngeneson,
    lrt = "groupTreatment:timeEoT",
    norm = TRUE,
    ebayes = TRUE,
    verbose = TRUE
) {
  t_start <- proc.time()[["elapsed"]]

  counts_mat <- as.matrix(SummarizedExperiment::assay(pb, "counts"))

  keep_samples <- colSums(counts_mat) > 0

  if (any(!keep_samples)) {
    if (verbose) {
      message("Removing ", sum(!keep_samples), " pseudobulk samples with zero total counts before MAST.")
    }

    pb <- pb[, keep_samples]
    counts_mat <- as.matrix(SummarizedExperiment::assay(pb, "counts"))
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(pb))

  metadata[] <- lapply(metadata, function(x) {
    if (is.factor(x)) droplevels(x)
    else if (is.character(x)) factor(x)
    else x
  })

  SummarizedExperiment::colData(pb) <- S4Vectors::DataFrame(metadata)

  if (norm) {
    if (verbose) message("Applying library-size log-normalization for MAST...")

    lib_size <- colSums(counts_mat)

    if (any(lib_size <= 0)) {
      stop("Some pseudobulk samples still have zero library size after filtering.")
    }

    expr_mat <- t(t(counts_mat) / lib_size) * 1e6
    expr_mat <- log2(expr_mat + 1)
  } else {
    if (verbose) message("Using log2(counts + 1) for MAST...")
    expr_mat <- log2(counts_mat + 1)
  }

  sca <- MAST::FromMatrix(
    exprsArray = expr_mat,
    cData = metadata
  )

  cdr <- colSums(expr_mat > 0)

  if (stats::sd(cdr) == 0) {
    SummarizedExperiment::colData(sca)$cngeneson <- 0
  } else {
    SummarizedExperiment::colData(sca)$cngeneson <- as.numeric(scale(cdr))
  }

  if (verbose) message("Fitting pseudobulk MAST model...")

  zlmfit <- tryCatch(
    MAST::zlm(
      formula = formula,
      sca = sca,
      method = "glm",
      ebayes = ebayes,
      silent = TRUE
    ),
    error = function(e) {
      message("MAST model failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(zlmfit)) return(NULL)

  s <- tryCatch(
    MAST::summary(zlmfit, doLRT = lrt)$datatable,
    error = function(e) {
      message("MAST summary failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(s)) return(NULL)

  fc <- s |>
    dplyr::filter(contrast == lrt, component == "logFC") |>
    dplyr::transmute(
      gene = primerid,
      log2FoldChange = coef,
      stat = z
    )

  pval <- s |>
    dplyr::filter(contrast == lrt, component == "H") |>
    dplyr::transmute(
      gene = primerid,
      pvalue = `Pr(>Chisq)`
    )

  cont_pcol <- intersect(c("Pr(>|z|)", "Pr(>z)"), names(s))[1]

  cont <- s |>
    dplyr::filter(contrast == lrt, component == "C") |>
    dplyr::transmute(
      gene = primerid,
      betaC = coef,
      zC = z,
      pvalueC = if (!is.na(cont_pcol)) .data[[cont_pcol]] else NA_real_
    )

  res <- Reduce(
    function(x, y) merge(x, y, by = "gene", all = TRUE),
    list(fc, pval, cont)
  )

  res$padj <- p.adjust(res$pvalue, method = "fdr")

  res <- res |>
    dplyr::arrange(padj) |>
    dplyr::mutate(
      lfcSE = ifelse(
        is.finite(stat) & stat != 0,
        abs(log2FoldChange / stat),
        NA_real_
      ),
      betaC_SE = ifelse(
        !is.na(zC) & zC != 0,
        abs(betaC / zC),
        NA_real_
      )
    )

  res$time <- proc.time()[["elapsed"]] - t_start

  if (verbose) {
    message("MAST pseudobulk completed: ", nrow(res), " genes tested.")
  }

  return(res)
}
## DESeq2 pseudobulk -----------------------------------------------------------
run_deseq_from_pseudobulk <- function(
    pb,
    design_formula,
    relevel = NULL,
    sfType = "ratio",
    lrt_term = NULL,
    alpha = 0.05,
    save_dds = FALSE,
    verbose = TRUE
) {
  t_start <- proc.time()[["elapsed"]]

  counts <- SummarizedExperiment::assay(pb, "counts")

  if (any(colSums(counts) <= 0)) {
    keep_samples <- colSums(counts) > 0
    if (verbose) {
      message("Removing ", sum(!keep_samples), " pseudobulk samples with zero total counts.")
    }
    pb <- pb[, keep_samples]
    counts <- SummarizedExperiment::assay(pb, "counts")
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(pb))
  metadata <- metadata[, colSums(!is.na(metadata)) > 0, drop = FALSE]

  metadata[] <- lapply(metadata, function(x) {
    if (is.character(x)) factor(x) else x
  })

  if (!is.null(relevel)) {
    for (var in names(relevel)) {
      ref <- relevel[[var]]

      if (var %in% colnames(metadata)) {
        metadata[[var]] <- stats::relevel(as.factor(metadata[[var]]), ref = ref)

        if (verbose) {
          message("Relevelled ", var, " with reference = ", ref)
        }
      }
    }
  }

  if (verbose) message("Creating DESeq2 dataset from shared pseudobulk object...")

  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = design_formula
    ),
    error = function(e) {
      message("DESeqDataSetFromMatrix failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(dds)) return(NULL)

  if (nrow(dds) < 2 || ncol(dds) < 2) {
    message("Not enough data for DESeq2.")
    return(NULL)
  }

  if (verbose) message("Running DESeq2...")

  dds <- tryCatch(
    DESeq2::DESeq(dds, sfType = sfType, quiet = !verbose),
    error = function(e) {
      message("DESeq() failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(dds)) return(NULL)

  dseq_names <- DESeq2::resultsNames(dds)

  if (verbose) {
    message("\nDESeq2 coefficient names:")
    print(dseq_names)
  }

  if (!is.null(lrt_term)) {
    coef_name <- grep(lrt_term, dseq_names, value = TRUE)

    if (length(coef_name) == 0) {
      stop(
        "Could not match lrt_term '", lrt_term,
        "' to any DESeq2 coefficient."
      )
    }

    if (length(coef_name) > 1) {
      stop(
        "lrt_term '", lrt_term,
        "' maps to multiple DESeq2 coefficients: ",
        paste(coef_name, collapse = ", ")
      )
    }

    res <- DESeq2::results(dds, name = coef_name, alpha = alpha)
  } else {
    res <- DESeq2::results(dds, alpha = alpha)
  }

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[order(df$padj, na.last = TRUE), ]

  df <- df |>
    dplyr::mutate(
      lfcSE = ifelse(
        is.finite(stat) & stat != 0,
        abs(log2FoldChange / stat),
        NA_real_
      )
    )

  df$time <- proc.time()[["elapsed"]] - t_start

  if (save_dds) {
    attr(df, "dds") <- dds
  }

  if (verbose) {
    message("DESeq2 completed: ", nrow(df), " genes tested.")
  }

  df
}
## Run selected methods --------------------------------------------------------
run_methods <- function(
    sce,
    relevel = list(Condition = "Placebo"),
    ids_vars = c("Sample", "Condition"),
    design_formula_mast = ~ Condition + cngeneson,
    design_formula_deseq = ~ Condition,
    design_formula_edger = design_formula_deseq,
    design_formula_limma = design_formula_deseq,
    lrt_term_mast = "ConditionTreatment",
    lrt_term_deseq = "ConditionTreatment",
    lrt_term_edger = lrt_term_deseq,
    lrt_term_limma = lrt_term_deseq,
    run_mast = TRUE,
    run_deseq = TRUE,
    run_edger = FALSE,
    run_limma = FALSE,
    norm = TRUE,
    min_counts = 10,
    ebayes = FALSE,
    verbose = TRUE
) {
  res_list <- list()

  # Prepare metadata before pseudobulk aggregation
  cd <- as.data.frame(SummarizedExperiment::colData(sce))

  if (!is.null(relevel)) {
    for (var in names(relevel)) {
      ref <- relevel[[var]]

      if (var %in% colnames(cd)) {
        cd[[var]] <- stats::relevel(as.factor(cd[[var]]), ref = ref)

        if (verbose) {
          message("Relevelled ", var, " with reference = ", ref)
        }
      }
    }
  }

  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd)

  # Build pseudobulk ids from the same metadata
  if (!all(ids_vars %in% colnames(cd))) {
    stop(
      "Some ids_vars are not present in colData(sce): ",
      paste(setdiff(ids_vars, colnames(cd)), collapse = ", ")
    )
  }

  ids <- SummarizedExperiment::colData(sce)[, ids_vars, drop = FALSE]

  if (verbose) {
    message(
      "Creating shared pseudobulk object using IDs: ",
      paste(ids_vars, collapse = ", ")
    )
  }

  # Create one shared pseudobulk object
  pb <- make_pseudobulk(
    sce = sce,
    ids = ids,
    min_counts = min_counts,
    verbose = verbose
  )

  # MAST pseudobulk
  if (run_mast) {
    if (verbose) message("Running MAST on shared pseudobulk object...")

    mast_res <- run_mast_pseudobulk(
      pb = pb,
      formula = design_formula_mast,
      lrt = lrt_term_mast,
      norm = norm,
      ebayes = ebayes,
      verbose = verbose
    )

    if (!is.null(mast_res)) {
      mast_res$method <- if (norm) {
        "MAST pseudobulk log-normalized"
      } else {
        "MAST pseudobulk log2(counts + 1)"
      }

      mast_res$formula <- deparse(design_formula_mast)
      mast_res$norm <- norm
      mast_res$min_counts <- min_counts

      res_list$mast <- mast_res
    }
  }
  # edgeR pseudobulk
  if (run_edger) {
    if (verbose) message("Running edgeR on shared pseudobulk object...")

    edger_res <- run_edger_from_pseudobulk(
      pb = pb,
      design_formula = design_formula_edger,
      coef_name = lrt_term_edger,
      relevel = NULL,
      verbose = verbose
    )

    if (!is.null(edger_res)) {
      edger_res$method <- "edgeR pseudobulk"
      edger_res$formula <- deparse(design_formula_edger)
      edger_res$norm <- FALSE
      edger_res$min_counts <- min_counts

      res_list$edger <- edger_res
    }
  }

  # DESeq2 pseudobulk
  if (run_deseq) {
    if (verbose) message("Running DESeq2 on shared pseudobulk object...")

    deseq_res <- run_deseq_from_pseudobulk(
      pb = pb,
      design_formula = design_formula_deseq,
      relevel = NULL,
      lrt_term = lrt_term_deseq,
      verbose = verbose
    )

    if (!is.null(deseq_res)) {
      deseq_res$method <- "DESeq2 pseudobulk"
      deseq_res$formula <- deparse(design_formula_deseq)
      deseq_res$norm <- FALSE
      deseq_res$min_counts <- min_counts

      res_list$deseq <- deseq_res
    }
  }
  # limma-voom pseudobulk
  if (run_limma) {
    if (verbose) message("Running limma-voom on shared pseudobulk object...")

    limma_res <- run_limma_from_pseudobulk(
      pb = pb,
      design_formula = design_formula_limma,
      coef_name = lrt_term_limma,
      relevel = NULL,
      verbose = verbose
    )

    if (!is.null(limma_res)) {
      limma_res$method <- "limma-voom pseudobulk"
      limma_res$formula <- deparse(design_formula_limma)
      limma_res$norm <- FALSE
      limma_res$min_counts <- min_counts

      res_list$limma <- limma_res
    }
  }

  if (verbose) message("Finished running selected methods.")

  return(res_list)
}


# Pipelines --------------------------------------------------------------------

pipeline_rescue_and_methods <- function(
    sce,
    run_mast = TRUE,
    run_deseq = TRUE,
    run_edger = FALSE,
    run_limma = FALSE,
    design_formula_mast = ~ group * time + cngeneson,
    design_formula_deseq = ~ group * time,
    design_formula_edger = design_formula_deseq,
    design_formula_limma = design_formula_deseq,
    lrt_term_mast = "groupTreatment:timeEoT",
    lrt_term_deseq = "groupTreatment.timeEoT",
    lrt_term_edger = "timeEoT:groupTreatment",
    lrt_term_limma = "timeEoT:groupTreatment",
    ids_vars = c("subjectID", "time"),
    relevel = list(group = "Placebo", time = "EoT"),
    norm = TRUE,
    min_counts = 10,
    verbose = TRUE
) {
  run_methods(
    sce = sce,
    relevel = relevel,
    ids_vars = ids_vars,
    design_formula_mast = design_formula_mast,
    design_formula_deseq = design_formula_deseq,
    design_formula_edger = design_formula_edger,
    design_formula_limma = design_formula_limma,
    lrt_term_mast = lrt_term_mast,
    lrt_term_deseq = lrt_term_deseq,
    lrt_term_edger = lrt_term_edger,
    lrt_term_limma = lrt_term_limma,
    run_mast = run_mast,
    run_deseq = run_deseq,
    run_edger = run_edger,
    run_limma = run_limma,
    norm = norm,
    min_counts = min_counts,
    verbose = verbose
  )
}

# Run Simulations  -------------------------------------------------------------
# Read pins with simulations
simulation <- pins::pin_read(board, "simulations_final")

sim_rescue      <- simulation$rescueSim
sim_rescue_123  <- simulation$rescueSim_123
sim_rescue_111  <- simulation$rescueSim_111
sim_rescue_222  <- simulation$rescueSim_222
sim_rescue_333  <- simulation$rescueSim_333


sim_rescue_norm      <- simulation[["rescueSim norm"]]
sim_rescue_123_norm  <- simulation[["rescueSim_123 norm"]]
sim_rescue_111_norm  <- simulation[["rescueSim_111 norm"]]
sim_rescue_222_norm  <- simulation[["rescueSim_222 norm"]]
sim_rescue_333_norm  <- simulation[["rescueSim_333 norm"]]




# ------------------------------------------------------------------------------
# RescueSim pseudobulk runs
# ------------------------------------------------------------------------------

prepare_rescue_sce <- function(
    sce,
    group_ref = "Placebo",
    time_ref = "EoT"
) {
  sce_ref <- sce

  cd <- as.data.frame(SummarizedExperiment::colData(sce_ref))

  cd$group <- stats::relevel(as.factor(cd$group), ref = group_ref)
  cd$time <- stats::relevel(as.factor(cd$time), ref = time_ref)

  SummarizedExperiment::colData(sce_ref) <- S4Vectors::DataFrame(cd)

  sce_ref
}


run_rescue_methods <- function(
    sce,
    analysis_name,
    ids_vars = c("subjectID", "time"),
    min_counts = 10,
    verbose = TRUE
) {
  sce_ref <- prepare_rescue_sce(
    sce = sce,
    group_ref = "Placebo",
    time_ref = "EoT"
  )

  out <- list()

  # MAST with CDR, normalized
  out$mast_cdr <- pipeline_rescue_and_methods(
    sce = sce_ref,
    run_mast = TRUE,
    run_deseq = FALSE,
    design_formula_mast = ~ time + group:time + subjectID + cngeneson,
    lrt_term_mast = "timeEoT:groupTreatment",
    ids_vars = ids_vars,
    relevel = list(group = "Placebo", time = "EoT"),
    norm = TRUE,
    min_counts = min_counts,
    verbose = verbose
  )

  out$limma_s <- pipeline_rescue_and_methods(
    sce = sce_ref,
    run_mast = FALSE,
    run_deseq = FALSE,
    run_edger = FALSE,
    run_limma = TRUE,
    design_formula_limma = ~ time + group:time,
    lrt_term_limma = "timeEoT:groupTreatment",
    ids_vars = ids_vars,
    relevel = list(group = "Placebo", time = "EoT"),
    min_counts = min_counts,
    verbose = verbose
  )

  # MAST with CDR, non-normalized
  out$mast_cdr_nn <- pipeline_rescue_and_methods(
    sce = sce_ref,
    run_mast = TRUE,
    run_deseq = FALSE,
    design_formula_mast = ~ time + group:time + subjectID + cngeneson,
    lrt_term_mast = "timeEoT:groupTreatment",
    ids_vars = ids_vars,
    relevel = list(group = "Placebo", time = "EoT"),
    norm = FALSE,
    min_counts = min_counts,
    verbose = verbose
  )

  # DESeq2 standard old formula
  out$deseq_s <- pipeline_rescue_and_methods(
    sce = sce_ref,
    run_mast = FALSE,
    run_deseq = TRUE,
    design_formula_deseq = ~ time + group:time,
    lrt_term_deseq = "timeEoT.groupTreatment",
    ids_vars = ids_vars,
    relevel = list(group = "Placebo", time = "EoT"),
    min_counts = min_counts,
    verbose = verbose
  )

  out$edger_s <- pipeline_rescue_and_methods(
    sce = sce_ref,
    run_mast = FALSE,
    run_deseq = FALSE,
    run_edger = TRUE,
    design_formula_edger = ~ time + group:time,
    lrt_term_edger = "timeEoT:groupTreatment",
    ids_vars = ids_vars,
    relevel = list(group = "Placebo", time = "EoT"),
    min_counts = min_counts,
    verbose = verbose
  )

  # Add labels
  for (variant in names(out)) {
    for (method in names(out[[variant]])) {
      out[[variant]][[method]]$analysis <- analysis_name
      out[[variant]][[method]]$analysis_variant <- variant
    }
  }

  out
}

rescue_sims <- list(
  "rescueSim" = sim_rescue,
  "rescueSim (123)" = sim_rescue_123,
  "rescueSim (111)" = sim_rescue_111,
  "rescueSim (222)" = sim_rescue_222,
  "rescueSim (333)" = sim_rescue_333,
  "rescueSim (norm)" = sim_rescue_norm,
  "rescueSim (123) (norm)" = sim_rescue_123_norm,
  "rescueSim (111) (norm)" = sim_rescue_111_norm,
  "rescueSim (222) (norm)" = sim_rescue_222_norm,
  "rescueSim (333) (norm)" = sim_rescue_333_norm
)

rescue_results <- lapply(
  names(rescue_sims),
  function(nm) {
    run_rescue_methods(
      sce = rescue_sims[[nm]],
      analysis_name = nm,
      ids_vars = c("subjectID", "time"),
      min_counts = 10,
      verbose = TRUE
    )
  }
)

names(rescue_results) <- names(rescue_sims)

flatten_method_results <- function(results_list) {
  flat <- list()

  for (analysis_name in names(results_list)) {
    analysis_results <- results_list[[analysis_name]]

    for (variant_name in names(analysis_results)) {
      variant_results <- analysis_results[[variant_name]]

      for (method_name in names(variant_results)) {
        key <- paste(analysis_name, variant_name, method_name, sep = "__")
        flat[[key]] <- variant_results[[method_name]]
      }
    }
  }

  flat
}

rescue_method_tables <- flatten_method_results(rescue_results)
method_data <- data.table::rbindlist(
  rescue_method_tables,
  fill = TRUE
)
# Pin data ---------------------------------------------------------------------

board <- pins::board_local("nrgn-data")

pins::pin_write(
  method_data,
  name = "method_data_paper",
  type = "rds",
  board = board,
  versioned = TRUE
)

# Merge ground_truth with method data ------------------------------------------

method_data <- pins::pin_read(board, "method_data_paper")
ground_truth <- pins::pin_read(board, "ground_truth_final")

make_analysis_join <- function(x) {
  x |>
    gsub(
      pattern = "^rescueSim \\((\\d+)\\)$",
      replacement = "rescueSim_\\1",
      x = _
    ) |>
    gsub(
      pattern = "^rescueSim \\(norm\\)$",
      replacement = "rescueSim norm",
      x = _
    ) |>
    gsub(
      pattern = "^rescueSim \\((\\d+)\\) \\(norm\\)$",
      replacement = "rescueSim_\\1 norm",
      x = _
    )
}

method_data <- method_data |>
  dplyr::mutate(
    analysis_join = make_analysis_join(analysis)
  )

rescue_truth_names <- c(
  "rescueSim",
  "rescueSim_123",
  "rescueSim_111",
  "rescueSim_222",
  "rescueSim_333",
  "rescueSim norm",
  "rescueSim_123 norm",
  "rescueSim_111 norm",
  "rescueSim_222 norm",
  "rescueSim_333 norm"
)

compare_data <- dplyr::bind_rows(
  lapply(rescue_truth_names, function(sim) {
    dplyr::left_join(
      method_data |> dplyr::filter(analysis_join == sim),
      ground_truth[[sim]],
      by = "gene"
    )
  })
)

pins::pin_write(
  compare_data,
  name = "compare_data_paper",
  type = "rds",
  board = board,
  versioned = TRUE
)

# CI coverage summary ----------------------------------------------------------

compare_data <- pins::pin_read(board, "compare_data_paper")

summarise_ci_coverage <- function(df, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)

  df <- df |>
    dplyr::filter(!is.na(trueDE), !is.na(true_logFC)) |>
    dplyr::mutate(
      lower = log2FoldChange - z * lfcSE,
      upper = log2FoldChange + z * lfcSE,
      in_ci = true_logFC >= lower & true_logFC <= upper,
      regulation = dplyr::case_when(
        true_logFC > 0 ~ "up",
        true_logFC < 0 ~ "down",
        TRUE ~ "not_sig"
      )
    )

  dplyr::tibble(
    total_genes = nrow(df),

    in_ci_total = sum(df$in_ci, na.rm = TRUE),

    total_up = sum(df$regulation == "up", na.rm = TRUE),
    in_ci_up = sum(df$in_ci & df$regulation == "up", na.rm = TRUE),

    total_down = sum(df$regulation == "down", na.rm = TRUE),
    in_ci_down = sum(df$in_ci & df$regulation == "down", na.rm = TRUE),

    total_not_sig = sum(df$regulation == "not_sig", na.rm = TRUE),
    in_ci_not_sig = sum(df$in_ci & df$regulation == "not_sig", na.rm = TRUE)
  )
}
compute_ci_summary <- function(
    compare_data,
    analyses,
    analysis_key = "analysis_join",
    ci_level = 0.95
) {
  compare_data |>
    dplyr::filter(.data[[analysis_key]] %in% analyses) |>
    dplyr::group_by(.data[[analysis_key]], method, formula) |>
    dplyr::group_modify(~ summarise_ci_coverage(.x, conf_level = ci_level)) |>
    dplyr::ungroup() |>
    dplyr::rename(analysis = .data[[analysis_key]]) |>
    dplyr::mutate(
      conf_level = ci_level,
      prop_ci_total = in_ci_total / total_genes,
      prop_ci_up = dplyr::if_else(total_up > 0, in_ci_up / total_up, NA_real_),
      prop_ci_down = dplyr::if_else(total_down > 0, in_ci_down / total_down, NA_real_),
      prop_ci_not_sig = dplyr::if_else(
        total_not_sig > 0,
        in_ci_not_sig / total_not_sig,
        NA_real_
      )
    )
}

ci_level <- 0.95

rescue_analyses <- c(
  "rescueSim",
  "rescueSim_123",
  "rescueSim_111",
  "rescueSim_222",
  "rescueSim_333"
)

rescue_norm_analyses <- c(
  "rescueSim norm",
  "rescueSim_123 norm",
  "rescueSim_111 norm",
  "rescueSim_222 norm",
  "rescueSim_333 norm"
)

ci_data <- dplyr::bind_rows(
  compute_ci_summary(
    compare_data = compare_data,
    analyses = rescue_analyses,
    analysis_key = "analysis_join",
    ci_level = ci_level
  ) |>
    dplyr::mutate(sim = "rescueSim"),

  compute_ci_summary(
    compare_data = compare_data,
    analyses = rescue_norm_analyses,
    analysis_key = "analysis_join",
    ci_level = ci_level
  ) |>
    dplyr::mutate(sim = "rescueSim norm")
)

pins::pin_write(
  ci_data,
  name = "rescueSim_ci_paper",
  type = "rds",
  board = board,
  versioned = TRUE
)
