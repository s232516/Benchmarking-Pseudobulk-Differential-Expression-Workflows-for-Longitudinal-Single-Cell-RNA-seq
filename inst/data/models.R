library(BiocParallel)
# Methods ----------------------------------------------------------------------
## MAST ------------------------------------------------------------------------
run_mast <- function(
    sce,
    formula = ~Condition,
    lrt = "ConditionTreatment",
    random_var = NULL,
    subject_var = NULL,
    random = FALSE, # switch between fixed and random effects
    method = NULL, # override if needed
    ebayes = TRUE,
    check_sanity = FALSE,
    verbose = TRUE # print progress messages
) {
  # Convert SCE to SingleCellAssay
  if (verbose) message("Converting SCE to SingleCellAssay...")
  sca <- MAST::SceToSingleCellAssay(sce, check_sanity = check_sanity)

  # Basic filtering
  expr <- as.matrix(SummarizedExperiment::assay(sca))

  #filter the genes that are not in all combinations

  if (random) {
    # Choose the sample-level variable
    sample_var <- random_var
    cd <- SummarizedExperiment::colData(sca)
    samples <- cd[[sample_var]]

    # Group cells by sample: split the cells per samples
    # Note: seq_len(ncol(expr)) gives cell indices: 1, 2, ..., #cells.
    sample_list <- split(seq_len(ncol(expr)), samples)

    # Determine, for each sample, whether each gene is detected
    detected_in_sample <- vapply(
      sample_list,
      function(idx) rowSums(expr[, idx, drop = FALSE] > 0) > 0,
      logical(nrow(expr))
    )

    # Ensure matrix shape: if only one sample, vapply returns a vector
    if (is.vector(detected_in_sample)) detected_in_sample <- matrix(
      detected_in_sample,
      ncol = 1
    )

    # Filter genes detected in ≥ minimum number of samples, in this case 5.
    min_samples <- 5
    keep_overall <- rowSums(detected_in_sample) >= min_samples

    # Require some expression in every group×time combo
    if (all(c("group", "time") %in% colnames(cd))) {
      combo <- interaction(cd$group, cd$time, drop = TRUE)
    } else if ("Condition" %in% colnames(cd)) {
      combo <- cd$Condition
    } else {
      combo <- factor(rep("All", ncol(expr)))
      message("No 'group' or 'time' columns found;
              using all cells as one combo.")
    }
    combo_list <- split(seq_len(ncol(expr)), combo)

    detected_in_combo <- vapply(
      combo_list,
      function(idx) rowSums(expr[, idx, drop = FALSE] > 0) > 0,
      logical(nrow(expr))
    )
    if (is.vector(detected_in_combo)) detected_in_combo <- matrix(
      detected_in_combo,
      ncol = 1
    )
    keep_combo <- rowSums(detected_in_combo) == ncol(detected_in_combo)

    # Ensure enough non-zero observations for random effects
    n_levels <- nlevels(samples)
    nz_total <- rowSums(expr > 0)
    keep_nz <- nz_total > n_levels  # glmer requires obs > levels

    # combine all random-effect filters
    keep <- keep_overall & keep_combo & keep_nz

    if (verbose) {
      message(
        "Keeping ", sum(keep),
        " genes after all random-effect filters (≥", min_samples,
        " samples, all combos, and > ", n_levels, " nonzero cells)."
      )
    }
  } else {
    # simpler case for fixed effects
    keep <- rowSums(expr > 0) > 0
    if (verbose) message("Keeping ", sum(keep), " non-zero genes out of ",
                         length(keep))
  }


  sca <- sca[keep, ]

  # Decide method automatically unless overridden
  if (is.null(method)) {
    method <- if (random) "glmer" else "glm"
  }

  # Fit the model
  if (verbose) {
    message(
      "Fitting MAST model (", method, ") with formula: ",
      deparse(formula)
    )
  }
  zlmfit <<- tryCatch(
    MAST::zlm(
      formula,
      sca,
      method = method,
      ebayes = ebayes
    ),
    error = function(e) {
      message("MAST model failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(zlmfit)) {
    return(NULL)
  }

  # Summarize and extract contrast results
  s <- tryCatch(
    MAST::summary(zlmfit, doLRT = lrt),
    error = function(e) {
      message("MAST summary failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(s)) {
    return(NULL)
  }

  dt <- s$datatable
  if (!"contrast" %in% names(dt)) {
    message("Summary did not include 'contrast' column — check your LRT name.")
    return(NULL)
  }
  View(dt |>
          dplyr::filter(contrast == lrt, component == "logFC"))

  fc <- dt |>
    dplyr::filter(contrast == lrt, component == "logFC") |>
    dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

  pval <- dt |>
    dplyr::filter(contrast == lrt, component == "H") |>
    dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

  res <- merge(fc, pval, by = "gene")
  res$padj <- p.adjust(res$pvalue, method = "fdr")
  res <- res[order(res$padj), ]
  res <- res |>
    dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

  if (verbose) message("MAST run completed: ", nrow(res), " genes tested.")
  return(res)
}

## DESeq2 pseudobulk -----------------------------------------------------------
run_deseq_pseudobulk <- function(
    sce,
    ids,
    design_formula,
    relevel = list(),
    # named list, e.g. list(Treatment = "Placebo", Time = "Baseline")
    min_counts = 10,
    sfType = "poscounts",
    contrast = NULL,
    # optional: DESeq2-style contrast or coefficient name
    alpha = 0.05,
    save_dds = FALSE,
    verbose = TRUE
) {
  # Aggregate to pseudobulk
  if (verbose) message("Aggregating to pseudobulk...")
  pb <- tryCatch(
    scater::aggregateAcrossCells(sce, ids = ids, statistics = "sum"),
    error = function(e) {
      message("AggregateAcrossCells failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(pb)) {
    return(NULL)
  }

  # Prepare metadata
  metadata <- as.data.frame(SummarizedExperiment::colData(pb))
  metadata <- metadata[, colSums(!is.na(metadata)) > 0, drop = FALSE]
  # drop all-NA cols
  metadata[] <- lapply(metadata, function(x) {
    if (is.character(x)) factor(x) else x
  })

  # Relevel factors if specified
  for (var in names(relevel)) {
    ref <- relevel[[var]]
    if (var %in% colnames(metadata)) {
      metadata[[var]] <- factor(metadata[[var]],
                                levels = unique(metadata[[var]])
      )
      metadata[[var]] <- stats::relevel(metadata[[var]], ref = ref)
      if (verbose) message("Relevelled ", var, " (reference = ", ref, ")")
    } else {
      message("Variable ", var, " not found in metadata; skipping relevel.")
    }
  }

  # Create DESeq dataset
  if (verbose) message("Creating DESeq2 dataset...")
  dds <- tryCatch(
    DESeq2::DESeqDataSetFromMatrix(
      countData = SummarizedExperiment::assay(pb, "counts"),
      colData = metadata,
      design = design_formula
    ),
    error = function(e) {
      message("DESeqDataSetFromMatrix failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dds)) {
    return(NULL)
  }

  # Filter low-count genes
  keep <- rowSums(DESeq2::counts(dds)) >= min_counts
  n_keep <- sum(keep)
  if (verbose) {
    message(
      "Keeping ", n_keep, " genes with ≥", min_counts,
      " counts."
    )
  }
  if (n_keep < 2 || ncol(dds) < 2) {
    message("Not enough data after filtering.")
    return(NULL)
  }
  dds <- dds[keep, ]

  # Fit model
  if (verbose) message("Running DESeq2...")
  dds <- tryCatch(
    DESeq2::DESeq(dds, sfType = sfType, quiet = !verbose),
    error = function(e) {
      message("DESeq() failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dds)) {
    return(NULL)
  }

  # Extract results
  if (verbose) message("Extracting results...")
  res <- tryCatch(
    if (is.null(contrast)) {
      DESeq2::results(dds, alpha = alpha)
    } else if (is.character(contrast)) {
      DESeq2::results(dds, name = contrast, alpha = alpha)
    } else {
      DESeq2::results(dds, contrast = contrast, alpha = alpha)
    },
    error = function(e) {
      message("results() failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(res)) {
    return(NULL)
  }

  # Clean and return results
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df[order(df$padj, na.last = TRUE), ]
  df <- df |>
    dplyr::mutate(
      lfcSE = ifelse(is.finite(stat), abs(log2FoldChange / stat), NA_real_)
    )

  if (save_dds) attr(df, "dds") <- dds

  if (verbose) message("DESeq2 run completed: ", nrow(df), " genes tested.")
  return(df)
}


## Combination of all of them:--------------------------------------------------
run_methods <- function(
    sce,
    condition_var = "Condition", # column for condition (e.g. Treatment)
    sample_var = "Sample", # column for random effect / pseudobulk aggregation
    design_formula_mast = ~ Condition,
    design_formula_mast_rand = ~ Condition + (1 | Sample),
    design_formula_deseq = ~ Condition,
    lrt_term = "ConditionTreatment", # term to test in MAST
    ids_vars = c("Sample", "Condition"), # columns for pseudobulk aggregation
    run_mast_fe = TRUE,
    run_mast_re = TRUE,
    run_deseq = TRUE
) {
  res_list <- list()

  # MAST (fixed effects)
  if (run_mast_fe) {
    message("Running MAST (fixed effects)...")
    mast_fe <- run_mast(
      sce,
      formula = design_formula_mast,
      lrt = lrt_term
    )

    if (!is.null(mast_fe)) {
      mast_fe$method <- "MAST (fixed)"
      mast_fe$formula <- deparse(design_formula_mast)
      res_list$mast_fe <- mast_fe
    }
  }

  # MAST (random effects)
  if (run_mast_re) {
    message("Running MAST (random effects)...")
    mast_re <- run_mast(
      sce,
      random_var = sample_var,
      formula = design_formula_mast_rand,
      lrt = lrt_term,
      random = TRUE,
      method = "glmer",
      ebayes = FALSE
    )

    if (!is.null(mast_re)) {
      mast_re$method <- "MAST (random)"
      mast_re$formula <- deparse(design_formula_mast_rand)
      res_list$mast_re <- mast_re
    }
  }

  # DESeq2 pseudobulk
  if (run_deseq) {
    message("Running DESeq2 (pseudobulk)...")
    ids <- SummarizedExperiment::colData(sce)[, ids_vars, drop = FALSE]
    deseq_res <- run_deseq_pseudobulk(
      sce,
      ids = ids,
      design_formula = design_formula_deseq,
      relevel = list(condition_var, sample_var)
    )
    if (!is.null(deseq_res)) {
      deseq_res$method <- "DESeq2"
      deseq_res$formula <- deparse(design_formula_deseq)
      res_list$deseq <- deseq_res
    }
  }

  message("Finished running selected methods.")
  return(res_list)
}


# Pipelines --------------------------------------------------------------------
## SplatPop pipeline with methods ----------------------------------------------
pipeline_splatpop_and_methods <- function(
    sce,
    n_samples = 20,
    n_snps = 1000,
    params = NULL,
    run_mast_fe = TRUE,
    run_mast_re = TRUE,
    run_deseq = TRUE,
    design_formula_mast = ~Condition,
    design_formula_mast_rand = ~ Condition + (1 | Sample),
    design_formula_deseq = ~Condition
) {

  splat_results <- run_methods(
    sce = sce,
    condition_var = "Condition",
    sample_var = "Sample",
    design_formula_mast = design_formula_mast,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = "ConditionTreatment",
    ids_vars = "Sample"
  )
}

## Jan pipeline with methods ---------------------------------------------------
pipeline_jan_and_methods <- function(sce,
                                     nSubjsPerGroup = 20,
                                     nCellsPerSubjectTimepoint = 100,
                                     nGenes = 100,
                                     run_mast_fe = TRUE,
                                     run_mast_re = TRUE,
                                     run_deseq = TRUE,
                                     design_formula = ~ group * time,
                                     design_formula_mast_rand =
                                       ~ group * time + (1 | sampleID),
                                     design_formula_deseq = ~ group * time) {
  jan_results <- run_methods(
    sce = sce,
    sample_var = "sampleID",
    design_formula_mast = design_formula,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = "groupTreatment:timeEoT",
    ids_vars = c("subjectID", "time"),
    run_mast_fe = run_mast_fe,
    run_mast_re = run_mast_re,
    run_deseq = run_deseq
  )
  jan_results
}

## rescueSim pipeline with methods ---------------------------------------------
pipeline_rescue_and_methods <- function(sce,
                                        nSubjsPerGroup = 5,
                                        run_mast = TRUE,
                                        run_mast_ran = TRUE,
                                        run_deseq = TRUE,
                                        design_formula_mast =
                                          ~ group * time,
                                        design_formula_mast_rand =
                                          ~  group * time + (1 | sampleID),
                                        design_formula_deseq = ~ group * time) {

  rescue_results <- run_methods(
    sce = sce,
    sample_var = "sampleID",
    design_formula_mast = design_formula_mast,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = "groupTreatment:timeEoT",
    ids_vars = c("subjectID", "time"),
    run_mast_fe = run_mast,
    run_mast_re = run_mast_ran,
    run_deseq = run_deseq
  )
  rescue_results
}





# Run Simulations  -------------------------------------------------------------
board <- pins::board_local("nrgn-data")

# Read pins with simulations
simulation <- pins::pin_read(board, "simulations")
sim_splatpop <- simulation$`splatPop (Treatment`
sim_jan <- simulation$`splatter (Treatment + Time`
sim_rescue <- simulation$`rescueSim (Treatment + Time`


## SplatPop --------------------------------------------------------------------
splatpop_out <- pipeline_splatpop_and_methods(sce = sim_splatpop)


## Jan -------------------------------------------------------------------------
jan_out <- pipeline_jan_and_methods(sce = sim_jan,
                                    run_mast_fe = FALSE,
                                    run_mast_re = TRUE,
                                    run_deseq = FALSE)
jan_out2 <- pipeline_jan_and_methods(sce = sim_jan,
                                     design_formula =
                                       ~ subjectID + time + group:time,
                                     design_formula_mast_rand =
                                       ~ subjectID + time + group:time +
                                       (1 | sampleID),
                                     design_formula_deseq =
                                      ~ subjectID + time + group:time)
jan_out3 <- pipeline_jan_and_methods(sce = sim_jan,
                                     design_formula =
                                       ~ subjectID + group*time,
                                     design_formula_mast_rand =
                                       ~ subjectID + time + group:time +
                                       (1 | subjectID:sampleID) + (1|sampleID),
                                     design_formula_deseq =
                                       ~ subjectID + group*time)
##  rescueSim ------------------------------------------------------------------
resc_out <- pipeline_rescue_and_methods(sce = sim_rescue)





# Pin data ---------------------------------------------------------------------
splatpop_out$mast_fe$analysis <- "splatpop"
splatpop_out$mast_re$analysis <- "splatpop"
splatpop_out$deseq$analysis   <- "splatpop"

resc_out$mast_fe$analysis <- "rescueSim"
resc_out$mast_re$analysis <- "rescueSim"
resc_out$deseq$analysis   <- "rescueSim"

jan_out$mast_fe$analysis  <- "jan"
jan_out$mast_re$analysis  <- "jan"
jan_out$deseq$analysis    <- "jan"

method_data <- rbind(
  splatpop_out$mast_fe, splatpop_out$mast_re, splatpop_out$deseq,
  resc_out$mast_fe, resc_out$mast_re, resc_out$deseq,
  jan_out$mast_fe, jan_out$mast_re, jan_out$deseq,
  fill = TRUE
)

board <- pins::board_local("nrgn-data")

pins::pin_write(method_data,
                name = "method_data",
                type = "rds",
                board = board,
                versioned = FALSE
)


#inspect structure of the data--------------------------------------------------
method_data  <- pins::pin_read(board, "method_data")
simulation <- pins::pin_read(board, "simulations")
summarize_sce_structure <- function(sce, group_var = NULL, time_var = NULL,
                                    subject_var = NULL) {
  df <- as.data.frame(SummarizedExperiment::colData(sce))

  out <- list()

  if (!is.null(group_var)) {
    out$groups <- df |>
      dplyr::count(!!dplyr::sym(group_var)) |>
      dplyr::arrange(desc(n))
  }

  if (!is.null(time_var)) {
    out$time <- df |>
      dplyr::count(!!dplyr::sym(time_var)) |>
      dplyr::arrange(desc(n))
  }

  if (!is.null(group_var) && !is.null(time_var)) {
    out$group_time <- df |>
      dplyr::count(!!dplyr::sym(group_var), !!dplyr::sym(time_var)) |>
      tidyr::pivot_wider(names_from = time_var, values_from = n)
  }

  if (!is.null(subject_var)) {
    out$subjects <- df |>
      dplyr::count(!!dplyr::sym(subject_var)) |>
      dplyr::arrange(desc(n))

    out$subjects_group <- df |>
      dplyr::count(!!dplyr::sym(subject_var), !!dplyr::sym(group_var))
  }

  return(out)
}
print_structure <- function(struct) {
  for (name in names(struct)) {
    cat("\n==== ", name, " ====\n")
    print(struct[[name]])
  }
}
sce <- simulation$Treatment

struct_splatpop <- summarize_sce_structure(
  sce,
  group_var = "Condition"
)

print_structure(struct_splatpop)

sce_rescue <- simulation$`Treatment + Time (rescueSim)`
struct_rescue <- summarize_sce_structure(
  sce_rescue,
  group_var = "group",
  time_var = "time",
  subject_var = "subjectID"
)
print_structure(struct_rescue)
sce_jan <- simulation$`Treatment + Time (splatter)`
struct_jan <- summarize_sce_structure(
  sce_jan,
  group_var = "group",
  time_var = "time",
  subject_var = "SubjectID"
)
print_structure(struct_jan)
