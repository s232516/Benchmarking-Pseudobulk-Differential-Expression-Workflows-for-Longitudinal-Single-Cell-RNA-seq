library(BiocParallel)
# Methods ----------------------------------------------------------------------
## MAST ------------------------------------------------------------------------
run_mast <- function(
  sce,
  formula = ~Condition,
  lrt = "ConditionTreatment",
  random_var = NULL,
  random = FALSE, # switch between fixed and random effects
  method = NULL, # override if needed
  ebayes = TRUE,

  check_sanity = FALSE,
  verbose = TRUE, # print progress messages
  min_samples = 5
) {
  # Convert SCE to SingleCellAssay
  if (verbose) message("Converting SCE to SingleCellAssay...")
  expr <- as.matrix(SummarizedExperiment::assay(sce))
  if (random) {
    # Choose the sample-level variable
    sample_var <- random_var
    cd <- SummarizedExperiment::colData(sce)
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
    if (is.vector(detected_in_sample)) {
      detected_in_sample <- matrix(
        detected_in_sample,
        ncol = 1
      )
    }

    # Filter genes detected in ≥ minimum number of samples, in this case 5.
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
    if (is.vector(detected_in_combo)) {
      detected_in_combo <- matrix(
        detected_in_combo,
        ncol = 1
      )
    }
    keep_combo <- rowSums(detected_in_combo) == ncol(detected_in_combo)

    # Ensure enough non-zero observations for random effects
    n_levels <- nlevels(samples)
    nz_total <- rowSums(expr > 0)
    keep_nz <- nz_total > n_levels # glmer requires obs > levels

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
    if (verbose) {
      message(
        "Keeping ", sum(keep), " non-zero genes out of ",
        length(keep)
      )
    }
  }

  sce <- sce[keep, ]

  # Log2(CPM + 1) normalization
  normcounts(sce) <- apply(SingleCellExperiment::counts(sce), 2, function(x) {
    log2((1e6 * x / sum(x)) + 1)
  })

  # Convert to SingleCellAssay
  sca <- suppressMessages(MAST::FromMatrix(as.matrix(normcounts(sce)),
                                           colData(sce)))
  # Calculate cellular detection rate
  cdr <- colSums(SummarizedExperiment::assay(sca) > 0)
  colData(sca)$cngeneson <- scale(cdr)


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
  zlmfit <- tryCatch(
    MAST::zlm(
      formula,
      sca,
      method = method,
      ebayes = ebayes,
      silent = TRUE
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
  relevel = NULL,
  # named list, e.g. list(Treatment = "Placebo", Time = "Baseline")
  min_counts = 10,
  sfType = "ratio",
  lrt_term = NULL,
  # optional: MAST-style contrast or coefficient name
  alpha = 0.05,
  save_dds = FALSE,
  verbose = TRUE,
  filter = FALSE
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
  if (!is.null(relevel)) {
    for (var in names(relevel)) {
      ref <- relevel[[var]]
      if (var %in% colnames(metadata)) {
        metadata[[var]] <- relevel(metadata[[var]], ref = ref)
        if (verbose) {
          message("Relevelled ", var, " (reference = ", ref, ")")
        }
      }
    }
  }

  if (filter) {
    if (verbose) message("Applying MAST-style filtering on pseudobulk...")

    # Extract pseudobulk counts
    expr <- SummarizedExperiment::assay(pb, "counts")

    # Metadata of pseudobulk samples
    md <- metadata

    # Identify condition variable automatically (first in design)
    cond_var <- all.vars(design_formula)[1]
    cond <- md[[cond_var]]
    cond <- droplevels(cond)

    # 1. DETECTION PER SAMPLE (gene expressed if count > 0)
    detected_by_sample <- expr > 0

    # 2. FILTER: gene detected in at least 1 sample per condition
    detected_by_condition <- sapply(levels(cond), function(gr) {
      idx <- which(cond == gr)
      rowSums(detected_by_sample[, idx, drop = FALSE]) > 0
    })

    # If only one condition, sapply returns a vector
    if (is.vector(detected_by_condition)) {
      detected_by_condition <- matrix(detected_by_condition, ncol = 1)
    }

    keep_condition <- rowSums(detected_by_condition) == ncol(detected_by_condition)

    # 3. FILTER: gene must have enough non-zero counts overall
    keep_nz <- rowSums(detected_by_sample) >= 2 # like MAST's min cells detected

    # 4. FILTER: gene must have reasonable total count
    keep_count <- rowSums(expr) >= min_counts

    # Combine all
    keep <- keep_condition & keep_nz & keep_count

    if (verbose) {
      message("Keeping ", sum(keep), " genes after filtering:")
      message(" - detected in ≥1 sample in EACH condition")
      message(" - detected in ≥2 pseudobulk samples total")
      message(" - total counts ≥ ", min_counts)
    }

    pb <- pb[keep, ]
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

  dseq_names <- DESeq2::resultsNames(dds)

  if (verbose) {
    message("\nDESeq2 coefficient names:")
    print(dseq_names)
  }

  if (!is.null(lrt_term)) {
    search_term <- lrt_term

    # find the match in DESeq2 coefficient names
    coef_name <- grep(search_term, dseq_names, value = TRUE)

    if (length(coef_name) == 0) {
      stop(
        "Could not match MAST term '", lrt_term,
        "' to any DESeq2 coefficient. Search term = ", search_term
      )
    }
    if (length(coef_name) > 1) {
      stop(
        "MAST term '", lrt_term, "' maps to multiple DESeq2 coefficients: ",
        paste(coef_name, collapse = ", ")
      )
    }


    if (verbose) {
      message("\nMAST term: ", lrt_term)
      message("Mapped DESeq2 coefficient: ", coef_name)
    }

    res <- DESeq2::results(dds, name = coef_name, alpha = alpha)
  } else {
    ### fallback (should rarely be used)
    res <- DESeq2::results(dds, alpha = alpha)
  }

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

  # In case we want to access to the DESeq object:
  if (save_dds) attr(df, "dds") <- dds

  if (verbose) message("DESeq2 run completed: ", nrow(df), " genes tested.")
  return(df)
}


## Combination of all of them:--------------------------------------------------
run_methods <- function(
  sce,
  relevel = list(Condition = "Placebo"),
  sample_var = "Sample", # column for random effect / pseudobulk aggregation
  design_formula_mast = ~Condition + cngeneson,
  design_formula_mast_rand = ~ Condition + cngeneson + (1 | Sample),
  design_formula_deseq = ~Condition,
  lrt_term = "ConditionTreatment", # term to test in MAST
  lrt_term_deseq = "Condition", # term to test in DESeq2
  ids_vars = c("Sample", "Condition"), # columns for pseudobulk aggregation
  run_mast_fe = TRUE,
  run_mast_re = TRUE,
  run_deseq = TRUE,
  filter = FALSE
) {
  res_list <- list()

  # MAST (fixed effects)
  if (run_mast_fe) {
    message("Running MAST (fixed effects)...")
    mast_fe <- run_mast(
      sce,
      formula = design_formula_mast,
      lrt = lrt_term,
      random_var = NULL,
      random = FALSE,
      ebayes = FALSE
    )
    message("MAST (fixed effects) completed")

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
    message("MAST (random effects) completed")

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
    message(
      "Using IDs for pseudobulk aggregation: ",
      paste(ids_vars, collapse = ", ")
    )
    deseq_res <- run_deseq_pseudobulk(
      sce,
      ids = ids,
      design_formula = design_formula_deseq,
      relevel = relevel,
      lrt_term = lrt_term_deseq,
      filter = filter
    )
    message("DESeq completed")
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
  run_mast_fe = TRUE,
  run_mast_re = TRUE,
  run_deseq = TRUE,
  design_formula_mast = ~Condition + cngeneson,
  design_formula_mast_rand = ~ Condition + cngeneson + (1 | Sample),
  design_formula_deseq = ~Condition
) {
  splat_results <- run_methods(
    sce = sce,
    sample_var = "Sample",
    design_formula_mast = design_formula_mast,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = "ConditionTreatment",
    lrt_term_deseq = "Condition_Treatment_vs_Placebo",
    ids_vars = "Sample",
    relevel = list(Condition = "Placebo"),
    run_mast_fe = run_mast_fe,
    run_mast_re = run_mast_re,
    run_deseq = run_deseq
  )
}

## Jan pipeline with methods ---------------------------------------------------
pipeline_jan_and_methods <- function(
  sce,
  run_mast_fe = TRUE,
  run_mast_re = TRUE,
  run_deseq = TRUE,
  design_formula_mast = ~ group * time + cngeneson,
  design_formula_mast_rand = ~ group * time + cngeneson + (1 | sampleID),
  design_formula_deseq = ~ group * time,
  lrt_term = "groupTreatment:timeEoT",
  lrt_term_deseq = "groupTreatment.timeEoT",
  relevel = list(group = "Placebo", time = "Baseline"),
  filter = FALSE
) {
  jan_results <- run_methods(
    sce = sce,
    sample_var = "sampleID",
    design_formula_mast = design_formula_mast,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = lrt_term,
    lrt_term_deseq = lrt_term_deseq,
    ids_vars = c("subjectID", "time"),
    relevel = list(group = "Placebo", time = "Baseline"),
    run_mast_fe = run_mast_fe,
    run_mast_re = run_mast_re,
    run_deseq = run_deseq,
    filter = filter
  )

  jan_results
}


## rescueSim pipeline with methods ---------------------------------------------
pipeline_rescue_and_methods <- function(
  sce,
  run_mast_fe = TRUE,
  run_mast_re = TRUE,
  run_deseq = TRUE,
  design_formula_mast = ~ group * time + cngeneson,
  design_formula_mast_rand = ~ group * time + cngeneson + (1 | sampleID),
  design_formula_deseq = ~ group * time,
  lrt_term = "groupTreatment:timeEoT",
  lrt_term_deseq = "groupTreatment.timeEoT",
  relevel = list(group = "Placebo", time = "Baseline"),
  filter = FALSE
) {
  rescue_results <- run_methods(
    sce = sce,
    sample_var = "sampleID",
    design_formula_mast = design_formula_mast,
    design_formula_mast_rand = design_formula_mast_rand,
    design_formula_deseq = design_formula_deseq,
    lrt_term = lrt_term,
    lrt_term_deseq = lrt_term_deseq,
    ids_vars = c("subjectID", "time"),
    relevel = list(group = "Placebo", time = "Baseline"),
    run_mast_fe = run_mast_fe,
    run_mast_re = run_mast_re,
    run_deseq = run_deseq,
    filter = filter
  )

  rescue_results
}


# Run Simulations  -------------------------------------------------------------
board <- pins::board_local("nrgn-data")

# Read pins with simulations
simulation <- pins::pin_read(board, "simulations")
sim_splatpop <- simulation$splatpop
sim_jan <- simulation$jan
sim_rescue <- simulation$rescueSim
sim_rescue_123 <- simulation$rescueSim_123
sim_rescue_111 <- simulation$rescueSim_111
sim_rescue_222 <- simulation$rescueSim_222
sim_rescue_333 <- simulation$rescueSim_333


## SplatPop --------------------------------------------------------------------
splatpop_out <- pipeline_splatpop_and_methods(sce = sim_splatpop)


## Jan -------------------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------

# WHEN THERE IS NO GROUP FACTOR
sim_jan_ref <- sim_jan
cd <- SummarizedExperiment::colData(sim_jan_ref)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_jan_ref) <- cd

# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_jan_mod <- sim_jan
cd <- as.data.frame(SummarizedExperiment::colData(sim_jan_mod))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_jan_mod) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------
jan_out <- pipeline_jan_and_methods(sce = sim_jan)

#### MAST fixed only -----------------------------------------------------------
jan_out_mast <- pipeline_jan_and_methods(
  sce = sim_jan_ref,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)

jan_out_mast_cdr <- pipeline_jan_and_methods(
  sce = sim_jan_ref,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "groupTreatment:timeEoT"
)
#### MAST random only ----------------------------------------------------------

jan_out_mast_ra <- pipeline_jan_and_methods(
  sce = sim_jan_ref,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID
      + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

jan_out_mast_ra_cdr <- pipeline_jan_and_methods(
  sce = sim_jan_ref,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson
  + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

#### DESeq only ----------------------------------------------------------------
jan_out_deseq <- pipeline_jan_and_methods(
  sce = sim_jan_mod,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

jan_out_deseq_s <- pipeline_jan_and_methods(
  sce = sim_jan_ref,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)


##  rescueSim ------------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------
# WHEN THERE IS NO GROUP FACTOR
sim_rescue_ref <- sim_rescue
cd <- SummarizedExperiment::colData(sim_rescue_ref)
cd$time <- relevel(cd$time, ref = "EoT")

SummarizedExperiment::colData(sim_rescue_ref) <- cd
# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_rescue_mod <- sim_rescue
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------

resc_out <- pipeline_rescue_and_methods(sce = sim_rescue)

#### MAST fixed only -----------------------------------------------------------
resc_out_mast <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)
resc_out_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment"
)

#### MAST random only ----------------------------------------------------------

resc_out_mast_ra <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
      (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
resc_out_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

#### DESeq only ----------------------------------------------------------------

resc_out_deseq <- pipeline_rescue_and_methods(
  sce = sim_rescue_mod,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

resc_out_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)




##  rescueSim 123 --------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------
# WHEN THERE IS NO GROUP FACTOR
sim_rescue_ref_123 <- sim_rescue_123
cd <- SummarizedExperiment::colData(sim_rescue_ref_123)
cd$time <- relevel(cd$time, ref = "EoT")

SummarizedExperiment::colData(sim_rescue_ref_123) <- cd
# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_rescue_mod_123 <- sim_rescue_123
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_123))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_123) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------

resc_out_123 <- pipeline_rescue_and_methods(sce = sim_rescue_123)

#### MAST fixed only -----------------------------------------------------------
resc_out_123_mast <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_123_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment"
)

#### MAST random only ----------------------------------------------------------

resc_out_123_mast_re <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
      (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
resc_out_123_mast_re_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_123,
  run_mast = FALSE,
  run_mast_ran = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ group * time + subjectID + cngeneson +
      (1 | sampleID),
  lrt_term = "groupTreatment:timeEoT"
)

#### DESeq only ----------------------------------------------------------------

resc_out_123_deseq <- pipeline_rescue_and_methods(
  sce = sim_rescue_mod_123,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

resc_out_123_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)


##  rescueSim 111 --------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------
# WHEN THERE IS NO GROUP FACTOR
sim_rescue_ref_111 <- sim_rescue_111
cd <- SummarizedExperiment::colData(sim_rescue_ref_111)
cd$time <- relevel(cd$time, ref = "EoT")

SummarizedExperiment::colData(sim_rescue_ref_111) <- cd
# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_rescue_mod_111 <- sim_rescue_111
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_111))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_111) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------

resc_out_111 <- pipeline_rescue_and_methods(sce = sim_rescue_111)

#### MAST fixed only -----------------------------------------------------------
resc_out_111_mast <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_111_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment"
)

#### MAST random only ----------------------------------------------------------
resc_out_111_mast_ra <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
      (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_111_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
#### DESeq only ----------------------------------------------------------------

resc_out_111_deseq <- pipeline_rescue_and_methods(
  sce = sim_rescue_mod_111,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

resc_out_111_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)


##  rescueSim 222 --------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------
# WHEN THERE IS NO GROUP FACTOR
sim_rescue_ref_222 <- sim_rescue_222
cd <- SummarizedExperiment::colData(sim_rescue_ref_222)
cd$time <- relevel(cd$time, ref = "EoT")

SummarizedExperiment::colData(sim_rescue_ref_222) <- cd
# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_rescue_mod_222 <- sim_rescue_222
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_222))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_222) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------

resc_out_222 <- pipeline_rescue_and_methods(sce = sim_rescue_222)

#### MAST fixed only -----------------------------------------------------------
resc_out_222_mast <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_222_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment"
)

#### MAST random only ----------------------------------------------------------
resc_out_222_mast_ra <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_222_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
#### DESeq only ----------------------------------------------------------------

resc_out_222_deseq <- pipeline_rescue_and_methods(
  sce = sim_rescue_mod_222,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

resc_out_222_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)

##  rescueSim 333 --------------------------------------------------------------
### Define other dataframes ( columns in another way) --------------------------
# WHEN THERE IS NO GROUP FACTOR
sim_rescue_ref_333 <- sim_rescue_333
cd <- SummarizedExperiment::colData(sim_rescue_ref_333)
cd$time <- relevel(cd$time, ref = "EoT")

SummarizedExperiment::colData(sim_rescue_ref_333) <- cd
# CODED TIME AND TIME.GROUP , 0NLY 0 AND 1S
sim_rescue_mod_333 <- sim_rescue_333
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_333))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_333) <- S4Vectors::DataFrame(cd)

### Run different versions of the pipeline -------------------------------------

resc_out_333 <- pipeline_rescue_and_methods(sce = sim_rescue_333)

#### MAST fixed only -----------------------------------------------------------
resc_out_333_mast <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID,
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_333_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment"
)

#### MAST random only ----------------------------------------------------------
resc_out_333_mast_ra <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_333_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
#### DESeq only ----------------------------------------------------------------

resc_out_333_deseq <- pipeline_rescue_and_methods(
  sce = sim_rescue_mod_333,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ is.time + is.grouptime + subjectID,
  lrt_term_deseq = "is.grouptime",
  relevel = NULL
)

resc_out_333_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = FALSE,
  run_mast_re = FALSE,
  run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline")
)

# Pin data ---------------------------------------------------------------------
splatpop_out$mast_fe$analysis <- "splatpop"
splatpop_out$mast_re$analysis <- "splatpop"
splatpop_out$deseq$analysis <- "splatpop"

resc_out$mast_fe$analysis <- "rescueSim" # 1
resc_out$mast_re$analysis <- "rescueSim"
resc_out$deseq$analysis <- "rescueSim"
resc_out_mast$mast_re$analysis <- "rescueSim" # 3
resc_out_mast_cdr$mast_fe$analysis <- "rescueSim" # 2
resc_out_mast_ra$mast_re$analysis <- "rescueSim" # 4
resc_out_mast_ra_cdr$mast_re$analysis <- "rescueSim" # 5
resc_out_deseq$deseq$analysis <- "rescueSim" # 6
resc_out_deseq_s$deseq$analysis <- "rescueSim" # 7

resc_out_123$mast_fe$analysis <- "rescueSim (123)" # 1
resc_out_123$mast_re$analysis <- "rescueSim (123)"
resc_out_123$deseq$analysis <- "rescueSim (123)"
resc_out_123_mast$mast_re$analysis <- "rescueSim (123)" # 3
resc_out_123_mast_cdr$mast_fe$analysis <- "rescueSim (123)" # 2
resc_out_123_mast_ra$mast_re$analysis <- "rescueSim (123)" # 4
resc_out_123_mast_ra_cdr$mast_re$analysis <- "rescueSim (123)" # 5
resc_out_123_deseq$deseq$analysis <- "rescueSim (123)" # 6
resc_out_123_deseq_s$deseq$analysis <- "rescueSim (123)" # 7


resc_out_111$mast_fe$analysis <- "rescueSim (111)" # 1
resc_out_111$mast_re$analysis <- "rescueSim (111)"
resc_out_111$deseq$analysis <- "rescueSim (111)"
resc_out_111_mast$mast_re$analysis <- "rescueSim (111)" # 3
resc_out_111_mast_cdr$mast_fe$analysis <- "rescueSim (111)" # 2
resc_out_111_mast_ra$mast_re$analysis <- "rescueSim (111)" # 4
resc_out_111_mast_ra_cdr$mast_re$analysis <- "rescueSim (111)" # 5
resc_out_111_deseq$deseq$analysis <- "rescueSim (111)" # 6
resc_out_111_deseq_s$deseq$analysis <- "rescueSim (111)" # 7

resc_out_222$mast_fe$analysis <- "rescueSim (222)" # 1
resc_out_222$mast_re$analysis <- "rescueSim (222)"
resc_out_222$deseq$analysis <- "rescueSim (222)"
resc_out_222_mast$mast_re$analysis <- "rescueSim (222)" # 3
resc_out_222_mast_cdr$mast_fe$analysis <- "rescueSim (222)" # 2
resc_out_222_mast_ra$mast_re$analysis <- "rescueSim (222)" # 4
resc_out_222_mast_ra_cdr$mast_re$analysis <- "rescueSim (222)" # 5
resc_out_222_deseq$deseq$analysis <- "rescueSim (222)" # 6
resc_out_222_deseq_s$deseq$analysis <- "rescueSim (222)" # 7

resc_out_333$mast_fe$analysis <- "rescueSim (333)" # 1
resc_out_333$mast_re$analysis <- "rescueSim (333)"
resc_out_333$deseq$analysis <- "rescueSim (333)"
resc_out_333_mast$mast_re$analysis <- "rescueSim (333)" # 3
resc_out_333_mast_cdr$mast_fe$analysis <- "rescueSim (333)" # 2
resc_out_333_mast_ra$mast_re$analysis <- "rescueSim (333)" # 4
resc_out_333_mast_ra_cdr$mast_re$analysis <- "rescueSim (333)" # 5
resc_out_333_deseq$deseq$analysis <- "rescueSim (333)" # 6
resc_out_333_deseq_s$deseq$analysis <- "rescueSim (333)" # 7

jan_out$mast_fe$analysis <- "jan"
jan_out$mast_re$analysis <- "jan"
jan_out$deseq$analysis <- "jan"
jan_out_mast$mast_re$analysis <- "jan"
jan_out_mast_cdr$mast_fe$analysis <- "jan"
jan_out_mast_ra$mast_re$analysis <- "jan"
jan_out_mast_ra_cdr$mast_re$analysis <- "jan"
jan_out_deseq$deseq$analysis <- "jan"
jan_out_deseq_s$deseq$analysis <- "jan"


method_data <- rbind(
  splatpop_out$mast_fe, splatpop_out$mast_re, splatpop_out$deseq,

  resc_out$mast_fe, resc_out$mast_re, resc_out$deseq,
  resc_out_mast$mast_re, resc_out_mast_cdr$mast_fe, resc_out_mast_ra$mast_re,
  resc_out_mast_ra_cdr$mast_re, resc_out_deseq$deseq, resc_out_deseq_s$deseq,

  jan_out$mast_fe, jan_out$mast_re, jan_out$deseq,
  jan_out_mast$mast_re, jan_out_mast_cdr$mast_fe, jan_out_mast_ra$mast_re,
  jan_out_mast_ra_cdr$mast_re, jan_out_deseq$deseq, jan_out_deseq_s$deseq,

  resc_out_123$mast_fe, resc_out_123$mast_re, resc_out_123$deseq,
  resc_out_123_mast$mast_re, resc_out_123_mast_cdr$mast_fe, resc_out_123_mast_ra$mast_re,
  resc_out_123_mast_ra_cdr$mast_re, resc_out_123_deseq$deseq, resc_out_123_deseq_s$deseq,

  resc_out_111$mast_fe, resc_out_111$mast_re, resc_out_111$deseq,
  resc_out_111_mast$mast_re, resc_out_111_mast_cdr$mast_fe, resc_out_111_mast_ra$mast_re,
  resc_out_111_mast_ra_cdr$mast_re, resc_out_111_deseq$deseq, resc_out_111_deseq_s$deseq,

  resc_out_222$mast_fe, resc_out_222$mast_re, resc_out_222$deseq,
  resc_out_222_mast$mast_re, resc_out_222_mast_cdr$mast_fe, resc_out_222_mast_ra$mast_re,
  resc_out_222_mast_ra_cdr$mast_re, resc_out_222_deseq$deseq, resc_out_222_deseq_s$deseq,

  resc_out_333$mast_fe, resc_out_333$mast_re, resc_out_333$deseq,
  resc_out_333_mast$mast_re, resc_out_333_mast_cdr$mast_fe, resc_out_333_mast_ra$mast_re,
  resc_out_333_mast_ra_cdr$mast_re, resc_out_333_deseq$deseq, resc_out_333_deseq_s$deseq,
  fill = TRUE
)


board <- pins::board_local("nrgn-data")

pins::pin_write(method_data,
  name = "method_data_final",
  type = "rds",
  board = board,
  versioned = TRUE
)

# Merge ground_truth with method data:

method_data <- pins::pin_read(board, "method_data_final")
ground_truth <- pins::pin_read(board, "ground_truth")

method_data$analysis <- gsub(
  pattern = "rescueSim \\((\\d+)\\)",
  replacement = "rescueSim_\\1",
  method_data$analysis
)
compare_data <- dplyr::bind_rows(
  lapply(names(ground_truth), function(sim) {
    dplyr::left_join(
      method_data[method_data$analysis == sim, ],
      ground_truth[[sim]],
      by = "gene"
    )
  })
)

pins::pin_write(
  compare_data,
  name = "compare_data_final",
  type = "rds",
  board = board,
  versioned = TRUE
)

# GROUND TRUTH HIST -----------------------------------------------------------

# Create a pin for how many of the genes where in the conf intervals for the
# rescuesim datasets
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
    total_up = sum(df$regulation == "up"),
    in_ci_up = sum(df$in_ci & df$regulation == "up", na.rm = TRUE),
    total_down = sum(df$regulation == "down"),
    in_ci_down = sum(df$in_ci & df$regulation == "down", na.rm = TRUE),
    total_not_sig = sum(df$regulation == "not_sig"),
    in_ci_not_sig = sum(df$in_ci & df$regulation == "not_sig", na.rm = TRUE)
  )
}

ci_level <- 0.95 # you can loop over several if you want

rescue_ci_summary <- compare_data |>
  dplyr::filter(grepl("^rescueSim", analysis)) |>
  dplyr::group_by(analysis, method, formula) |>
  dplyr::group_modify(~ summarise_ci_coverage(.x, conf_level = ci_level)) |>
  dplyr::ungroup() |>
  dplyr::mutate(conf_level = ci_level)

rescue_ci_summary <- rescue_ci_summary |>
  dplyr::mutate(
    prop_ci_total    = in_ci_total / total_genes,
    prop_ci_up       = in_ci_up / total_up,
    prop_ci_down     = in_ci_down / total_down,
    prop_ci_not_sig  = in_ci_not_sig / total_not_sig
  )

board <- pins::board_local("nrgn-data")

pins::pin_write(
  rescue_ci_summary,
  name = "rescueSim_ci_new",
  type = "rds",
  board = board,
  versioned = TRUE
)
