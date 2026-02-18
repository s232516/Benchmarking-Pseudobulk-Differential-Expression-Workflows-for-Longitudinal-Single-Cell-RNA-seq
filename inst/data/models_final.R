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
    norm = TRUE, # whether to do log-CPM normalization or log-raw counts
    verbose = TRUE, # print progress messages
    min_samples = 5
    ) {
      t_start <- proc.time()[["elapsed"]]
      # Convert SCE to SingleCellAssay
      if (verbose) message("Converting SCE to SingleCellAssay...")
      expr <- as.matrix(SingleCellExperiment::counts(sce))
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

      # Build expression matrix for MAST input
      counts_mat <- as.matrix(SingleCellExperiment::counts(sce))

      if (norm) {
        # Log2(CPM + 1) normalization
        expr_mat <- apply(counts_mat, 2, function(x) {
          s <- sum(x)
          if (s == 0) return(rep(0, length(x)))
          log2((1e6 * x / s) + 1)
        })
      } else {
        # log2(counts + 1): log-transformed but NOT normalized
        expr_mat <- log2(counts_mat + 1)
        # (numerically stable alternative: log1p(counts_mat)/log(2))
      }

      # Convert to SingleCellAssay using the chosen expression matrix
      sca <- suppressMessages(
        MAST::FromMatrix(as.matrix(expr_mat), SummarizedExperiment::colData(sce))
      )

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

      # Continuous-part coefficient (what you asked for)
      # p-value column for z-tests can vary by MAST version; grab the common one if present.
      cont_pcol <- c("Pr(>|z|)", "Pr(>z)")[
        c("Pr(>|z|)", "Pr(>z)") %in% names(dt)
      ][1]

      cont <- dt |>
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
      res <- res[order(res$padj), ]
      res <- res |>
        dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

      res$betaC_SE <- with(res, ifelse(!is.na(zC) & zC != 0, abs(betaC / zC), NA_real_))

      t_end <- proc.time()[["elapsed"]]
      timing <- t_end - t_start

      res$time <- timing
      if (verbose) {
        message(
          "MAST run completed: ", nrow(res), " genes tested. ",
          "Time (s): total=", round(timing, 2)
        )
      }
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
    verbose = TRUE
    ) {
      t_start <- proc.time()[["elapsed"]]
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
            metadata[[var]] <- stats::relevel(metadata[[var]], ref = ref)
            if (verbose) {
              message("Relevelled ", var, " (reference = ", ref, ")")
            }
          }
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
      t_end <- proc.time()[["elapsed"]]
      timing <- t_end - t_start
      df$time <- timing

      # In case we want to access to the DESeq object:
      if (save_dds) attr(df, "dds") <- dds

      if (verbose) message("DESeq2 run completed: ", nrow(df), " genes tested. Time (s): ", round(timing, 2))
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
    norm = TRUE
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
          ebayes = FALSE,
          norm = norm
        )
        message("MAST (fixed effects) completed")

        if (!is.null(mast_fe)) {
          mast_fe$method <- "MAST (fixed)"
          mast_fe$formula <- deparse(design_formula_mast)
          mast_fe$norm <- norm
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
          ebayes = FALSE,
          norm = norm
        )
        message("MAST (random effects) completed")

        if (!is.null(mast_re)) {
          mast_re$method <- "MAST (random)"
          mast_re$formula <- deparse(design_formula_mast_rand)
          mast_re$norm <- norm
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
          lrt_term = lrt_term_deseq
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
    design_formula_deseq = ~Condition,
    norm = TRUE
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
        run_deseq = run_deseq,
        norm = norm)
    }

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
    norm = TRUE
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
        norm = norm
      )

      rescue_results
    }

# Run Simulations  -------------------------------------------------------------
# Read pins with simulations
simulation <- pins::pin_read(board, "simulations_final")

sim_splatpop    <- simulation$splatpop
sim_splatpop_11 <- simulation$splatpop_11
sim_splatpop_22 <- simulation$splatpop_22
sim_splatpop_33 <- simulation$splatpop_33
sim_splatpop_44 <- simulation$splatpop_44


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



## SplatPop --------------------------------------------------------------------
splatpop_out <- pipeline_splatpop_and_methods(sce = sim_splatpop)
splatpop_out_no_cdr <- pipeline_splatpop_and_methods(sce = sim_splatpop,
                                                     design_formula_mast =
                                                       ~Condition,
                                                     design_formula_mast_rand =
                                                       ~ Condition + (1 | Sample),
                                                     run_deseq = TRUE)
splatpop_out_no_norm <- pipeline_splatpop_and_methods(sce = sim_splatpop,
                                                      norm = FALSE,
                                                      run_deseq = FALSE)


# Seed 11
splatpop_out_11 <- pipeline_splatpop_and_methods(sce = sim_splatpop_11)

# Seed 22
splatpop_out_22 <- pipeline_splatpop_and_methods(sce = sim_splatpop_22)

# Seed 33
splatpop_out_33 <- pipeline_splatpop_and_methods(sce = sim_splatpop_33)

# Seed 44
splatpop_out_44 <- pipeline_splatpop_and_methods(sce = sim_splatpop_44)


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

#resc_out <- pipeline_rescue_and_methods(sce = sim_rescue)

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
resc_out_mast_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

resc_out_mast_ra_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

#resc_out_123 <- pipeline_rescue_and_methods(sce = sim_rescue_123)

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

resc_out_123_mast_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
)

#### MAST random only ----------------------------------------------------------

resc_out_123_mast_re <- pipeline_rescue_and_methods(
  sce = sim_rescue_123,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID +
      (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)
resc_out_123_mast_re_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
      (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment"
)

resc_out_123_mast_re_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_123,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ group * time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "groupTreatment:timeEoT",
  norm = FALSE
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

# resc_out_111 <- pipeline_rescue_and_methods(sce = sim_rescue_111)

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

resc_out_111_mast_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

resc_out_111_mast_ra_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

# resc_out_222 <- pipeline_rescue_and_methods(sce = sim_rescue_222)

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
resc_out_222_mast_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

resc_out_222_mast_ra_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

# resc_out_333 <- pipeline_rescue_and_methods(sce = sim_rescue_333)

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

resc_out_333_mast_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = TRUE,
  run_mast_re = FALSE,
  run_deseq = FALSE,
  design_formula_mast =
    ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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

resc_out_333_mast_ra_cdr_nn <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333,
  run_mast_fe = FALSE,
  run_mast_re = TRUE,
  run_deseq = FALSE,
  design_formula_mast_rand =
    ~ time + group:time + subjectID + cngeneson +
    (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = FALSE
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


## rescueSim NORM ----------------------------------------------------------------
## Base (rescueSim norm)
sim_rescue_ref_norm <- sim_rescue_norm
cd <- SummarizedExperiment::colData(sim_rescue_ref_norm)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_rescue_ref_norm) <- cd

sim_rescue_mod_norm <- sim_rescue_norm
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_norm))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_norm) <- S4Vectors::DataFrame(cd)

resc_out_norm_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_norm,
  run_mast_fe = TRUE, run_mast_re = FALSE, run_deseq = FALSE,
  design_formula_mast = ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_norm_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_norm,
  run_mast_fe = FALSE, run_mast_re = TRUE, run_deseq = FALSE,
  design_formula_mast_rand = ~ time + group:time + subjectID + cngeneson + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_norm_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_norm,
  run_mast_fe = FALSE, run_mast_re = FALSE, run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline"),
  norm = TRUE
)


## rescueSim_123 NORM -----------------------------------------------------------
sim_rescue_ref_123_norm <- sim_rescue_123_norm
cd <- SummarizedExperiment::colData(sim_rescue_ref_123_norm)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_rescue_ref_123_norm) <- cd

sim_rescue_mod_123_norm <- sim_rescue_123_norm
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_123_norm))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_123_norm) <- S4Vectors::DataFrame(cd)

resc_out_123_norm_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123_norm,
  run_mast_fe = TRUE, run_mast_re = FALSE, run_deseq = FALSE,
  design_formula_mast = ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_123_norm_mast_re_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123_norm,
  run_mast_fe = FALSE, run_mast_re = TRUE, run_deseq = FALSE,
  design_formula_mast_rand = ~ time + group:time + subjectID + cngeneson + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_123_norm_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_123_norm,
  run_mast_fe = FALSE, run_mast_re = FALSE, run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline"),
  norm = TRUE
)


## rescueSim_111 NORM -----------------------------------------------------------
sim_rescue_ref_111_norm <- sim_rescue_111_norm
cd <- SummarizedExperiment::colData(sim_rescue_ref_111_norm)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_rescue_ref_111_norm) <- cd

sim_rescue_mod_111_norm <- sim_rescue_111_norm
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_111_norm))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_111_norm) <- S4Vectors::DataFrame(cd)

resc_out_111_norm_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111_norm,
  run_mast_fe = TRUE, run_mast_re = FALSE, run_deseq = FALSE,
  design_formula_mast = ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_111_norm_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111_norm,
  run_mast_fe = FALSE, run_mast_re = TRUE, run_deseq = FALSE,
  design_formula_mast_rand = ~ time + group:time + subjectID + cngeneson + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_111_norm_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_111_norm,
  run_mast_fe = FALSE, run_mast_re = FALSE, run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline"),
  norm = TRUE
)


## rescueSim_222 NORM -----------------------------------------------------------
sim_rescue_ref_222_norm <- sim_rescue_222_norm
cd <- SummarizedExperiment::colData(sim_rescue_ref_222_norm)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_rescue_ref_222_norm) <- cd

sim_rescue_mod_222_norm <- sim_rescue_222_norm
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_222_norm))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_222_norm) <- S4Vectors::DataFrame(cd)

resc_out_222_norm_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222_norm,
  run_mast_fe = TRUE, run_mast_re = FALSE, run_deseq = FALSE,
  design_formula_mast = ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_222_norm_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222_norm,
  run_mast_fe = FALSE, run_mast_re = TRUE, run_deseq = FALSE,
  design_formula_mast_rand = ~ time + group:time + subjectID + cngeneson + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_222_norm_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_222_norm,
  run_mast_fe = FALSE, run_mast_re = FALSE, run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline"),
  norm = TRUE
)


## rescueSim_333 NORM -----------------------------------------------------------
sim_rescue_ref_333_norm <- sim_rescue_333_norm
cd <- SummarizedExperiment::colData(sim_rescue_ref_333_norm)
cd$time <- relevel(cd$time, ref = "EoT")
SummarizedExperiment::colData(sim_rescue_ref_333_norm) <- cd

sim_rescue_mod_333_norm <- sim_rescue_333_norm
cd <- as.data.frame(SummarizedExperiment::colData(sim_rescue_mod_333_norm))
cd$is.time <- ifelse(cd$time == "EoT", 1, 0)
cd$is.grouptime <- ifelse(cd$group == "Treatment" & cd$time == "EoT", 1, 0)
SummarizedExperiment::colData(sim_rescue_mod_333_norm) <- S4Vectors::DataFrame(cd)

resc_out_333_norm_mast_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333_norm,
  run_mast_fe = TRUE, run_mast_re = FALSE, run_deseq = FALSE,
  design_formula_mast = ~ time + group:time + subjectID + cngeneson,
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_333_norm_mast_ra_cdr <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333_norm,
  run_mast_fe = FALSE, run_mast_re = TRUE, run_deseq = FALSE,
  design_formula_mast_rand = ~ time + group:time + subjectID + cngeneson + (1 | sampleID),
  lrt_term = "timeEoT:groupTreatment",
  norm = TRUE
)

resc_out_333_norm_deseq_s <- pipeline_rescue_and_methods(
  sce = sim_rescue_ref_333_norm,
  run_mast_fe = FALSE, run_mast_re = FALSE, run_deseq = TRUE,
  design_formula_deseq = ~ time + group:time,
  lrt_term_deseq = "timeEoT.groupTreatment",
  relevel = list(group = "Placebo", time = "Baseline"),
  norm = TRUE
)


# Pin data ---------------------------------------------------------------------
splatpop_out$mast_fe$analysis <- "splatpop"
splatpop_out$mast_re$analysis <- "splatpop"
splatpop_out$deseq$analysis <- "splatpop"

splatpop_out_no_cdr$mast_fe$analysis <- "splatpop"
splatpop_out_no_cdr$mast_re$analysis <- "splatpop"

splatpop_out_no_norm$mast_fe$analysis <- "splatpop"
splatpop_out_no_norm$mast_re$analysis <- "splatpop"

#seed 1

splatpop_out_11$mast_fe$analysis <- "splatpop (11)"
splatpop_out_11$mast_re$analysis <- "splatpop (11)"
splatpop_out_11$deseq$analysis <- "splatpop (11)"

#seed 2

splatpop_out_22$mast_fe$analysis <- "splatpop (22)"
splatpop_out_22$mast_re$analysis <- "splatpop (22)"
splatpop_out_22$deseq$analysis <- "splatpop (22)"

#seed 3

splatpop_out_33$mast_fe$analysis <- "splatpop (33)"
splatpop_out_33$mast_re$analysis <- "splatpop (33)"
splatpop_out_33$deseq$analysis <- "splatpop (33)"

#seed 4

splatpop_out_44$mast_fe$analysis <- "splatpop (44)"
splatpop_out_44$mast_re$analysis <- "splatpop (44)"
splatpop_out_44$deseq$analysis <- "splatpop (44)"


resc_out_mast$mast_fe$analysis <- "rescueSim"
resc_out_mast_cdr$mast_fe$analysis <- "rescueSim"
resc_out_mast_cdr_nn$mast_fe$analysis <- "rescueSim"
resc_out_mast_ra$mast_re$analysis <- "rescueSim"
resc_out_mast_ra_cdr$mast_re$analysis <- "rescueSim"
resc_out_mast_ra_cdr_nn$mast_re$analysis <- "rescueSim"
resc_out_deseq$deseq$analysis <- "rescueSim"
resc_out_deseq_s$deseq$analysis <- "rescueSim"

resc_out_123_mast$mast_fe$analysis <- "rescueSim (123)"
resc_out_123_mast_cdr$mast_fe$analysis <- "rescueSim (123)"
resc_out_123_mast_cdr_nn$mast_fe$analysis <- "rescueSim (123)"
resc_out_123_mast_re$mast_re$analysis <- "rescueSim (123)"
resc_out_123_mast_re_cdr$mast_re$analysis <- "rescueSim (123)"
resc_out_123_mast_re_cdr_nn$mast_re$analysis <- "rescueSim (123)"
resc_out_123_deseq$deseq$analysis <- "rescueSim (123)"
resc_out_123_deseq_s$deseq$analysis <- "rescueSim (123)"


resc_out_111_mast$mast_fe$analysis <- "rescueSim (111)"
resc_out_111_mast_cdr$mast_fe$analysis <- "rescueSim (111)"
resc_out_111_mast_cdr_nn$mast_fe$analysis <- "rescueSim (111)"
resc_out_111_mast_ra$mast_re$analysis <- "rescueSim (111)"
resc_out_111_mast_ra_cdr$mast_re$analysis <- "rescueSim (111)"
resc_out_111_mast_ra_cdr_nn$mast_re$analysis <- "rescueSim (111)"
resc_out_111_deseq$deseq$analysis <- "rescueSim (111)"
resc_out_111_deseq_s$deseq$analysis <- "rescueSim (111)"

resc_out_222_mast$mast_fe$analysis <- "rescueSim (222)"
resc_out_222_mast_cdr$mast_fe$analysis <- "rescueSim (222)"
resc_out_222_mast_cdr_nn$mast_fe$analysis <- "rescueSim (222)"
resc_out_222_mast_ra$mast_re$analysis <- "rescueSim (222)"
resc_out_222_mast_ra_cdr$mast_re$analysis <- "rescueSim (222)"
resc_out_222_mast_ra_cdr_nn$mast_re$analysis <- "rescueSim (222)"
resc_out_222_deseq$deseq$analysis <- "rescueSim (222)"
resc_out_222_deseq_s$deseq$analysis <- "rescueSim (222)"

resc_out_333_mast$mast_fe$analysis <- "rescueSim (333)"
resc_out_333_mast_cdr$mast_fe$analysis <- "rescueSim (333)"
resc_out_333_mast_cdr_nn$mast_fe$analysis <- "rescueSim (333)"
resc_out_333_mast_ra$mast_re$analysis <- "rescueSim (333)"
resc_out_333_mast_ra_cdr$mast_re$analysis <- "rescueSim (333)"
resc_out_333_mast_ra_cdr_nn$mast_re$analysis <- "rescueSim (333)"
resc_out_333_deseq$deseq$analysis <- "rescueSim (333)"
resc_out_333_deseq_s$deseq$analysis <- "rescueSim (333)"

## rescueSim norm --------------------------------------------------------------
resc_out_norm_mast_cdr$mast_fe$analysis      <- "rescueSim (norm)"
resc_out_norm_mast_ra_cdr$mast_re$analysis   <- "rescueSim (norm)"
resc_out_norm_deseq_s$deseq$analysis         <- "rescueSim (norm)"

## rescueSim (123) norm --------------------------------------------------------
resc_out_123_norm_mast_cdr$mast_fe$analysis    <- "rescueSim (123) (norm)"
resc_out_123_norm_mast_re_cdr$mast_re$analysis <- "rescueSim (123) (norm)"
resc_out_123_norm_deseq_s$deseq$analysis       <- "rescueSim (123) (norm)"

## rescueSim (111) norm --------------------------------------------------------
resc_out_111_norm_mast_cdr$mast_fe$analysis    <- "rescueSim (111) (norm)"
resc_out_111_norm_mast_ra_cdr$mast_re$analysis <- "rescueSim (111) (norm)"
resc_out_111_norm_deseq_s$deseq$analysis       <- "rescueSim (111) (norm)"

## rescueSim (222) norm --------------------------------------------------------
resc_out_222_norm_mast_cdr$mast_fe$analysis    <- "rescueSim (222) (norm)"
resc_out_222_norm_mast_ra_cdr$mast_re$analysis <- "rescueSim (222) (norm)"
resc_out_222_norm_deseq_s$deseq$analysis       <- "rescueSim (222) (norm)"

## rescueSim (333) norm --------------------------------------------------------
resc_out_333_norm_mast_cdr$mast_fe$analysis    <- "rescueSim (333) (norm)"
resc_out_333_norm_mast_ra_cdr$mast_re$analysis <- "rescueSim (333) (norm)"
resc_out_333_norm_deseq_s$deseq$analysis       <- "rescueSim (333) (norm)"



method_data <- rbind(
  splatpop_out$mast_fe, splatpop_out$mast_re, splatpop_out$deseq,
  splatpop_out_no_cdr$mast_fe, splatpop_out_no_cdr$mast_re,
  splatpop_out_no_norm$mast_fe, splatpop_out_no_norm$mast_re,
  splatpop_out_11$mast_fe, splatpop_out_11$mast_re, splatpop_out_11$deseq,
  splatpop_out_22$mast_fe, splatpop_out_22$mast_re, splatpop_out_22$deseq,
  splatpop_out_33$mast_fe, splatpop_out_33$mast_re, splatpop_out_33$deseq,
  splatpop_out_44$mast_fe, splatpop_out_44$mast_re, splatpop_out_44$deseq,


  # resc_out$mast_fe, resc_out$mast_re, resc_out$deseq,
  resc_out_mast$mast_fe, resc_out_mast_cdr$mast_fe,
  resc_out_mast_cdr_nn$mast_fe,
  resc_out_mast_ra$mast_re,resc_out_mast_ra_cdr$mast_re,
  resc_out_deseq$deseq, resc_out_deseq_s$deseq,

# resc_out_123$mast_fe, resc_out_123$mast_re, resc_out_123$deseq,
  resc_out_123_mast$mast_fe, resc_out_123_mast_cdr$mast_fe,
resc_out_123_mast_cdr_nn$mast_fe,
resc_out_123_mast_re$mast_re,resc_out_123_mast_re_cdr$mast_re,
resc_out_123_deseq$deseq, resc_out_123_deseq_s$deseq,

  # resc_out_111$mast_fe, resc_out_111$mast_re, resc_out_111$deseq,
  resc_out_111_mast$mast_fe, resc_out_111_mast_cdr$mast_fe,
resc_out_111_mast_cdr_nn$mast_fe,
resc_out_111_mast_ra$mast_re, resc_out_111_mast_ra_cdr$mast_re,
# resc_out_111_mast_ra_cdr_nn$mast_re,
resc_out_111_deseq$deseq, resc_out_111_deseq_s$deseq,

  # resc_out_222$mast_fe, resc_out_222$mast_re, resc_out_222$deseq,
  resc_out_222_mast$mast_fe, resc_out_222_mast_cdr$mast_fe,
resc_out_222_mast_cdr_nn$mast_fe,
resc_out_222_mast_ra$mast_re, resc_out_222_mast_ra_cdr$mast_re,
# resc_out_222_mast_ra_cdr_nn$mast_re,
resc_out_222_deseq$deseq, resc_out_222_deseq_s$deseq,

  # resc_out_333$mast_fe, resc_out_333$mast_re, resc_out_333$deseq,
  resc_out_333_mast$mast_fe, resc_out_333_mast_cdr$mast_fe,
resc_out_333_mast_cdr_nn$mast_fe,
resc_out_333_mast_ra$mast_re, resc_out_333_mast_ra_cdr$mast_re,
# resc_out_333_mast_ra_cdr_nn$mast_re,
resc_out_333_deseq$deseq, resc_out_333_deseq_s$deseq,

resc_out_norm_mast_cdr$mast_fe, resc_out_norm_mast_ra_cdr$mast_re,
resc_out_norm_deseq_s$deseq, resc_out_123_norm_mast_cdr$mast_fe,
resc_out_123_norm_mast_re_cdr$mast_re,resc_out_123_norm_deseq_s$deseq,
resc_out_111_norm_mast_cdr$mast_fe, resc_out_111_norm_mast_ra_cdr$mast_re,
resc_out_111_norm_deseq_s$deseq, resc_out_222_norm_mast_cdr$mast_fe,
resc_out_222_norm_mast_ra_cdr$mast_re, resc_out_222_norm_deseq_s$deseq,
resc_out_333_norm_mast_cdr$mast_fe, resc_out_333_norm_mast_ra_cdr$mast_re,
resc_out_333_norm_deseq_s$deseq,

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

ground_truth <- pins::pin_read(board, "ground_truth_final")

method_data <- method_data |>
  dplyr::mutate(
    analysis_join = analysis
  ) |>
  dplyr::mutate(
    # rescueSim (123) -> rescueSim_123
    analysis_join = gsub(
      pattern = "^rescueSim \\((\\d+)\\)$",
      replacement = "rescueSim_\\1",
      x = analysis_join
    ),
    # splatpop (11) -> splatpop_11
    analysis_join = gsub(
      pattern = "^splatpop \\((\\d+)\\)$",
      replacement = "splatpop_\\1",
      x = analysis_join
    ),
    # rescueSim (norm) -> rescueSim norm
    analysis_join = gsub(
      pattern = "^rescueSim \\(norm\\)$",
      replacement = "rescueSim norm",
      x = analysis_join
    ),
    # rescueSim (123) (norm) -> rescueSim_123 norm
    analysis_join = gsub(
      pattern = "^rescueSim \\((\\d+)\\) \\(norm\\)$",
      replacement = "rescueSim_\\1 norm",
      x = analysis_join
    )
  )

compare_data <- dplyr::bind_rows(
  lapply(names(ground_truth), function(sim) {
    dplyr::left_join(
      method_data |> dplyr::filter(analysis_join == sim),
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

## rescuesim ----------------------------------------------------------------
# Create a pin for how many of the genes where in the conf intervals for the
# rescuesim datasets
compare_data <- pins::pin_read(board, "compare_data_final")
summarise_ci_coverage <- function(df, conf_level = 0.99) {
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

ci_level <- 0.99

compute_ci_summary <- function(compare_data, analyses, analysis_key = "analysis", ci_level = 0.99) {

  compare_data |>
    dplyr::filter(.data[[analysis_key]] %in% analyses) |>
    dplyr::group_by(.data[[analysis_key]], method, formula) |>
    dplyr::group_modify(~ summarise_ci_coverage(.x, conf_level = ci_level)) |>
    dplyr::ungroup() |>
    dplyr::rename(analysis = .data[[analysis_key]]) |>
    dplyr::mutate(
      conf_level = ci_level,
      prop_ci_total   = in_ci_total / total_genes,
      prop_ci_up      = dplyr::if_else(total_up > 0,
                                       in_ci_up / total_up,
                                       NA_real_),
      prop_ci_down    = dplyr::if_else(total_down > 0,
                                       in_ci_down / total_down,
                                       NA_real_),
      prop_ci_not_sig = dplyr::if_else(total_not_sig > 0,
                                       in_ci_not_sig / total_not_sig,
                                       NA_real_)
    )
}

# if you used the join-key names like "splatpop_11", "rescueSim_123 norm", etc.
splat_analyses <- c("splatpop", "splatpop_11", "splatpop_22", "splatpop_33",
                    "splatpop_44")

rescue_analyses <- c("rescueSim", "rescueSim_123", "rescueSim_111",
                     "rescueSim_222", "rescueSim_333")

rescue_norm_analyses <- c(
  "rescueSim norm",
  "rescueSim_123 norm",
  "rescueSim_111 norm",
  "rescueSim_222 norm",
  "rescueSim_333 norm"
)

splat_ci_summary <- compute_ci_summary(compare_data,
                                             splat_analyses,
                                             analysis_key = "analysis_join",
                                             ci_level = ci_level)
splat_ci_summary$sim <- "splatpop"

rescue_ci_summary <- compute_ci_summary(compare_data, rescue_analyses,
                                             analysis_key = "analysis_join",
                                             ci_level = ci_level)
rescue_ci_summary$sim <- "rescueSim"
rescue_norm_ci_summary <- compute_ci_summary(compare_data,
                                             rescue_norm_analyses,
                                             analysis_key = "analysis_join",
                                             ci_level = ci_level)


rescue_norm_ci_summary$sim <- "rescuesim norm"

## merge dataframes ------------------------------------------------------------
ci_data <- rbind(splat_ci_summary, rescue_ci_summary, rescue_norm_ci_summary)

board <- pins::board_local("nrgn-data")

pins::pin_write(
  ci_data,
  name = "rescueSim_ci_real",
  type = "rds",
  board = board,
  versioned = TRUE
)



