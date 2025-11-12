# Simulation functions----------------------------------------------------------
## TRT -------------------------------------------------------------------------
### SplatPop -------------------------------------------------------------------
#### SplatPop params -----------------------------------------------------------
make_splatpop_params <- function(eqtl.n = 0.5,
                                 batchCells = 100,
                                 similarity.scale = 6,
                                 condition.prob = c(0.5, 0.5),
                                 eqtl.condition.specific = 0.5,
                                 cde.facLoc = 1,
                                 cde.facScale = 1,
                                 dropout.type = "experiment",
                                 dropout.shape = -1,
                                 dropout.mid = 5) {
  p <- splatter::newSplatPopParams(
    eqtl.n = eqtl.n,
    # Proportion of genes with genetic effects (eQTLs); 50% of genes affected
    batchCells = batchCells,
    # Number of cells per simulated batch/group (here 100)
    similarity.scale = similarity.scale,
    # How similar individuals/cells are; higher = more similar
    condition.prob = condition.prob,
    # Probability of assigning cells to each condition (50% Placebo,
    # 50% Treatment)
    eqtl.condition.specific = eqtl.condition.specific,
    # Fraction of eQTLs that are condition-specific (50% only in one condition)
    cde.facLoc = cde.facLoc,
    # Location factor for condition-specific differential expression
    # (average log2FC)
    cde.facScale = cde.facScale
    # Scale factor for variability of condition effects (log2FC variance)
  )
  # Set dropout type and parameters to simulate realistic zero-inflation
  p <- splatter::setParam(p, "dropout.type", dropout.type)
  # Dropout depends on gene expression like in real experiments
  p <- splatter::setParam(p, "dropout.shape", dropout.shape)
  # Shape of dropout curve; negative = more zeros in low-expression genes
  p <- splatter::setParam(p, "dropout.mid", dropout.mid)
  # Midpoint of dropout curve;
  # expression ~5 (log-scale) has 50% dropout probability
  print("params")
  p
}

#### SplatPop data -------------------------------------------------------------
simulate_splatpop <- function(n_samples = 20,
                              n_snps = 1000,
                              params = NULL,
                              sparsify = TRUE,
                              seed = NULL) {
  set.seed(42)
  if (is.null(params)) params <- make_splatpop_params()
  vcf <- splatter::mockVCF(n.snps = n_snps, n.samples = n_samples)
  gff <- splatter::mockGFF()
  sim <- splatter::splatPopSimulate(
    vcf = vcf,
    gff = gff,
    params = params,
    sparsify = sparsify # to make it more realisic, for dropout
  )
  print("simulated")
  sim <- scater::logNormCounts(sim)
  sim <- scater::runPCA(sim, ncomponents = 10)
  sim <- scater::runUMAP(sim, dimred = "PCA")
  print("dim red")
  ## rename conditions for clarity
  SummarizedExperiment::colData(sim)$Condition <- factor(
    SummarizedExperiment::colData(sim)$Condition,
    levels = c("Condition1", "Condition2"),
    labels = c("Placebo", "Treatment")
  )
  sim
}

## TRT + time ------------------------------------------------------------------
### rescueSim ------------------------------------------------------------------
#### rescueSim params ----------------------------------------------------------
make_rescuesim_params <- function(RecAM_sce,
                                  nTimepoints = 2,
                                  twoGroupDesign = TRUE,
                                  nSubjsPerGroup = 5,
                                  propDE = 0.20,
                                  deLog2FC = 1,
                                  sampleVariable = "sampleID",
                                  subjectVariable = "subjectID",
                                  timepointVariable = "time") {
  p <- rescueSim::RescueSimParams()
  p <- rescueSim::estRescueSimParams(RecAM_sce,
    sampleVariable = sampleVariable,
    subjectVariable = subjectVariable,
    timepointVariable = timepointVariable,
    groupVariable = NULL,
    nonDEGs = NULL,
    cellParamsByCondition = FALSE
  )

  ## For the 2nd time point, simulate +/- 1 log2FC in 20% of genes
  group0_log2FC<-sample(c(-1, 0, 1),
                       nrow(RecAM_sce),
                       replace=T,
                       prob = c(.1, .8, .1))
  ## For final time point, multiply log2FC by 3 to get +/- 3 log2FC in 20% of genes
  group1_log2FC <- 3*group0_log2FC
  p@nTimepoints <- nTimepoints
  p@twoGroupDesign <- twoGroupDesign
  p@nSubjsPerGroup <- nSubjsPerGroup
  p@propDE <- propDE
  p@deLog2FC <- list(time1_group0= group0_log2FC, time1_group1 = group1_log2FC)
  p
}
#### rescueSim data ------------------------------------------------------------
simulate_rescue <- function(params, seed = NULL) {
  set.seed(42)
  sim <- rescueSim::simRescueData(params)
  sim <- scater::logNormCounts(sim)
  sim <- scater::runPCA(sim, ncomponents = 10)
  sim <- scater::runUMAP(sim, dimred = "PCA")
  # harmonize colData
  SummarizedExperiment::colData(sim)$Sample <- SummarizedExperiment::colData(
    sim
  )$SampleID
  SummarizedExperiment::colData(sim)$group <- factor(
    SummarizedExperiment::colData(sim)$group,
    levels = c("group0", "group1"),
    labels = c("Placebo", "Treatment")
  )
  SummarizedExperiment::colData(sim)$time <- factor(
    SummarizedExperiment::colData(sim)$time,
    levels = c("time0", "time1"),
    labels = c("Baseline", "EoT")
  )
  SummarizedExperiment::colData(sim)$Condition <- paste0(
    SummarizedExperiment::colData(sim)$group, "_",
    SummarizedExperiment::colData(sim)$time
  )
  SummarizedExperiment::colData(sim)$sampleID <- factor(
    SummarizedExperiment::colData(sim)$sampleID
  )
  SummarizedExperiment::colData(sim)$subjectID <- factor(
    SummarizedExperiment::colData(sim)$subjectID
  )
  SummarizedExperiment::colData(sim)$CellID <- factor(
    SummarizedExperiment::rownames(SummarizedExperiment::colData(sim))
  )
  ground_truth <- SummarizedExperiment::rowData(sim)
  simulation <-  list(
    sce = sim,
    ground_truth = ground_truth)

  simulation
}

## Jan -------------------------------------------------------------------------
simulate_jan <- function(nSubjsPerGroup = 20,
                         nCellsPerSubjectTimepoint = 100,
                         nGenes = 100) {
  set.seed(42)

  # Jan-style simulation setup
  n_subjects <- nSubjsPerGroup * 2
  n_timepoints <- 2

  # Experimental design
  design_metadata <- expand.grid(
    Subject = 1:n_subjects,
    Time = 0:(n_timepoints - 1)
  ) |>
    dplyr::slice(rep(1:dplyr::n(), each = nCellsPerSubjectTimepoint)) |>
    dplyr::mutate(
      Treatment = ifelse(Subject <= nSubjsPerGroup, 0, 1),
      CellID = paste0("Cell_", 1:dplyr::n())
    )

  # Splatter base counts
  params <- splatter::newSplatParams()
  params <- splatter::setParam(params, "nGenes", nGenes)
  params <- splatter::setParam(params, "batchCells", nrow(design_metadata))
  params <- splatter::setParam(params, "group.prob", 1)
  params <- splatter::setParam(params, "de.prob", 0)
  params <- splatter::setParam(params, "dropout.type", "experiment")

  sce <- splatter::splatSimulate(params, verbose = FALSE)

  # Add effects manually (like Jan example)
  counts_matrix <- SingleCellExperiment::counts(sce)
  log_counts <- log2(counts_matrix + 1)

  # Placebo and treatment effects
  placebo_genes <- 1:20
  treatment_genes <- 15:25
  both_effects_genes <- intersect(placebo_genes, treatment_genes)

  for (gene_idx in 1:nGenes) {
    for (cell_idx in 1:nrow(design_metadata)) {
      subject <- design_metadata$Subject[cell_idx]
      time <- design_metadata$Time[cell_idx]
      treatment <- design_metadata$Treatment[cell_idx]
      base_expr <- log_counts[gene_idx, cell_idx]
      if (gene_idx %in% placebo_genes) base_expr <- base_expr + time * 1.0
      if (gene_idx %in% treatment_genes) {
        base_expr <- base_expr +
          time * treatment * 1.5
      }
      new_count <- round(2^base_expr - 1)
      counts_matrix[gene_idx, cell_idx] <- ifelse(new_count < 2, 0, new_count)
    }
  }

  # Update SCE
  SingleCellExperiment::counts(sce) <- as(counts_matrix, "dgCMatrix")
  SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce)
  + 1)
  SummarizedExperiment::colData(sce) <- cbind(
    SummarizedExperiment::colData(
      sce
    ),
    design_metadata
  )
  colnames(sce) <- design_metadata$CellID
  rownames(sce) <- paste0("Gene_", 1:nGenes)

  # Factor columns for wrappers
  sce$CellID <- factor(sce$CellID)
  sce$SubjectID <- factor(sce$Subject)
  sce$group <- factor(ifelse(sce$Treatment == 0, "Placebo", "Treatment"))
  sce$time <- factor(ifelse(sce$Time == 0, "Baseline", "EoT"))
  sce$Condition <- paste0(
    sce$group, "_", sce$time
  )

  sim_jan <- sce

  sim_jan <- scater::runPCA(sim_jan, ncomponents = 10)
  sim_jan <- scater::runUMAP(sim_jan, dimred = "PCA")
  ground_truth <- data.frame(
    Gene = paste0("Gene_", 1:nGenes),
    PlaceboEffect = 1:nGenes %in% placebo_genes,
    TreatmentEffect = 1:nGenes %in% treatment_genes,
    BothEffects = 1:nGenes %in% both_effects_genes
  )

  simulation <-  list(
    sce = sim_jan,
    ground_truth = ground_truth)

  simulation
}




# Methods ----------------------------------------------------------------------
## MAST ---------------------- -------------------------------------------------
run_mast <- function(
  sce,
  formula = ~Condition,
  lrt = "ConditionTreatment",
  random_var = NULL,
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
    sample_var <- random_var
    cd <- SummarizedExperiment::colData(sca)
    samples <- cd[[sample_var]]

    # --- 1. basic per-sample detection ---
    sample_list <- split(seq_len(ncol(expr)), samples)
    detected_in_sample <- vapply(
      sample_list,
      function(idx) rowSums(expr[, idx, drop = FALSE] > 0) > 0,
      logical(nrow(expr))
    )
    if (is.vector(detected_in_sample)) detected_in_sample <- matrix(
      detected_in_sample,
      ncol = 1
      )
    min_samples <- 5
    keep_overall <- rowSums(detected_in_sample) >= min_samples

    # --- 2. require some expression in every group×time combo ---
    if (all(c("group", "time") %in% colnames(cd))) {
      combo <- interaction(cd$group, cd$time, drop = TRUE)
    } else if ("Condition" %in% colnames(cd)) {
      combo <- cd$Condition
    } else {
      combo <- factor(rep("All", ncol(expr)))
      message("No 'group' or 'time' columns found; using all cells as one combo.")
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

    # --- 3. NEW: ensure enough non-zero observations for random effects ---
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
  zlmfit <- tryCatch(
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
  # Step 1: Aggregate to pseudobulk
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

  # Step 2: Prepare metadata
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

  # Step 3: Create DESeq dataset
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

  # Step 4: Filter low-count genes
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

  # Step 5: Fit model
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

  # Step 6: Extract results
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

  # Step 7: Clean and return results
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
  design_formula_mast = ~Condition,
  design_formula_mast_rand = ~ Condition + (1 | Sample),
  design_formula_deseq = ~Condition,
  lrt_term = "ConditionTreatment", # term to test in MAST
  ids_vars = c("Sample", "Condition"), # columns for pseudobulk aggregation
  run_mast_fe = TRUE,
  run_mast_re = TRUE,
  run_deseq = TRUE
) {
  res_list <- list(simulation = sce)

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
      res_list$deseq <- deseq_res
    }
  }

  message("Finished running selected methods.")
  return(res_list)
}


# Pipelines --------------------------------------------------------------------
## SplatPop pipeline with methods ----------------------------------------------
pipeline_splatpop_and_methods <- function(n_samples = 20,
                                          n_snps = 1000,
                                          params = NULL,
                                          run_mast_fe = TRUE,
                                          run_mast_re = TRUE,
                                          run_deseq = TRUE) {
  set.seed(42)
  if (is.null(params)) params <- make_splatpop_params()
  sim <- simulate_splatpop(
    n_samples = n_samples,
    n_snps = n_snps,
    params = params
  )
  splat_results <- run_methods(
    sce = sim,
    condition_var = "Condition",
    sample_var = "Sample",
    design_formula_mast = ~Condition,
    design_formula_mast_rand = ~ Condition + (1 | Sample),
    design_formula_deseq = ~Condition,
    lrt_term = "ConditionTreatment",
    ids_vars = c("Sample", "Condition")
  )
}

## Jan pipeline with methods ---------------------------------------------------
pipeline_jan_and_methods <- function(nSubjsPerGroup = 20,
                                     nCellsPerSubjectTimepoint = 100,
                                     nGenes = 100,
                                     run_mast_fe = TRUE,
                                     run_mast_re = TRUE,
                                     run_deseq = TRUE) {
  data <- simulate_jan(
    nSubjsPerGroup = nSubjsPerGroup,
    nCellsPerSubjectTimepoint = nCellsPerSubjectTimepoint,
    nGenes = nGenes
  )
  sim_jan <- data$sce


  ground_truth <- data$ground_truth
  jan_results <- run_methods(
    sce = sim_jan,
    sample_var = "SubjectID",
    design_formula = ~ group * time,
    design_formula_mast_rand = ~ group * time + (1 | SubjectID),
    design_formula_deseq = ~ group * time,
    lrt_term = "groupTreatment:timeEoT",
    ids_vars = c("SubjectID", "group", "time"),
    run_mast_fe = run_mast_fe,
    run_mast_re = run_mast_re,
    run_deseq = run_deseq
  )
  jan_results$ground_truth <- ground_truth
  jan_results
}

## rescueSim pipeline with methods ---------------------------------------------
pipeline_rescue_and_methods <- function(RecAM_sce,
                                        nSubjsPerGroup = 5,
                                        run_mast = TRUE,
                                        run_mast_ran = TRUE,
                                        run_deseq = TRUE) {
  params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = nSubjsPerGroup)
  data <- simulate_rescue(params)
  sim <- data$sce
  ground_truth <- data$ground_truth

  rescue_results <- run_methods(
    sce = sim,
    sample_var = "sampleID",
    design_formula_mast = ~ subjectID + group * time,
    design_formula_mast_rand = ~  group * time + (1 | sampleID),
    design_formula_deseq = ~ group * time,
    lrt_term = "groupTreatment:timeEoT",
    ids_vars = c("subjectID", "group", "time"),
    run_mast_fe = run_mast,
    run_mast_re = run_mast_ran,
    run_deseq = run_deseq
  )
  rescue_results$ground_truth <- ground_truth
  rescue_results
}





# Run Simulations  -------------------------------------------------------------

## SplatPop --------------------------------------------------------------------
splatpop_out <- pipeline_splatpop_and_methods()
## Jan -------------------------------------------------------------------------
jan_out <- pipeline_jan_and_methods()
##  rescueSim ------------------------------------------------------------------
data("RecAM_sce", package = "rescueSim")
resc_out <- pipeline_rescue_and_methods(RecAM_sce)




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
simulation <- c(
  "Treatment" = splatpop_out$simulation,
  "Treatment + Time (rescueSim)" = resc_out$simulation,
  "Treatment + Time (splatter)" = jan_out$simulation
)

ground_truth <- simulation <- c(
  # "Treatment" = splatpop_out$simulation,
  "Treatment + Time (rescueSim)" = resc_out$ground_truth,
  "Treatment + Time (splatter)" = jan_out$ground_truth
)

board <- pins::board_local("nrgn-data")

pins::pin_write(simulation,
  name = "simulations",
  type = "rds",
  board = board,
  versioned = FALSE
)

pins::pin_write(method_data,
  name = "method_data",
  type = "rds",
  board = board,
  versioned = FALSE
)

pins::pin_write(ground_truth,
                name = "ground_truth",
                type = "rds",
                board = board,
                versioned = FALSE
)
