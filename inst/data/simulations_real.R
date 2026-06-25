# Simulation functions ----------------------------------------------------------
## Fair rescueSim benchmark -----------------------------------------------------

# ------------------------------------------------------------------------------
# Fair rescueSim benchmark simulations
#
# Memory-safe version
#
# What this does:
#   - creates 30 simulations per scenario and nSubjsPerGroup
#   - saves each SingleCellExperiment as an individual .rds file
#   - saves each ground truth table as an individual .rds file
#   - pins only the small metadata/path tables
#
# Total simulations:
#   5 scenarios * 3 n values * 30 seeds = 450 simulations
#
# IMPORTANT:
#   This file does NOT create fair_rescue_simulations <- list()
#   because that is what was causing the memory crash.
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(S4Vectors)
library(scater)
library(pins)

# ------------------------------------------------------------------------------
# Create fair rescueSim parameters
# ------------------------------------------------------------------------------

make_rescuesim_fair_params <- function(
    RecAM_sce,
    nTimepoints = 2,
    twoGroupDesign = TRUE,
    nSubjsPerGroup = 5,
    prop_time = 0.50,
    prop_trt = 0.20,
    time_sd = 0.30,
    trt_lfc_mean = 0.60,
    trt_lfc_sd = 0.15,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  p <- rescueSim::RescueSimParams()

  p <- rescueSim::estRescueSimParams(
    RecAM_sce,
    sampleVariable = "sampleID",
    subjectVariable = "subjectID",
    timepointVariable = "time",
    groupVariable = NULL,
    nonDEGs = NULL,
    cellParamsByCondition = FALSE
  )

  G <- nrow(RecAM_sce)

  # --------------------------------------------------------------------------
  # 1) Shared time effect
  #
  # This affects placebo and treatment equally.
  # It should NOT be detected as treatment effect.
  # --------------------------------------------------------------------------

  time_effect <- rep(0, G)

  if (prop_time > 0 && time_sd > 0) {
    n_time <- round(prop_time * G)

    time_idx <- sample.int(
      n = G,
      size = n_time,
      replace = FALSE
    )

    time_effect[time_idx] <- rnorm(
      n = n_time,
      mean = 0,
      sd = time_sd
    )
  }

  # --------------------------------------------------------------------------
  # 2) Extra treatment-by-time effect
  #
  # This is the true DE signal that methods should detect.
  # --------------------------------------------------------------------------

  trt_extra <- rep(0, G)

  if (prop_trt > 0 && trt_lfc_mean > 0) {
    n_trt <- round(prop_trt * G)

    trt_idx <- sample.int(
      n = G,
      size = n_trt,
      replace = FALSE
    )

    trt_sign <- sample(
      x = c(-1, 1),
      size = n_trt,
      replace = TRUE
    )

    trt_magnitude <- abs(rnorm(
      n = n_trt,
      mean = trt_lfc_mean,
      sd = trt_lfc_sd
    ))

    trt_extra[trt_idx] <- trt_sign * trt_magnitude
  }

  # Placebo EoT - Placebo Baseline
  group0_log2FC <- time_effect

  # Treatment EoT - Treatment Baseline
  group1_log2FC <- time_effect + trt_extra

  p@nTimepoints <- nTimepoints
  p@twoGroupDesign <- twoGroupDesign
  p@nSubjsPerGroup <- nSubjsPerGroup
  p@propDE <- prop_trt

  p@deLog2FC <- list(
    time1_group0 = group0_log2FC,
    time1_group1 = group1_log2FC
  )

  p
}

# ------------------------------------------------------------------------------
# Simulate one fair rescueSim dataset
# ------------------------------------------------------------------------------

simulate_rescue_fair <- function(
    params,
    scenario_name,
    scenario,
    nSubjsPerGroup,
    seed,
    compute_dimred = FALSE
) {
  if (!is.null(seed)) set.seed(seed)

  sim <- rescueSim::simRescueData(params)

  # For benchmarking methods, PCA/UMAP are not needed.
  # Keeping this FALSE saves memory and disk space.
  if (compute_dimred) {
    sim <- scater::logNormCounts(sim)
    sim <- scater::runPCA(sim, ncomponents = 10)
    sim <- scater::runUMAP(sim, dimred = "PCA")
  }

  cd <- as.data.frame(SummarizedExperiment::colData(sim))

  # rescueSim sometimes has both SampleID and sampleID depending on object version.
  if ("SampleID" %in% colnames(cd)) {
    cd$Sample <- cd$SampleID
  } else if ("sampleID" %in% colnames(cd)) {
    cd$Sample <- cd$sampleID
  }

  cd$group <- factor(
    cd$group,
    levels = c("group0", "group1"),
    labels = c("Placebo", "Treatment")
  )

  cd$time <- factor(
    cd$time,
    levels = c("time0", "time1"),
    labels = c("Baseline", "EoT")
  )

  cd$group <- stats::relevel(factor(cd$group), ref = "Placebo")
  cd$time <- stats::relevel(factor(cd$time), ref = "Baseline")

  cd$Condition <- paste0(cd$group, "_", cd$time)

  cd$sampleID <- factor(cd$sampleID)
  cd$subjectID <- factor(cd$subjectID)
  cd$CellID <- factor(rownames(cd))

  # Useful for the fixed paired design later
  cd$time_num <- as.integer(cd$time == "EoT")
  cd$treat_num <- as.integer(cd$group == "Treatment")
  cd$trt_time <- cd$time_num * cd$treat_num

  SummarizedExperiment::colData(sim) <- S4Vectors::DataFrame(cd)

  # --------------------------------------------------------------------------
  # Ground truth
  # --------------------------------------------------------------------------

  g0 <- SummarizedExperiment::rowData(sim)$deLog2FC.time1_group0
  g1 <- SummarizedExperiment::rowData(sim)$deLog2FC.time1_group1

  true_logFC <- g1 - g0

  tol <- 1e-12

  ground_truth <- data.frame(
    gene = rownames(sim),
    true_logFC = true_logFC,
    trueDE = abs(true_logFC) > tol,
    placebo_effect = g0,
    treatment_effect = g1,
    scenario = scenario,
    simulation = scenario_name,
    nSubjsPerGroup = nSubjsPerGroup,
    seed = seed,
    effect_type = dplyr::case_when(
      abs(true_logFC) > tol ~ "treatment_time",
      abs(true_logFC) <= tol & abs(g0) > tol ~ "time_only",
      TRUE ~ "null"
    ),
    direction = dplyr::case_when(
      true_logFC > tol ~ "up",
      true_logFC < -tol ~ "down",
      TRUE ~ "none"
    )
  )

  list(
    sce = sim,
    ground_truth = ground_truth
  )
}

# ------------------------------------------------------------------------------
# Benchmark grid
# ------------------------------------------------------------------------------

n_sims_per_setting <- 30
base_seed <- 202600

fair_rescue_grid <- tidyr::crossing(
  scenario = c(
    "global_null",
    "time_only_null",
    "weak",
    "medium",
    "strong"
  ),
  nSubjsPerGroup = c(5, 10, 20),
  seed = base_seed + seq_len(n_sims_per_setting)
) |>
  dplyr::mutate(
    prop_time = dplyr::case_when(
      scenario == "global_null" ~ 0.00,
      TRUE ~ 0.50
    ),

    prop_trt = dplyr::case_when(
      scenario %in% c("global_null", "time_only_null") ~ 0.00,
      TRUE ~ 0.20
    ),

    time_sd = dplyr::case_when(
      scenario == "global_null" ~ 0.00,
      TRUE ~ 0.30
    ),

    trt_lfc_mean = dplyr::case_when(
      scenario %in% c("global_null", "time_only_null") ~ 0.00,
      scenario == "weak" ~ 0.30,
      scenario == "medium" ~ 0.60,
      scenario == "strong" ~ 1.00
    ),

    trt_lfc_sd = dplyr::case_when(
      scenario %in% c("global_null", "time_only_null") ~ 0.00,
      scenario == "weak" ~ 0.10,
      scenario == "medium" ~ 0.15,
      scenario == "strong" ~ 0.25
    ),

    sim_name = paste0(
      "fair_",
      scenario,
      "_n",
      nSubjsPerGroup,
      "_seed",
      seed
    )
  )

# Check that the grid is correct
print(
  fair_rescue_grid |>
    dplyr::count(scenario, nSubjsPerGroup)
)

print(nrow(fair_rescue_grid))
# Should be 450

# ------------------------------------------------------------------------------
# Run fair rescueSim simulations safely
#
# This is the key part.
# It saves one simulation, removes it from memory, then continues.
# ------------------------------------------------------------------------------

data("RecAM_sce", package = "rescueSim")

board <- pins::board_local("nrgn-data")

base_dir <- "nrgn-data"

sim_dir <- file.path(base_dir, "fair_rescue_simulation_files")
truth_dir <- file.path(base_dir, "fair_rescue_ground_truth_files")

dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)

fair_rescue_sim_paths <- list()

# If FALSE, already-created files are skipped.
# This lets you resume after a crash or after stopping the script.
force_rerun <- FALSE

for (i in seq_len(nrow(fair_rescue_grid))) {
  row_i <- fair_rescue_grid[i, ]

  sim_file <- file.path(
    sim_dir,
    paste0(row_i$sim_name, ".rds")
  )

  truth_file <- file.path(
    truth_dir,
    paste0(row_i$sim_name, "_truth.rds")
  )

  message(
    "Running ",
    i,
    " / ",
    nrow(fair_rescue_grid),
    ": ",
    row_i$sim_name
  )

  # --------------------------------------------------------------------------
  # Skip if already saved
  # --------------------------------------------------------------------------

  if (!force_rerun && file.exists(sim_file) && file.exists(truth_file)) {
    message("Files already exist. Skipping: ", row_i$sim_name)

    fair_rescue_sim_paths[[row_i$sim_name]] <- data.frame(
      simulation = row_i$sim_name,
      file = sim_file,
      ground_truth_file = truth_file,
      scenario = row_i$scenario,
      nSubjsPerGroup = row_i$nSubjsPerGroup,
      seed = row_i$seed,
      prop_time = row_i$prop_time,
      prop_trt = row_i$prop_trt,
      time_sd = row_i$time_sd,
      trt_lfc_mean = row_i$trt_lfc_mean,
      trt_lfc_sd = row_i$trt_lfc_sd,
      stringsAsFactors = FALSE
    )

    next
  }

  # --------------------------------------------------------------------------
  # Simulate
  # --------------------------------------------------------------------------

  params_i <- make_rescuesim_fair_params(
    RecAM_sce = RecAM_sce,
    nSubjsPerGroup = row_i$nSubjsPerGroup,
    prop_time = row_i$prop_time,
    prop_trt = row_i$prop_trt,
    time_sd = row_i$time_sd,
    trt_lfc_mean = row_i$trt_lfc_mean,
    trt_lfc_sd = row_i$trt_lfc_sd,
    seed = row_i$seed
  )

  sim_i <- simulate_rescue_fair(
    params = params_i,
    scenario_name = row_i$sim_name,
    scenario = row_i$scenario,
    nSubjsPerGroup = row_i$nSubjsPerGroup,
    seed = row_i$seed,
    compute_dimred = FALSE
  )

  # --------------------------------------------------------------------------
  # Save to disk immediately
  # --------------------------------------------------------------------------

  saveRDS(
    sim_i$sce,
    file = sim_file,
    compress = "gzip"
  )

  saveRDS(
    sim_i$ground_truth,
    file = truth_file,
    compress = "gzip"
  )

  # --------------------------------------------------------------------------
  # Store only paths and metadata in memory
  # --------------------------------------------------------------------------

  fair_rescue_sim_paths[[row_i$sim_name]] <- data.frame(
    simulation = row_i$sim_name,
    file = sim_file,
    ground_truth_file = truth_file,
    scenario = row_i$scenario,
    nSubjsPerGroup = row_i$nSubjsPerGroup,
    seed = row_i$seed,
    prop_time = row_i$prop_time,
    prop_trt = row_i$prop_trt,
    time_sd = row_i$time_sd,
    trt_lfc_mean = row_i$trt_lfc_mean,
    trt_lfc_sd = row_i$trt_lfc_sd,
    stringsAsFactors = FALSE
  )

  # Remove heavy objects
  rm(params_i, sim_i)
  gc()
}

fair_rescue_sim_paths <- dplyr::bind_rows(fair_rescue_sim_paths)

# ------------------------------------------------------------------------------
# Ground truth summary
#
# This reads one ground truth file at a time.
# It does not load all SCEs.
# ------------------------------------------------------------------------------

fair_truth_summary <- dplyr::bind_rows(
  lapply(seq_len(nrow(fair_rescue_sim_paths)), function(i) {
    truth_i <- readRDS(fair_rescue_sim_paths$ground_truth_file[i])

    out <- truth_i |>
      dplyr::summarise(
        simulation = unique(simulation),
        scenario = unique(scenario),
        nSubjsPerGroup = unique(nSubjsPerGroup),
        seed = unique(seed),
        n_genes = dplyr::n(),
        n_true_de = sum(trueDE),
        prop_true_de = mean(trueDE),
        n_treatment_time = sum(effect_type == "treatment_time"),
        n_time_only = sum(effect_type == "time_only"),
        n_null = sum(effect_type == "null"),
        median_abs_true_lfc = median(abs(true_logFC), na.rm = TRUE),
        max_abs_true_lfc = max(abs(true_logFC), na.rm = TRUE),
        .groups = "drop"
      )

    rm(truth_i)
    gc()

    out
  })
)

print(
  fair_truth_summary |>
    dplyr::arrange(scenario, nSubjsPerGroup, seed)
)

# Check: should be 30 for every scenario and n
print(
  fair_truth_summary |>
    dplyr::count(scenario, nSubjsPerGroup)
)

# ------------------------------------------------------------------------------
# Save fair benchmark metadata
#
# Important:
#   We pin paths, not the full list of SingleCellExperiment objects.
# ------------------------------------------------------------------------------

pins::pin_write(
  fair_rescue_sim_paths,
  name = "fair_rescue_sim_paths",
  type = "rds",
  board = board,
  versioned = TRUE
)

pins::pin_write(
  fair_rescue_grid,
  name = "fair_rescue_grid",
  type = "rds",
  board = board,
  versioned = TRUE
)

pins::pin_write(
  fair_truth_summary,
  name = "fair_rescue_truth_summary",
  type = "rds",
  board = board,
  versioned = TRUE
)

# ------------------------------------------------------------------------------
# Optional helper: read one simulation later
# ------------------------------------------------------------------------------

read_fair_simulation <- function(simulation_name, board = pins::board_local("nrgn-data")) {
  paths <- pins::pin_read(board, "fair_rescue_sim_paths")

  row_i <- paths |>
    dplyr::filter(simulation == simulation_name)

  if (nrow(row_i) != 1) {
    stop("Could not find exactly one simulation called: ", simulation_name)
  }

  readRDS(row_i$file)
}

read_fair_ground_truth <- function(simulation_name, board = pins::board_local("nrgn-data")) {
  paths <- pins::pin_read(board, "fair_rescue_sim_paths")

  row_i <- paths |>
    dplyr::filter(simulation == simulation_name)

  if (nrow(row_i) != 1) {
    stop("Could not find exactly one simulation called: ", simulation_name)
  }

  readRDS(row_i$ground_truth_file)
}

# ------------------------------------------------------------------------------
# End
# ------------------------------------------------------------------------------
