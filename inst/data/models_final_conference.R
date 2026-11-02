# Simulation functions ----------------------------------------------------------
## Fair rescueSim benchmark -----------------------------------------------------

# ------------------------------------------------------------------------------
# Fair rescueSim benchmark simulations
#
# Memory-safe version:
#   - 100 seeds per scenario and nSubjsPerGroup
#   - saves each SingleCellExperiment to disk immediately
#   - does NOT keep all 1500 simulations in RAM
#
# Total simulations:
#   5 scenarios * 3 n values * 100 seeds = 1500 simulations
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
  # This is the real treatment effect that methods should detect.
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

  # For benchmarking methods, PCA/UMAP is not needed.
  # Keeping this FALSE saves a lot of memory.
  if (compute_dimred) {
    sim <- scater::logNormCounts(sim)
    sim <- scater::runPCA(sim, ncomponents = 10)
    sim <- scater::runUMAP(sim, dimred = "PCA")
  }

  cd <- as.data.frame(SummarizedExperiment::colData(sim))

  cd$Sample <- cd$SampleID

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

  # Useful for fixed paired designs later
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
#
# This creates:
#   100 seeds per scenario and per nSubjsPerGroup
# ------------------------------------------------------------------------------

n_sims_per_setting <- 100
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

# Check grid size
print(
  fair_rescue_grid |>
    dplyr::count(scenario, nSubjsPerGroup)
)

print(nrow(fair_rescue_grid))
# Should be 1500

# ------------------------------------------------------------------------------
# Run fair rescueSim simulations safely
#
# Important:
#   We save each SCE immediately as an .rds file.
#   We do NOT store all simulations in fair_rescue_simulations.
# ------------------------------------------------------------------------------

data("RecAM_sce", package = "rescueSim")

board <- pins::board_local("nrgn-data")

sim_dir <- file.path("nrgn-data", "fair_rescue_simulation_files")
truth_dir <- file.path("nrgn-data", "fair_rescue_ground_truth_files")

dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(truth_dir, recursive = TRUE, showWarnings = FALSE)

fair_rescue_sim_paths <- list()

# Set to TRUE if you want to overwrite already created files
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

  # If the files already exist, skip simulation.
  # This lets you resume if R crashes or you stop the script.
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

  rm(params_i, sim_i)
  gc()
}

fair_rescue_sim_paths <- dplyr::bind_rows(fair_rescue_sim_paths)

# ------------------------------------------------------------------------------
# Read ground truth files into one list
# ------------------------------------------------------------------------------

fair_rescue_ground_truth <- stats::setNames(
  lapply(fair_rescue_sim_paths$ground_truth_file, readRDS),
  fair_rescue_sim_paths$simulation
)

# ------------------------------------------------------------------------------
# Ground truth summary
# ------------------------------------------------------------------------------

fair_truth_summary <- dplyr::bind_rows(
  lapply(names(fair_rescue_ground_truth), function(nm) {
    fair_rescue_ground_truth[[nm]] |>
      dplyr::summarise(
        simulation = nm,
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
  })
)

print(
  fair_truth_summary |>
    dplyr::arrange(scenario, nSubjsPerGroup, seed)
)

# Check that we really have 100 simulations per scenario and n
print(
  fair_truth_summary |>
    dplyr::count(scenario, nSubjsPerGroup)
)

# ------------------------------------------------------------------------------
# Save fair benchmark metadata
#
# Important:
#   We pin the file paths, not the full SCE list.
# ------------------------------------------------------------------------------

pins::pin_write(
  fair_rescue_sim_paths,
  name = "fair_rescue_sim_paths",
  type = "rds",
  board = board,
  versioned = TRUE
)

pins::pin_write(
  fair_rescue_ground_truth,
  name = "fair_rescue_ground_truth",
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
# End
# ------------------------------------------------------------------------------
