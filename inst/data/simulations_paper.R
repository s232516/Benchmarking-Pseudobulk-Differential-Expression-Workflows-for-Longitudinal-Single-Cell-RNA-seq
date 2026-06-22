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
  p
}

#### SplatPop data -------------------------------------------------------------
simulate_splatpop <- function(n_samples = 20,
                              n_snps = 1000,
                              params = NULL,
                              sparsify = TRUE) {
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
  sim <- scater::logNormCounts(sim)
  sim <- scater::runPCA(sim, ncomponents = 10)
  sim <- scater::runUMAP(sim, dimred = "PCA")
  ## rename conditions for clarity
  SummarizedExperiment::colData(sim)$Condition <- factor(
    SummarizedExperiment::colData(sim)$Condition,
    levels = c("Condition1", "Condition2"),
    labels = c("Placebo", "Treatment")
  )
  ## extract ground truth
  real_data <- SummarizedExperiment::rowData(sim)

  cond1 <- real_data$ConditionDE.Condition1
  cond2 <- real_data$ConditionDE.Condition2

  # Some SplatPop runs store these as numeric; others as characters
  cond1 <- as.numeric(cond1)
  cond2 <- as.numeric(cond2)

  true_logFC <- log2(cond2 / cond1)

  ground_truth <- data.frame(
    gene = SummarizedExperiment::rownames(sim),
    trueDE = true_logFC != 0,
    true_logFC = true_logFC,
    source = "splatPop"
  )
  simulation <- list(
    sce = sim,
    ground_truth = ground_truth,
    simulation = "splatPop"
  )

  simulation
}

## TRT + time ------------------------------------------------------------------
### rescueSim ------------------------------------------------------------------
#### rescueSim params ----------------------------------------------------------
make_rescuesim_params <- function(RecAM_sce,
                                  nTimepoints = 2,
                                  twoGroupDesign = TRUE,
                                  nSubjsPerGroup = 5,
                                  propDE = 0.20) {
  p <- rescueSim::RescueSimParams()
  p <- rescueSim::estRescueSimParams(RecAM_sce,
    sampleVariable = "sampleID",
    subjectVariable = "subjectID",
    timepointVariable = "time",
    groupVariable = NULL,
    nonDEGs = NULL,
    cellParamsByCondition = FALSE
  )

  ## For the 2nd time point, simulate +/- 1 log2FC in 20% of genes
  group0_log2FC <- sample(c(-1, 0, 1),
    nrow(RecAM_sce),
    replace = T,
    prob = c(.1, .8, .1)
  )
  # For final time point, multiply log2FC by 3 to get +/- 3 log2FC in 20% of
  # genes
  group1_log2FC <- 3 * group0_log2FC
  p@nTimepoints <- nTimepoints
  p@twoGroupDesign <- twoGroupDesign
  p@nSubjsPerGroup <- nSubjsPerGroup
  p@propDE <- propDE
  p@deLog2FC <- list(time1_group0 = group0_log2FC, time1_group1 = group1_log2FC)
  p
}

# rescueSim norm lfc -----------------------------------------------------------
make_rescuesim_norm_params <- function(RecAM_sce,
                                       nTimepoints = 2,
                                       twoGroupDesign = TRUE,
                                       nSubjsPerGroup = 5,
                                       prop_time = 0.50, # % genes with time effect (placebo + treatment)
                                       prop_trt = 0.20, # % genes with EXTRA treatment effect
                                       sd_time = 0.40, # size of placebo time drift
                                       sd_trt = 0.80) {
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

  # ----- 1) shared TIME effect (moves placebo_EoT away from baseline) -----
  n_time <- round(prop_time * G)
  time_idx <- sample.int(G, n_time, replace = FALSE)

  time_effect <- rep(0, G)
  time_effect[time_idx] <- rnorm(n_time, mean = 0, sd = sd_time)

  # ----- 2) EXTRA treatment effect (moves treatment_EoT further) -----
  n_trt <- round(prop_trt * G)
  trt_idx <- sample.int(G, n_trt, replace = FALSE)

  trt_extra <- rep(0, G)
  trt_extra[trt_idx] <- rnorm(n_trt, mean = 0, sd = sd_trt)

  # placebo vs baseline
  group0_log2FC <- time_effect

  # treatment vs baseline = time + extra treatment
  group1_log2FC <- time_effect + trt_extra


  p@nTimepoints <- nTimepoints
  p@twoGroupDesign <- twoGroupDesign
  p@nSubjsPerGroup <- nSubjsPerGroup

  # RescueSim expects these names for the final timepoint effects by group
  p@deLog2FC <- list(
    time1_group0 = group0_log2FC,
    time1_group1 = group1_log2FC
  )

  p@propDE <- prop_trt

  p
}


#### rescueSim data ------------------------------------------------------------
simulate_rescue <- function(params) {
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
  cd <- SummarizedExperiment::colData(sim)
  cd$group <- relevel(factor(cd$group), ref = "Placebo")
  cd$time <- relevel(factor(cd$time), ref = "Baseline")
  SummarizedExperiment::colData(sim) <- cd
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

  # compute ground truth

  g0 <- SummarizedExperiment::rowData(sim)$deLog2FC.time1_group0
  g1 <- SummarizedExperiment::rowData(sim)$deLog2FC.time1_group1

  true_logFC <- g1 - g0

  ground_truth <- data.frame(
    gene = rownames(sim),
    true_logFC = true_logFC,
    trueDE = true_logFC != 0,
    placebo_effect = g0,
    treatment_effect = g1,
    simulation = "rescueSim"
  )

  simulation <- list(
    sce = sim,
    ground_truth = ground_truth
  )

  simulation
}


## Run simulations -------------------------------------------------------------
## SplatPop --------------------------------------------------------------------
set.seed(42)
params <- make_splatpop_params()
data <- simulate_splatpop(
  n_samples = 20,
  n_snps = 1000,
  params = params
)
splatpop_sim <- data$sce
splatpop_ground_truth <- data$ground_truth


set.seed(111)
params_11 <- make_splatpop_params()
data_11 <- simulate_splatpop(
  n_samples = 20,
  n_snps = 1000,
  params = params_11
)
splatpop_sim_11 <- data_11$sce
splatpop_ground_truth_11 <- data_11$ground_truth

set.seed(222)
params_22 <- make_splatpop_params()
data_22 <- simulate_splatpop(
  n_samples = 20,
  n_snps = 1000,
  params = params_22
)
splatpop_sim_22 <- data_22$sce
splatpop_ground_truth_22 <- data_22$ground_truth

set.seed(333)
params_33 <- make_splatpop_params()
data_33 <- simulate_splatpop(
  n_samples = 20,
  n_snps = 1000,
  params = params_33
)
splatpop_sim_33 <- data_33$sce
splatpop_ground_truth_33 <- data_33$ground_truth

set.seed(444)
params_44 <- make_splatpop_params()
data_44 <- simulate_splatpop(
  n_samples = 20,
  n_snps = 1000,
  params = params_11
)
splatpop_sim_44 <- data_44$sce
splatpop_ground_truth_44 <- data_44$ground_truth

## rescueSim -------------------------------------------------------------------
data("RecAM_sce", package = "rescueSim")
set.seed(42)
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(42)
data <- simulate_rescue(params)
rescuesim_sim <- data$sce
rescuesim_ground_truth <- data$ground_truth
# seed 123
data("RecAM_sce", package = "rescueSim")
set.seed(123)
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(123)
data <- simulate_rescue(params)
rescuesim_sim_2 <- data$sce
rescuesim_ground_truth_2 <- data$ground_truth
# seed 111
data("RecAM_sce", package = "rescueSim")
set.seed(111)
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(111)
data <- simulate_rescue(params)
rescuesim_sim_3 <- data$sce
rescuesim_ground_truth_3 <- data$ground_truth
# seed 222
data("RecAM_sce", package = "rescueSim")
set.seed(222)
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(222)
data <- simulate_rescue(params)
rescuesim_sim_4 <- data$sce
rescuesim_ground_truth_4 <- data$ground_truth
# seed 333
data("RecAM_sce", package = "rescueSim")
set.seed(333)
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(333)
data <- simulate_rescue(params)
rescuesim_sim_5 <- data$sce
rescuesim_ground_truth_5 <- data$ground_truth


# rescuesim norm lfc -----------------------------------------------------------

data("RecAM_sce", package = "rescueSim")
set.seed(42)
params <- make_rescuesim_norm_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(42)
data <- simulate_rescue(params)
rescuesim_sim_norm <- data$sce
rescuesim_ground_truth_norm <- data$ground_truth
# seed 123
data("RecAM_sce", package = "rescueSim")
set.seed(123)
params <- make_rescuesim_norm_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(123)
data <- simulate_rescue(params)
rescuesim_sim_2_norm <- data$sce
rescuesim_ground_truth_2_norm <- data$ground_truth
# seed 111
data("RecAM_sce", package = "rescueSim")
set.seed(111)
params <- make_rescuesim_norm_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(111)
data <- simulate_rescue(params)
rescuesim_sim_3_norm <- data$sce
rescuesim_ground_truth_3_norm <- data$ground_truth
# seed 222
data("RecAM_sce", package = "rescueSim")
set.seed(222)
params <- make_rescuesim_norm_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(222)
data <- simulate_rescue(params)
rescuesim_sim_4_norm <- data$sce
rescuesim_ground_truth_4_norm <- data$ground_truth
# seed 333
data("RecAM_sce", package = "rescueSim")
set.seed(333)
params <- make_rescuesim_norm_params(RecAM_sce, nSubjsPerGroup = 5)
set.seed(333)
data <- simulate_rescue(params)
rescuesim_sim_5_norm <- data$sce
rescuesim_ground_truth_5_norm <- data$ground_truth

# Pin simulations --------------------------------------------------------------
board <- pins::board_local("nrgn-data")
simulation <- c(
  "splatpop" = splatpop_sim,
  "splatpop_11" = splatpop_sim_11,
  "splatpop_22" = splatpop_sim_22,
  "splatpop_33" = splatpop_sim_33,
  "splatpop_44" = splatpop_sim_44,
  "rescueSim" = rescuesim_sim,
  "rescueSim_123" = rescuesim_sim_2,
  "rescueSim_111" = rescuesim_sim_3,
  "rescueSim_222" = rescuesim_sim_4,
  "rescueSim_333" = rescuesim_sim_5,
  "rescueSim norm" = rescuesim_sim_norm,
  "rescueSim_123 norm" = rescuesim_sim_2_norm,
  "rescueSim_111 norm" = rescuesim_sim_3_norm,
  "rescueSim_222 norm" = rescuesim_sim_4_norm,
  "rescueSim_333 norm" = rescuesim_sim_5_norm
)

ground_truth <- list(
  "splatpop" = splatpop_ground_truth,
  "splatpop_11" = splatpop_ground_truth_11,
  "splatpop_22" = splatpop_ground_truth_22,
  "splatpop_33" = splatpop_ground_truth_33,
  "splatpop_44" = splatpop_ground_truth_44,
  "rescueSim" = rescuesim_ground_truth,
  "rescueSim_123" = rescuesim_ground_truth_2,
  "rescueSim_111" = rescuesim_ground_truth_3,
  "rescueSim_222" = rescuesim_ground_truth_4,
  "rescueSim_333" = rescuesim_ground_truth_5,
  "rescueSim norm" = rescuesim_ground_truth_norm,
  "rescueSim_123 norm" = rescuesim_ground_truth_2_norm,
  "rescueSim_111 norm" = rescuesim_ground_truth_3_norm,
  "rescueSim_222 norm" = rescuesim_ground_truth_4_norm,
  "rescueSim_333 norm" = rescuesim_ground_truth_5_norm
)


pins::pin_write(simulation,
  name = "simulations_final",
  type = "rds",
  board = board,
  versioned = FALSE
)

pins::pin_write(ground_truth,
  name = "ground_truth_final",
  type = "rds",
  board = board,
  versioned = FALSE
)
