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


#### rescueSim data ------------------------------------------------------------
simulate_rescue <- function(params) {
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
  cd <- SummarizedExperiment::colData(sim)
  cd$group <- relevel(factor(cd$group), ref = "Placebo")
  cd$time  <- relevel(factor(cd$time),  ref = "Baseline")
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

## Jan -------------------------------------------------------------------------
simulate_jan <- function(nSubjsPerGroup = 20,
                         nCellsPerSubjectTimepoint = 100,
                         nGenes = 100) {
  set.seed(42)

  # Jan-style simulation setup
  n_subjects <- nSubjsPerGroup * 2
  n_timepoints <- 2

  # Experimental design
  design_metadata <- df <- expand.grid(
    Subject = 1:n_subjects,
    Time = 0:(n_timepoints - 1)
  ) |>
    dplyr::mutate(Sample = dplyr::row_number()) |>
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

  # Add effects manually
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
  sce$subjectID <- factor(sce$Subject)
  sce$sampleID <- factor(sce$Sample)
  sce$group <- factor(ifelse(sce$Treatment == 0, "Placebo", "Treatment"))
  sce$time <- factor(ifelse(sce$Time == 0, "Baseline", "EoT"))
  sce$Condition <- paste0(
    sce$group, "_", sce$time
  )
  cd <- SummarizedExperiment::colData(sce)
  cd$group <- relevel(factor(cd$group), ref = "Placebo")
  cd$time  <- relevel(factor(cd$time),  ref = "Baseline")
  SummarizedExperiment::colData(sce) <- cd



  sim_jan <- sce

  sim_jan <- scater::runPCA(sim_jan, ncomponents = 10)
  sim_jan <- scater::runUMAP(sim_jan, dimred = "PCA")
  # ground truth

  placebo_shift <- rep(0, nGenes)
  treatment_shift <- rep(0, nGenes)

  placebo_shift[placebo_genes] <- 1.0      # +1 at EoT for placebo genes
  treatment_shift[treatment_genes] <- 1.5  # +1.5 at EoT for treatment genes

  true_logFC <- (placebo_shift + treatment_shift) - placebo_shift

  ground_truth <- data.frame(
    gene = paste0("Gene_", 1:nGenes),
    true_logFC = true_logFC,
    trueDE = true_logFC != 0,
    placebo_effect = placebo_shift,
    treatment_effect = placebo_shift + treatment_shift,
    simulation = "janSim"
  )


  simulation <- list(
    sce = sim_jan,
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
splatpo_sim <- data$sce
splatpop_ground_truth <- data$ground_truth


## Jan -------------------------------------------------------------------------
data <- simulate_jan(
  nSubjsPerGroup = 20,
  nCellsPerSubjectTimepoint = 100,
  nGenes = 100
)
jan_sim <- data$sce
jan_ground_truth <- data$ground_truth


## rescueSim -------------------------------------------------------------------
data("RecAM_sce", package = "rescueSim")
params <- make_rescuesim_params(RecAM_sce, nSubjsPerGroup = 5)
data <- simulate_rescue(params)
rescuesim_sim <- data$sce
rescuesim_ground_truth <- data$ground_truth


# Pin simulations --------------------------------------------------------------
board <- pins::board_local("nrgn-data")
simulation <- c(
  "splatpop" = splatpop_sim,
  "rescueSim" = rescuesim_sim,
  "jan" = jan_sim
)

ground_truth <- list(
  "splatpop"   = splatpop_ground_truth,
  "rescueSim"  = rescuesim_ground_truth,
  "jan"        = jan_ground_truth
)


pins::pin_write(simulation,
  name = "simulations",
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
