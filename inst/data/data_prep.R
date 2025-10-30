# Create/read a board
board <- pins::board_local("nrgn-data")

# This script simulates single-cell RNA-seq data with two conditions

# TRT -----------------------------------------------------------------------
set.seed(42)


# Parameters (kept as in your original script)
n_samples <- 20
n_snps <- 1000
batchcells <- 100


default_params <- splatter::newSplatPopParams(
  eqtl.n = 0.5,
  batchCells = batchcells,
  similarity.scale = 6,
  condition.prob = c(0.5, 0.5),
  eqtl.condition.specific = 0.5,
  cde.facLoc = 1,
  cde.facScale = 1
)
default_params <- splatter::setParam(
  default_params, "dropout.type",
  "experiment"
)
default_params <- splatter::setParam(default_params, "dropout.shape", -1)
default_params <- splatter::setParam(default_params, "dropout.mid", 5)


vcf <- splatter::mockVCF(n.snps = n_snps, n.samples = n_samples)
gff <- splatter::mockGFF()


sim_pop_cond <- splatter::splatPopSimulate(
  vcf = vcf, gff = gff,
  params = default_params,
  sparsify = TRUE # to make it more realisic, for dropout
)
sim_pop_cond <- scater::logNormCounts(sim_pop_cond)
sim_pop_cond <- scater::runPCA(sim_pop_cond, ncomponents = 10)
sim_pop_cond <- scater::runUMAP(sim_pop_cond, dimred = "PCA")

# rename conditions
SummarizedExperiment::colData(sim_pop_cond)$Condition <- factor(
  SummarizedExperiment::colData(sim_pop_cond)$Condition,
  levels = c("Condition1", "Condition2"),
  labels = c("Placebo", "Treatment")
)

## MAST (try-catch since it can be heavy) -------------------------------------
mast_res <- NULL
try({
  sca <- MAST::SceToSingleCellAssay(sim_pop_cond, check_sanity = FALSE)
  SummarizedExperiment::colData(sca)$Condition <- factor(
    SummarizedExperiment::colData(sca)$Condition,
    levels = c("Placebo", "Treatment")
  )
  zlmcond <- MAST::zlm(~Condition, sca)
  summarycond <- MAST::summary(zlmcond, doLRT = "ConditionTreatment")
  summarydt <- summarycond$datatable

  # Using dplyr to avoid .() / data.table issues
  fc <- summarydt |>
    dplyr::filter(contrast == "ConditionTreatment", component == "logFC") |>
    dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

  pval <- summarydt |>
    dplyr::filter(contrast == "ConditionTreatment", component == "H") |>
    dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

  mast_res <- merge(fc, pval, by = "gene")
  mast_res$padj <- p.adjust(mast_res$pvalue, method = "fdr")
  mast_res <- mast_res[order(mast_res$padj), ] |>
    dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

  mast_res$analysis <- "trt"
  mast_res$method <- "MAST"
})
## MAST (random) ---------------------------------------------------------------
mast_rand <- NULL
try({
  sca <- MAST::SceToSingleCellAssay(sim_pop_cond, check_sanity = FALSE)
  SummarizedExperiment::colData(sca)$Condition <- factor(
    SummarizedExperiment::colData(sca)$Condition,
    levels = c("Placebo", "Treatment")
  )

  SummarizedExperiment::colData(sca)$Sample <- factor(
    SummarizedExperiment::colData(sca)$Sample
  )
  table(colData(sca)$Sample)
  zlmcond <- MAST::zlm(~ Condition + (1 | Sample), sca,
    method = "glmer",
    ebayes = FALSE,
    strictConvergence = TRUE
  )
  summarycond <- MAST::summary(zlmcond, doLRT = "ConditionTreatment")
  summarydt <- summarycond$datatable

  # Using dplyr to avoid .() / data.table issues
  fc <- summarydt |>
    dplyr::filter(contrast == "ConditionTreatment", component == "logFC") |>
    dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

  pval <- summarydt |>
    dplyr::filter(contrast == "ConditionTreatment", component == "H") |>
    dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

  mast_rand <- merge(fc, pval, by = "gene")
  mast_rand$padj <- p.adjust(mast_rand$pvalue, method = "fdr")
  mast_rand <- mast_rand[order(mast_rand$padj), ] |>
    dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

  mast_rand$analysis <- "trt"
  mast_rand$method <- "MAST with random effects"
})


## DESeq2 pseudo-bulk (use raw counts) ----------------------------------------

pseudobulk <- scater::aggregateAcrossCells(
  sim_pop_cond,
  ids = SummarizedExperiment::colData(sim_pop_cond)$Sample,
  statistics = "sum"
)


# 2. Create metadata (one row per sample)
metadata <- as.data.frame(metadata <- as.data.frame(colData(pseudobulk)) |>
  dplyr::select(Sample, Condition))
metadata$Condition <- factor(metadata$Condition)
metadata$Sample <- factor(metadata$Sample)


# 3. DESeq2 analysis (account for batch effects)
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = SummarizedExperiment::assay(
    pseudobulk,
    "counts"
  ),
  colData = metadata,
  design = ~Condition
)

keep <- rowSums(SingleCellExperiment::counts(dds)) >= 10
dds <- dds[keep, ]
if (ncol(dds) > 1 && nrow(dds) > 1) {
  dds$Condition <- stats::relevel(dds$Condition, ref = "Placebo")
  dds <- DESeq2::DESeq(dds, sfType = "poscounts")
  res <- DESeq2::results(dds)
  deseq_res <- as.data.frame(res)
  deseq_res$gene <- rownames(deseq_res)
  deseq_res <- deseq_res[order(deseq_res$padj), ]
  deseq_res$analysis <- "trt"
  deseq_res$method <- "DESeq"
} else {
  deseq_res <- NULL
}

# TRT + Time ------------------------------------------------------------------
## Simulate with SplatPop -----------------------------------------------------
# sce <- sim_pop_cond
#
# # Duplicate placebo baseline to make treatment_baseline & placebo_eot
# placebo_idx <- which(sce$Condition == "Placebo")
# treatment_baseline <- sce[, placebo_idx]  # copy
# placebo_eot <- sce[, placebo_idx]          # copy
#
# # Combine all into a new SCE
# sce_extended <- SummarizedExperiment::cbind(sce[, placebo_idx],
#                                             treatment_baseline,
#                                             placebo_eot,
#                       sce[, sce$Condition == "Treatment"])  # Treatment EOT
#
# sce_extended$Condition <- factor(
#   c(rep("Placebo_baseline", length(placebo_idx)),
#     rep("Treatment_baseline", length(placebo_idx)),
#     rep("Placebo_EOT", length(placebo_idx)),
#     rep("Treatment_EOT", length(which(sce$Condition == "Treatment"))))
# )
# sce_extended <- scater::logNormCounts(sce_extended)
# sce_extended <- scater::runPCA(sce_extended, ncomponents = 10)
# scater::plotPCA(sce_extended, colour_by = "Condition")
# sce_extended <- scater::runUMAP(sce_extended, dimred = "PCA")
#

## Simulate rescue data with rescueSim ----------------------------------------

set.seed(42)
data("RecAM_sce")
## Subset data to 40 genes
RecAM_sce_small <- RecAM_sce[1:40, ]

# Create empty parameter object
myParams <- rescueSim::RescueSimParams()

## Estimate parameters using full RecAM_sce object
myParams <- rescueSim::estRescueSimParams(RecAM_sce,
  sampleVariable = "sampleID",
  subjectVariable = "subjectID",
  timepointVariable = "time",
  groupVariable = NULL, nonDEGs = NULL,
  cellParamsByCondition = FALSE
)
# Fill in basic experimental design
myParams@nTimepoints <- 2 # T0 and T1
myParams@twoGroupDesign <- TRUE # two groups (treatment vs control)
myParams@nSubjsPerGroup <- 5 # 10 subjects per group

# baseline (time1) is the reference
myParams@propDE <- 0.20 # 20% of genes are DE at any change

# assign to parameters
myParams@deLog2FC <- 1

sim <- rescueSim::simRescueData(myParams)
sim <- scater::logNormCounts(sim)
sim <- scater::runPCA(sim, ncomponents = 10)
sim <- scater::runUMAP(sim, dimred = "PCA")

SummarizedExperiment::colData(sim)$Sample <- SummarizedExperiment::colData(sim)$SampleID
# Rename group
SummarizedExperiment::colData(sim)$group <- factor(
  SummarizedExperiment::colData(sim)$group,
  levels = c("group0", "group1"),
  labels = c("Placebo", "Treatment")
)

# Rename time
SummarizedExperiment::colData(sim)$time <- factor(
  SummarizedExperiment::colData(sim)$time,
  levels = c("time0", "time1"),
  labels = c("Baseline", "EoT")
)

# Create a combined condition column
SummarizedExperiment::colData(sim)$Condition <- paste0(
  SummarizedExperiment::colData(sim)$group, "_",
  SummarizedExperiment::colData(sim)$time
)

### MAST (try-catch since it can be heavy) -------------------------------------
mast_recam <- NULL
try({
  sca <- MAST::SceToSingleCellAssay(sim, check_sanity = FALSE)

  zlm <- MAST::zlm(~ group * time, sca)
  summary <- MAST::summary(zlm, doLRT = "groupTreatment:timeEoT")
  summarydt <- summary$datatable

  fc <- summarydt |>
    dplyr::filter(contrast == "groupTreatment:timeEoT", component == "logFC") |>
    dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

  pval <- summarydt |>
    dplyr::filter(contrast == "groupTreatment:timeEoT", component == "H") |>
    dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

  mast_recam <- merge(fc, pval, by = "gene")
  mast_recam$padj <- p.adjust(mast_recam$pvalue, method = "fdr")
  mast_recam <- mast_recam[order(mast_recam$padj), ] |>
    dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

  mast_recam$analysis <- "trt*time"
  mast_recam$method <- "MAST (recam)"
})

mast_recam




### DESeq2 pseudo-bulk (use raw counts) ----------------------------------------


# 1. Aggregate counts to pseudobulk level (one row per sample)
pseudobulk <- pseudobulk <- scater::aggregateAcrossCells(
  sim,
  ids = SummarizedExperiment::colData(
    sim
  )[, c("subjectID", "sampleID")],
  statistics = "sum"
)


# 2. Create clean metadata (one row per sample)
metadata <- data.frame(
  sampleID = factor(SummarizedExperiment::colData(pseudobulk)$sampleID),
  subjectID = factor(SummarizedExperiment::colData(pseudobulk)$subjectID),
  Treatment = factor(SummarizedExperiment::colData(pseudobulk)$group),
  Time = factor(SummarizedExperiment::colData(pseudobulk)$time)
)


# 3. Build DESeq2 dataset with interaction
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = SummarizedExperiment::assay(pseudobulk, "counts"),
  colData = metadata,
  design = ~ Treatment + Time + Treatment:Time
)

# 4. Filter low-count genes
keep <- rowSums(DESeq2::counts(dds)) >= 10
dds <- dds[keep, ]

# 5. Relevel reference groups (so baseline is Placebo at Baseline)
dds$Treatment <- relevel(dds$Treatment, ref = "Placebo")
dds$Time <- relevel(dds$Time, ref = "Baseline")

keep <- rowSums(SingleCellExperiment::counts(dds)) >= 10
dds <- dds[keep, ]
if (ncol(dds) > 1 && nrow(dds) > 1) {
  dds <- DESeq2::DESeq(dds, sfType = "poscounts")

  recam <- DESeq2::results(dds)
  deseq_recam <- as.data.frame(recam)
  deseq_recam$gene <- rownames(deseq_recam)
  deseq_recam <- deseq_recam[order(deseq_recam$padj), ]
  deseq_recam$analysis <- "trt*time"
  deseq_recam$method <- "DESeq (recam)"
} else {
  deseq_recam <- NULL
}

## MAST (trt * time) NOT WORKING ---------------------------------------
mast_recam_rand <- NULL
try({
  # Convert to SingleCellAssay
  sca <- MAST::SceToSingleCellAssay(sim, check_sanity = FALSE)

  # Ensure factors are set correctly
  SummarizedExperiment::colData(sca)$group <- factor(
    SummarizedExperiment::colData(sca)$group,
    levels = c("Placebo", "Treatment")
  )
  SummarizedExperiment::colData(sca)$time <- factor(
    SummarizedExperiment::colData(sca)$time,
    levels = c("Baseline", "EoT")
  )
  SummarizedExperiment::colData(sca)$subjectID <- factor(
    SummarizedExperiment::colData(sim)$subjectID
  )
  SummarizedExperiment::colData(sca)$sampleID <- factor(
    SummarizedExperiment::colData(sim)$sampleID
  )

  # Fit MAST model with random intercepts for subjectID and sampleID
  zlm_recam_rand <- MAST::zlm(
    ~ group * time + (1 | subjectID) + (1 | sampleID),
    sca,
    method = "glmer",
    ebayes = FALSE,
    strictConvergence = TRUE,
    silent = TRUE
  )

  # Perform likelihood ratio test for the interaction term
  summary_rand <- MAST::summary(zlm_recam_rand, doLRT = "groupTreatment:timeEoT")
  summarydt_rand <- summary_rand$datatable

  # Extract logFC and stats
  fc_rand <- summarydt_rand |>
    dplyr::filter(contrast == "groupTreatment:timeEoT", component == "logFC") |>
    dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

  # Extract p-values
  pval_rand <- summarydt_rand |>
    dplyr::filter(contrast == "groupTreatment:timeEoT", component == "H") |>
    dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

  # Merge results
  mast_recam_rand <- merge(fc_rand, pval_rand, by = "gene")
  mast_recam_rand$padj <- p.adjust(mast_recam_rand$pvalue, method = "fdr")
  mast_recam_rand <- mast_recam_rand[order(mast_recam_rand$padj), ] |>
    dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

  mast_recam_rand$analysis <- "trt*time"
  mast_recam_rand$method <- "MAST (trt*time with random effects)"
})

mast_recam_rand


# Jan Simulation --------------------------------------------------------------
# Set parameters
set.seed(42)
n_subjects_per_arm <- 20
n_cells_per_subject_timepoint <- 50
n_genes <- 100
n_cell_types <- 3


# UPDATED: genes 1-20 have placebo effect, genes 15-25 have treatment effect
placebo_genes <- 1:20
treatment_genes <- 15:25 # Now 15-20 have BOTH effects
both_effects_genes <- intersect(placebo_genes, treatment_genes) # genes 15-20


# Create experimental design
# First 20 subjects Placebo
design_metadata <- expand.grid(
  Subject = 1:40,
  Time = c(0, 1)
) |>
  dplyr::mutate(Treatment = ifelse(Subject <= 20, 0, 1)) |>
  dplyr::slice(rep(1:dplyr::n(), each = n_cells_per_subject_timepoint)) |>
  # 50 replicates of each (50 cells)
  dplyr::mutate(CellID = paste0("Cell_", 1:dplyr::n()))


# print(paste("Total cells:", nrow(design_metadata)))


# Use splatter just for realistic count structure, then modify
params <- splatter::newSplatParams()
params <- splatter::setParam(params, "nGenes", n_genes)
params <- splatter::setParam(params, "batchCells", nrow(design_metadata))
params <- splatter::setParam(params, "group.prob", 1) # Single group
params <- splatter::setParam(params, "de.prob", 0) # No DE from splatter


# Generate base expression
sce <- splatter::splatSimulate(params, verbose = FALSE)


# Now manually add the effects we want
counts_matrix <- SingleCellExperiment::counts(sce)
log_counts <- log2(counts_matrix + 1)

# Create mask ONCE before loop

modification_mask <- counts_matrix > 0 # TRUE for expressed cells
zero_mask <- counts_matrix == 0
induction_probs <- matrix(runif(prod(dim(zero_mask))) < 0.2,
  nrow = nrow(zero_mask)
)

modification_mask <- modification_mask | (zero_mask & induction_probs)

for (gene_idx in 1:n_genes) {
  gene_name <- paste0("Gene_", gene_idx)
  for (cell_idx in 1:nrow(design_metadata)) {
    subject <- design_metadata$Subject[cell_idx]
    time <- design_metadata$Time[cell_idx]
    treatment <- design_metadata$Treatment[cell_idx]

    if (!modification_mask[gene_idx, cell_idx]) {
      counts_matrix[gene_idx, cell_idx] <- 0 # Keep as zero
      next # Skip to next iteration
    }
    # Only reaches here if mask allows modification
    base_expr <- log_counts[gene_idx, cell_idx]

    if (gene_idx %in% placebo_genes) {
      base_expr <- base_expr + time * 1.0
    }
    if (gene_idx %in% treatment_genes) {
      base_expr <- base_expr + time * treatment * 1.5
    }

    new_count <- round(2^base_expr - 1)
    counts_matrix[gene_idx, cell_idx] <- ifelse(new_count < 2, 0, new_count)
  }
}


# Update the SCE object
SingleCellExperiment::counts(sce) <- as(counts_matrix, "dgCMatrix")
SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 1)


# Add metadata
colData(sce) <- cbind(SummarizedExperiment::colData(sce), design_metadata)
colnames(sce) <- design_metadata$CellID
rownames(sce) <- paste0("Gene_", 1:n_genes)

s <- as.data.frame(as.matrix(SingleCellExperiment::logcounts(sce)))

# Ground truth
ground_truth <- data.frame(
  Gene = paste0("Gene_", 1:n_genes),
  PlaceboEffect = 1:n_genes %in% placebo_genes,
  TreatmentEffect = 1:n_genes %in% treatment_genes,
  BothEffects = 1:n_genes %in% both_effects_genes
)

# Make a grouping column for visualization
sce$Group <- ifelse(sce$Treatment == 0, "Placebo", "Drug")
sce$Condition <- paste0(sce$Group, "_", ifelse(sce$Time == 0, "Baseline", "End"))

# PCA ----
sce <- scater::runPCA(sce, ncomponents = 10)

# UMAP (using PCA as input) ----
sce <- scater::runUMAP(sce, dimred = "PCA")


# Create final object
simulation_results <- list(
  sce = sce,
  ground_truth = ground_truth,
  parameters = list(
    placebo_genes = placebo_genes,
    treatment_genes = treatment_genes,
    both_effects_genes = both_effects_genes
  )
)

simulation_results
## MAST --------------------------------------------------------------
# Prepare data for MAST
# Convert to SingleCellAssay object

sca <- MAST::SceToSingleCellAssay(sce)

# Filter lowly expressed genes
sca <- sca[rowSums(assay(sca)) > 0, ] # Remove all-zero genes

# Create metadata dataframe with proper factors

SummarizedExperiment::colData(sca)$Subject <- factor(
  SummarizedExperiment::colData(sca)$Subject
)
SummarizedExperiment::colData(sca)$Time <- factor(
  SummarizedExperiment::colData(sca)$Time
)
SummarizedExperiment::colData(sca)$Treatment <- factor(
  SummarizedExperiment::colData(sca)$Treatment,
  levels = c(0, 1),
  labels = c("Placebo", "Drug")
)

# Model specification (accounts for within-subject pairing)
formula <- ~ Time + Treatment + Time:Treatment + (1 | Subject)

# Fit hurdle model
zlm.fit <- MAST::zlm(formula, sca,
  method = "glmer",
  ebayes = FALSE,
  strictConvergence = TRUE,
  silent = TRUE
)

# Test for treatment-time interaction

contrast <- c("Time1:TreatmentDrug") # The interaction term
summary.fit <- summary(zlm.fit, doLRT = contrast)

# Extract results
results <- summary.fit$datatable

# Using dplyr to avoid .() / data.table issues
fc <- results |>
  dplyr::filter(contrast == "Time1:TreatmentDrug", component == "logFC") |>
  dplyr::transmute(gene = primerid, log2FoldChange = coef, stat = z)

pval <- results |>
  dplyr::filter(contrast == "Time1:TreatmentDrug", component == "H") |>
  dplyr::transmute(gene = primerid, pvalue = `Pr(>Chisq)`)

mast_jan <- merge(fc, pval, by = "gene")
mast_jan$padj <- p.adjust(mast_jan$pvalue, method = "fdr")
mast_jan <- mast_jan[order(mast_jan$padj), ] |>
  dplyr::mutate(lfcSE = abs(log2FoldChange / stat))

mast_jan$analysis <- "trt*time (splatter)"
mast_jan$method <- "MAST"


## DESeq -----------------------------------------------------------
# Create pseudobulk counts
pseudobulk <- scater::aggregateAcrossCells(sce,
  ids = SummarizedExperiment::colData(
    sce
  )[, c("Subject", "Time", "Treatment")],
  statistics = "sum"
)


# Create metadata for pseudobulk
metadata <- as.data.frame(SummarizedExperiment::colData(pseudobulk))
metadata$Time <- factor(metadata$Time)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Subject <- factor(metadata$Subject)


# DESeq2 dataset
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = SummarizedExperiment::assay(
    pseudobulk,
    "counts"
  ),
  colData = metadata,
  design = ~ Time + Treatment + Time:Treatment
)


# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run analysis
dds <- DESeq2::DESeq(dds)


# Extract results
res <- DESeq2::results(dds, name = "Time1.Treatment1", alpha = 0.05)
deseq_jan <- as.data.frame(res)
deseq_jan$gene <- rownames(deseq_jan)
deseq_jan <- deseq_jan[order(deseq_jan$padj), ]
deseq_jan$analysis <- "trt*time(Splatter)"
deseq_jan$method <- "DESeq"

#  Pin data ----------------------------------------------------------
method_data <- rbind(mast_res, mast_rand, deseq_res, mast_recam,
  deseq_recam, deseq_jan,
  fill = TRUE
)
simultation <- c(
  "Treatment" = sim_pop_cond,
  "Treatment + Time (rescueSim)" = sim,
  "Treatment + Time (splatter)" = sce
)
board |> pins::pin_write(simultation,
  name = "simulations", type = "rds",
  versioned = FALSE
)
board |> pins::pin_write(method_data,
  name = "method_data", type = "rds",
  versioned = FALSE
)
