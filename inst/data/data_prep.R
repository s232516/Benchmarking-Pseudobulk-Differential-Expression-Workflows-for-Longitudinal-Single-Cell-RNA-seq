# Create/read a board
board <- pins::board_local("nrgn-data")

# This script simulates single-cell RNA-seq data with two conditions

# TRT -----------------------------------------------------------------------

set.seed(42)


# Parameters (kept as in your original script)
n_samples <- 20
n_snps <- 1000
batchcells <- 50


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
  sparsify = FALSE
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

# MAST (try-catch since it can be heavy)
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



# DESeq2 pseudo-bulk (use raw counts)
rawdata <- as.matrix(SummarizedExperiment::assay(sim_pop_cond, "counts"))
coldata <- as.data.frame(SummarizedExperiment::colData(sim_pop_cond))
samplename <- colnames(rawdata)
condition <- coldata$Condition
meta_data <- data.frame(samplename, condition)
rownames(meta_data) <- meta_data$samplename


dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = rawdata, colData = meta_data,
  design = ~condition
)
keep <- rowSums(SingleCellExperiment::counts(dds)) >= 10
dds <- dds[keep, ]
if (ncol(dds) > 1 && nrow(dds) > 1) {
  dds$condition <- stats::relevel(dds$condition, ref = "Placebo")
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


method_data <- rbind(mast_res, deseq_res, fill=TRUE)


# Pin your data objects
board |> pins::pin_write(sim_pop_cond, name = "sim_pop_cond", type = "rds",
                         versioned = FALSE)
board |> pins::pin_write(method_data, name = "method_data", type = "rds",
                         versioned = FALSE)
# Time ------------------------------------------------------------------------
# TRT + Time ------------------------------------------------------------------
