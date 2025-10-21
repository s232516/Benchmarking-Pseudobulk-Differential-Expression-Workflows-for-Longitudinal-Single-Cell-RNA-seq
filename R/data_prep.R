# if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
# BiocManager::install(c('splatter','scater','MAST','DESeq2','ROTS','scde'), ask = FALSE, update = FALSE)
# install.packages(c('data.table','tibble'))


library(splatter)
library(scater)
library(MAST)
library(DESeq2)
library(SummarizedExperiment)
library(data.table)


set.seed(42)


# Parameters (kept as in your original script)
n.samples <- 20
n.snps <- 1000
batchCells <- 50


default_params <- splatter::newSplatPopParams(
  eqtl.n = 0.5,
  batchCells = batchCells,
  similarity.scale = 6,
  condition.prob = c(0.5,0.5),
  eqtl.condition.specific = 0.5,
  cde.facLoc = 1,
  cde.facScale = 1
)


default_params <- splatter::setParam(default_params, 'dropout.type', 'experiment')
default_params <- splatter::setParam(default_params, 'dropout.shape', -1)
default_params <- splatter::setParam(default_params, 'dropout.mid', 5)


vcf <- splatter::mockVCF(n.snps = n.snps, n.samples = n.samples)
gff <- splatter::mockGFF()


sim.pop.cond <- splatter::splatPopSimulate(vcf = vcf, gff = gff, params = default_params, sparsify = FALSE)
sim.pop.cond <- scater::logNormCounts(sim.pop.cond)


# rename conditions
colData(sim.pop.cond)$Condition <- factor(colData(sim.pop.cond)$Condition, levels = c('Condition1','Condition2'), labels = c('Placebo','Treatment'))


# PCA/UMAP
sim.pop.cond <- scater::runPCA(sim.pop.cond, ncomponents = 10)
sim.pop.cond <- scater::runUMAP(sim.pop.cond, dimred = 'PCA')


# MAST (try-catch since it can be heavy)
mast_res <- NULL
try({
  sca <- MAST::SceToSingleCellAssay(sim.pop.cond, check_sanity = FALSE)
  colData(sca)$Condition <- factor(colData(sca)$Condition, levels = c('Placebo','Treatment'))
  zlmCond <- MAST::zlm(~ Condition, sca)
  summaryCond <- MAST::summary(zlmCond, doLRT = 'ConditionTreatment')
  summaryDt <- summaryCond$datatable
  fc <- summaryDt[contrast == 'ConditionTreatment' & component == 'logFC', .(gene = primerid, log2FoldChange = coef, stat = z)]
  pval <- summaryDt[contrast == 'ConditionTreatment' & component == 'H', .(gene = primerid, pvalue = `Pr(>Chisq)`)]
  mast_res <- merge(fc, pval, by = 'gene')
  mast_res$padj <- p.adjust(mast_res$pvalue, method = 'fdr')
  mast_res <- mast_res[order(mast_res$padj), ]
})


# DESeq2 pseudo-bulk (use raw counts)
rawData <- as.matrix(SummarizedExperiment::assay(sim.pop.cond, 'counts'))
coldata <- as.data.frame(SummarizedExperiment::colData(sim.pop.cond))
SampleName <- colnames(rawData)
condition <- coldata$Condition
meta_data <- data.frame(SampleName, condition)
rownames(meta_data) <- meta_data$SampleName


dds <- DESeq2::DESeqDataSetFromMatrix(countData = rawData, colData = meta_data, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
if (ncol(dds) > 1 && nrow(dds) > 1) {
  dds$condition <- relevel(dds$condition, ref = 'Placebo')
  dds <- DESeq2::DESeq(dds,   sfType = "poscounts")
  res <- DESeq2::results(dds)
  deseq_res <- as.data.frame(res)
  deseq_res$gene <- rownames(deseq_res)
  deseq_res <- deseq_res[order(deseq_res$padj), ]
} else {
  deseq_res <- NULL
}


# Save everything (fixed)
if (!dir.exists('data')) dir.create('data')
save(sim.pop.cond, mast_res, deseq_res, file = 'data/simData.rda')
cat('Saved data/simData.rda (fixed seed = 42)
')
