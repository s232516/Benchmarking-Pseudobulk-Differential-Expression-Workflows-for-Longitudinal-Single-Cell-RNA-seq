# Benchmarking Pseudobulk Differential Expression Workflows for Longitudinal Single-Cell RNA-seq

This repository contains an **R Shiny application**, built with **{golem}**, for exploring a simulation-based benchmark of differential expression workflows in **paired longitudinal single-cell RNA-seq (scRNA-seq) treatment studies**. The app provides interactive views of simulated data and benchmarking results, including dimensionality reduction, volcano plots, ROC curves, and confidence interval coverage.

The accompanying study compares **DESeq2, edgeR, limma-voom**, and two **MAST-on-pseudobulk** variants (using log2(CPM + 1) or log2(aggregated counts + 1)). Its main question is whether a workflow can reliably identify **treatment-specific changes over time**, rather than changes shared by both the Treatment and Placebo groups.

## Benchmark design

Each simulated subject belongs to either the **Treatment** or **Placebo** group and is sampled at **Baseline** and **End of Treatment (EoT)**. Counts are aggregated across cells within each subject–timepoint sample to create pseudobulk profiles. The target treatment effect is the difference between the longitudinal changes in the two groups (the treatment-by-time interaction).

- **Sample sizes:** 5, 10, or 20 subjects per group
- **Scenarios:** global null, time-only null, weak, medium, and strong treatment effects
- **Replicates:** 30 random seeds (202601–202630) per scenario and sample size
- **Total:** 5 scenarios × 3 sample sizes × 30 seeds = **450 simulated datasets**, with **1,940 genes** each
- **Simulation framework:** `rescueSim` (v0.99.0), calibrated from the `RecAM_sce` reference dataset

| Scenario | Genes with shared time effects | Genes with treatment effects | Mean absolute treatment effect (log2FC) |
| --- | ---: | ---: | ---: |
| Global null | 0% | 0% | 0.00 |
| Time-only null | 50% | 0% | 0.00 |
| Weak | 50% | 20% | 0.30 |
| Medium | 50% | 20% | 0.60 |
| Strong | 50% | 20% | 1.00 |

Importantly, **time-only genes are true negatives** for the treatment-by-time interaction, even though their expression changes between Baseline and EoT.

## Data-generation details

Baseline expression, dispersion, library sizes, and heterogeneity were estimated from `RecAM_sce`. The simulations preserve the paired design through a **shared subject-specific expression factor** across both visits, together with **independent sample-specific variation**.

| Characteristic | Fitted setting / summary |
| --- | --- |
| Cells per subject–timepoint sample | 67–136 |
| Log library size | Mean 7.4689; SD 0.4437 |
| Library-size heterogeneity parameter | 7.264 × 10⁻⁴ |
| Gene-level mean expression | Median 0.2833; IQR 0.0624–0.7121 |
| Gene dispersion | Median 0.07931; IQR 0.00010–0.37534 |
| Between-subject factor variability | Location −4.9584; SD 1.1944 |
| Between-sample factor variability | Location −5.3935; SD 0.9591 |

The final two rows describe **parameters of variability-factor distributions**, not directly observed correlations or biological variances. No single within-subject correlation coefficient was fixed.

**Illustrative dataset:** In the medium-effect scenario with 10 subjects per group and seed `202601`, there were 20 subjects, 40 paired pseudobulk samples, and 4,113 cells. The median was 106.5 cells per sample (range 72–130), 1,761 counts per cell (240–16,563), and 204,624 counts per pseudobulk sample (108,157–389,534). These are values from **one simulation**, not ranges pooled over all 450 datasets.

## Repository structure and running the app

- `app.R` — launches the Shiny application
- `R/` — app modules, plotting, and server logic
- `data/` — package data
- `dev/` — development scripts
- `tests/` — tests

After installing the R dependencies declared in `DESCRIPTION`, launch from the project root:

```r
shiny::runApp(".")
```

**Local data setup is required.** The current `R/app_server.R` reads precomputed results from a local **{pins} board** and contains a machine-specific `board_path`. Update that path to your own data location and supply the corresponding pins (including `simulations_final`, `ground_truth_final`, the comparison results, and confidence-interval results). Without these data objects, the app will not run as-is after cloning.

## Associated research

*Benchmarking Pseudobulk Differential Expression Workflows for Longitudinal Single-Cell RNA-seq Treatment Studies* — Nerea Gondra Zarain, Amaya Zaratiegui Pedrosa, and Jan Christian Refsgaard (IEEE BIBM 2026 workshop paper).
