# **Methylation heterogeneity reveals hidden epigenome–phenotype associations**

**Authors:** Shengwei Li<sup>&#42;</sup>, Chaorong Chen<sup>&#42;</sup>, Shuting Zhou, Shengying Wang, Ya Allen Cui, Wei Li<sup>✉</sup>, Jiejun Shi<sup>✉</sup>

<sup>*</sup> These authors contributed equally. <sup>✉</sup> Corresponding authors.

## Description

Epigenome-wide association studies (EWAS) have linked DNA methylation to many human diseases and traits, but most array-based studies rely on the average methylation level at each CpG. However, average methylation obscures the distribution of methylation states across DNA molecules, a disease-relevant feature captured by methylation heterogeneity but largely inaccessible to conventional array-based EWAS. Here we introduce heterogeneity-adjusted EWAS (haEWAS), which uses average-to-heterogeneity conversion models trained on deep methylomes from 205 reference samples to infer CpG-level methylation heterogeneity from standard methylation arrays and link this hidden methylation layer to human phenotypes at population scale. Applied to 572 cohorts comprising 139,869 samples, haEWAS identified 59,832 phenotype-associated loci across 46 phenotypes, including 42,768 missed by conventional EWAS. haEWAS-derived loci improved systemic lupus erythematosus case–control classification across independent cohorts, strengthened cis-regulatory interpretation in Alzheimer’s disease, improved survival prediction across 15 cancer types and resolved cell-type-specific methylation alterations in hepatocellular adenoma. Together, these findings establish methylation heterogeneity as a disease-relevant epigenetic signal inferable from methylation arrays, enabling epigenome–phenotype discovery beyond average methylation.

## Dependency

[![optparse](https://img.shields.io/badge/R-optparse-276DC3)](https://github.com/trevorld/r-optparse) [![data.table](https://img.shields.io/badge/R-data.table-276DC3)](https://github.com/Rdatatable/data.table) [![minfi](https://img.shields.io/badge/Bioconductor-minfi-3FA34D)](https://github.com/hansenlab/minfi) [![ChAMP](https://img.shields.io/badge/Bioconductor-ChAMP-3FA34D)](https://github.com/bioc/ChAMP) [![impute](https://img.shields.io/badge/Bioconductor-impute-3FA34D)](https://github.com/bioc/impute) [![sva](https://img.shields.io/badge/Bioconductor-sva-3FA34D)](https://github.com/bioc/sva) [![limma](https://img.shields.io/badge/Bioconductor-limma-3FA34D)](https://github.com/bioc/limma)

[![meta](https://img.shields.io/badge/R-meta-276DC3)](https://github.com/guido-s/meta) [![EpiDISH](https://img.shields.io/badge/R-EpiDISH-276DC3)](https://github.com/sjczheng/EpiDISH) [![EpiSCORE](https://img.shields.io/badge/R-EpiSCORE-276DC3)](https://github.com/aet21/EpiSCORE) [![CETS](https://img.shields.io/badge/R-CETS-276DC3)](https://www.r-forge.r-project.org/scm/?group_id=1511) [![survival](https://img.shields.io/badge/R-survival-276DC3)](https://github.com/therneau/survival) [![pROC](https://img.shields.io/badge/R-pROC-276DC3)](https://github.com/xrobin/pROC)

[![numpy](https://img.shields.io/badge/Python-numpy-3776AB)](https://github.com/numpy/numpy) [![pandas](https://img.shields.io/badge/Python-pandas-3776AB)](https://github.com/pandas-dev/pandas) [![SciPy](https://img.shields.io/badge/Python-SciPy-3776AB)](https://github.com/scipy/scipy) [![Matplotlib](https://img.shields.io/badge/Python-Matplotlib-3776AB)](https://github.com/matplotlib/matplotlib) [![seaborn](https://img.shields.io/badge/Python-seaborn-3776AB)](https://github.com/mwaskom/seaborn)

[![scikit--learn](https://img.shields.io/badge/Python-scikit--learn-3776AB)](https://github.com/scikit-learn/scikit-learn) [![lifelines](https://img.shields.io/badge/Python-lifelines-3776AB)](https://github.com/CamDavidsonPilon/lifelines) [![joblib](https://img.shields.io/badge/Python-joblib-3776AB)](https://github.com/joblib/joblib) [![matplotlib--venn](https://img.shields.io/badge/Python-matplotlib--venn-3776AB)](https://github.com/konstantint/matplotlib-venn)

## Repository structure

```text
Reproduction/
├── code/
│   ├── main_analysis/       # Main statistical analysis workflows
│   └── figures/             # Figure plotting scripts
├── data/                    # Panel-level CSV source data
├── figures/                 # Current PDF figure exports
├── figures_reproduction.ipynb
└── README.md
```

## Main data analysis

| Analysis | Main scripts |
| --- | --- |
| Methylation quality control, normalization, imputation, and batch correction | `code/main_analysis/01_preprocessing/methylation_preprocessing.R` |
| A2H transformation of methylation-average matrices | `code/main_analysis/02_A2H_transform/A2H_transform.py` |
| Single-cohort EWAS/haEWAS and cross-cohort meta-analysis | `code/main_analysis/03_association/single_cohort_association.R`; `cross_cohort_meta_analysis.R` |
| SLE methylation risk-score construction and validation | `code/main_analysis/04_SLE_MRS/sle_mrs_construction_validation.R` |
| AD Braak-stage association, meta-analysis, eQTM analysis, and blood MRS validation | `code/main_analysis/05_AD_Braak/` |
| TCGA single-CpG Cox association and repeated cross-validated survival prediction | `code/main_analysis/06_TCGA_survival/` |
| Cell-type-specific haEWAS with CellDMC | `code/main_analysis/07_CellDMC/cellDMC_analysis.R` |

## Figure guide

| Figure | Description | Figure code |
| --- | --- | --- |
| Figure 1 | Overview of the haEWAS framework. | Created with BioRender |
| Figure 2 | haEWAS expands the landscape of methylation–phenotype associations across 46 human phenotypes. | `code/figures/plot_fig2a_landscape_pvalue_pooled.py`; `code/figures/plot_fig2b_landscape_traits_pooled.py`; `code/figures/plot_fig2c_gain_score_pooled.py` |
| Figure 3 | haEWAS identifies SLE-associated loci and improves cross-cohort disease-status classification. | `code/figures/plot_fig3a_sle_venn_pooled.py`; `code/figures/plot_fig3b_sle_signed_manhattan_pooled.py`; `code/figures/plot_fig3c.py`; `code/figures/plot_fig3d.py`; `code/figures/plot_fig3e_g.py` |
| Figure 4 | haEWAS reveals neural functional features, stronger transcriptional regulatory associations, and blood-based predictive potential at AD Braak stage-associated loci. | `code/figures/plot_fig4a.py`–`plot_fig4f.py` |
| Figure 5 | haEWAS identifies pan-cancer overall-survival-associated loci and improves prognostic stratification. | `code/figures/Fig5a.py`; `code/figures/plot_fig5b.py`; `code/figures/Fig5cde.py` |
| Figure 6 | haEWAS uncovers additional cell type-specific differentially methylated loci in hepatocellular adenoma. | `code/figures/plot_fig6a.py`–`plot_fig6f.py` |
| Supplementary Figure 1 | Participant-level validation of A2H predictions in PGP-UK blood samples. | `code/figures/plot_figs1ab.py` |
| Supplementary Figure 2 | Comparison of A2H with *m*² and *m*³ in matched CpG groups. | `code/figures/plot_figs2ad.py` |
| Supplementary Figure 3 | CpG-specific prediction accuracy of A2H. | `code/figures/plot_figs3ad.py` |
| Supplementary Figure 4 | Robustness of A2H prediction performance to sequencing depth and read length. | `code/figures/plot_figs4ab.py`                               |
| Supplementary Figure 5 | Effect of cell-type exclusion on A2H prediction performance in the PGP-UK blood cohort. | `code/figures/plot_figs5ab.py` |
| Supplementary Figure 6 | Cell-type holdout validation of A2H in the EGA reference atlas. | `code/figures/plot_figs6.py` |
| Supplementary Figure 7 | Distributional properties and cross-cohort consistency of average methylation and methylation heterogeneity. | `code/figures/plot_figs7.py` |
| Supplementary Figure 8 | Distribution of haEWAS-specific, EWAS-specific, and common association loci across human phenotypes. | `code/figures/plot_figs8.py` |
| Supplementary Figure 9 | Contributions of hm(CAMDA) and hm(CHALM) to associations across 46 phenotypes. | `code/figures/plot_figs9.py` |
| Supplementary Figure 10 | Genomic inflation and functional annotation of SLE-associated loci identified by EWAS and haEWAS. | `code/figures/plot_figs10.py` |
| Supplementary Figure 11 | Sensitivity of genomic inflation to residual-PC adjustment in five SLE cohorts. | `code/figures/plot_figs11.py` |
| Supplementary Figure 12 | Sensitivity of SLE associations to residual-PC adjustment. | `code/figures/plot_figs12ab.py` |
| Supplementary Figure 13 | Leave-one-cohort-out validation of haEWAS associations in SLE. | `code/figures/plot_figs13.py` |
| Supplementary Figure 14 | Discovery-cohort feature optimization and external-validation workflow for SLE MRS construction and disease-status classification. | Created with BioRender |
| Supplementary Figure 15 | Cross-cohort consistency and functional enrichment of SLE-associated MRS features. | `code/figures/plot_figs15.py` |
| Supplementary Figure 16 | SLE classification using EWAS, hm(CAMDA), hm(CHALM), haEWAS, and combined MRS. | `code/figures/plot_figs16.py` |
| Supplementary Figure 17 | Comparison of EWAS/haEWAS loci and genomic inflation in the AD Braak stage meta-analysis. | `code/figures/plot_figs17.py` |
| Supplementary Figure 18 | Leave-one-cohort-out validation of haEWAS associations with AD Braak stage. | `code/figures/plot_figs18.py` |
| Supplementary Figure 19 | Classification of advanced Braak stage in the London blood cohort. | `code/figures/plot_figs19.py` |
| Supplementary Figure 20 | Contributions of hm(CAMDA) and hm(CHALM) to survival associations across 15 TCGA cancer types. | `code/figures/plot_figs20.py` |
| Supplementary Figure 21 | Pan-cancer prognostic stratification using OS-associated loci identified by haEWAS and EWAS. | `code/figures/plot_figs21.py` |
| Supplementary Figure 22 | Contributions of hm(CAMDA) and hm(CHALM) to DMCTs in hepatocellular adenoma. | `code/figures/plot_figs22.py` |

## Reproducing figures

Run `figures_reproduction.ipynb` in Jupyter to reproduce the figures. It runs the plotting scripts in `code/figures/`, reads the panel tables from `data/`, displays the plots, and saves PDFs in `figures/`.
