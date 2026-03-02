#' Processed TCGA PAAD dataset (pdac)
#'
#' A reduced and cleaned subset of the TCGA pancreatic ductal adenocarcinoma (PAAD)
#' dataset, derived from The Cancer Genome Atlas (TCGA) PAAD cohort. This version,
#' `pdac`, is smaller and simplified for practical analyses and package examples.
#'
#' @details
#' This dataset was originally compiled and curated in the open-source \code{pdacR}
#' package by Torre-Healy et al. (2023), which harmonized and integrated the TCGA
#' PAAD gene expression and clinical data. The current version further reduces and
#' simplifies the data for efficient modeling demonstrations and survival analyses.
#'
#' The data frame includes:
#' \itemize{
#'   \item \strong{time}: Overall survival time in months.
#'   \item \strong{status}: Event indicator; 1 = event occurred, 0 = censored.
#'   \item \strong{treatment}: Binary treatment indicator; 1 = radiation therapy, 0 = control.
#'   \item \strong{age}: Age at initial pathologic diagnosis (numeric).
#'   \item \strong{sex}: Binary sex indicator; 1 = male, 0 = female.
#'   \item \strong{grade}: Tumor differentiation grade (ordinal; 1 = well, 2 = moderate, 3 = poor, 4 = undifferentiated).
#'   \item \strong{tumor.cellularity}: Tumor cellularity estimate (numeric).
#'   \item \strong{tumor.purity}: Tumor purity class (binary; 1 = high, 0 = low).
#'   \item \strong{absolute.purity}: Absolute purity estimate (numeric).
#'   \item \strong{moffitt.cluster}: Moffitt transcriptional subtype (binary; 1 = basal-like, 0 = classical).
#'   \item \strong{meth.leukocyte.percent}: DNA methylation leukocyte estimate (numeric).
#'   \item \strong{meth.purity.mode}: DNA methylation purity mode (numeric).
#'   \item \strong{stage}: Nodal stage indicator (binary; 1 = n1, 0 = n0).
#'   \item \strong{lymph.nodes}: Number of lymph nodes examined (numeric).
#'   \item \strong{Driver gene columns}: Expression values of key driver genes (e.g., KRAS, TP53, CDKN2A, SMAD4, BRCA1, BRCA2).
#'   \item \strong{Other gene columns}: Expression values of ~3,000 most variable non-driver genes (based on median absolute deviation).
#' }
#'
#' @format A data frame with rows corresponding to patients and columns as described above.
#'
#' @source \doi{10.1016/j.ccell.2017.07.007}
#'
#' @references
#' \itemize{
#'   \item Raphael BJ, et al. "Integrated genomic characterization of pancreatic ductal adenocarcinoma." Cancer Cell. 2017 Aug 14;32(2):185–203.e13. PMID: 28810144.
#'   \item Torre-Healy LA, Kawalerski RR, Oh K, et al. "Open-source curation of a pancreatic ductal adenocarcinoma gene expression analysis platform (pdacR) supports a two-subtype model." Communications Biology. 2023; https://doi.org/10.1038/s42003-023-04461-6.
#'   \item The Cancer Genome Atlas (TCGA), PAAD project, DbGaP: phs000178.
#' }
"pdac"
