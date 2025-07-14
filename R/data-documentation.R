#' TCGA PAAD pancreatic adenocarcinoma dataset (Raphael et al., 2017)
#'
#' This dataset contains RNA-seq gene expression data and associated sample annotations
#' from pancreatic ductal adenocarcinoma (PAAD) tumors in The Cancer Genome Atlas (TCGA),
#' as compiled and distributed in the \code{pdacR} package.
#'
#' The expression data (`ex`) includes 20,531 genes across 181 primary tumors.
#' The sample information (`sampInfo`) includes 55 clinical and molecular variables
#' (e.g., tumor purity, mutations, subtypes, survival outcomes).
#' Additional metadata and gene feature annotations are provided.
#'
#' @format A list with 4 elements:
#' \describe{
#'   \item{ex}{Numeric matrix with 20,531 genes (rows) and 181 tumor samples (columns).}
#'   \item{sampInfo}{Data frame with 181 rows and 55 columns, containing sample-level annotations.}
#'   \item{metadata}{List with experiment details, including log-transformation flag,
#'    reference citation, accession IDs, description, survival column names, and default plot settings.}
#'   \item{featInfo}{Data frame with 20,531 rows and 2 columns, containing gene symbols and Entrez IDs.}
#' }
#'
#' @details
#' This dataset is derived from the TCGA PAAD cohort (accession: DbGaP phs000178, BroadGDAC: PAAD).
#' It was processed and summarized in Raphael et al. (2017), and further curated, integrated, and
#' distributed via the open-source \code{pdacR} package as described by Torre-Healy et al. (2023).
#'
#' @source 
#' Raphael BJ, et al. "Integrated genomic characterization of pancreatic ductal adenocarcinoma."
#' Cancer Cell. 2017 Aug 14;32(2):185–203.e13. PMID: 28810144.
#' 
#' Torre-Healy LA, Kawalerski RR, Oh K, et al. "Open-source curation of a pancreatic ductal
#' adenocarcinoma gene expression analysis platform (pdacR) supports a two-subtype model."
#' Communications Biology. 2023; https://doi.org/10.1038/s42003-023-04461-6.
#'
#' @references
#' \itemize{
#'   \item Raphael BJ, et al., Cancer Cell, 2017, PMID: 28810144.
#'   \item The Cancer Genome Atlas (TCGA), PAAD project, DbGaP: phs000178.
#'   \item Torre-Healy LA, et al., Communications Biology, 2023, https://doi.org/10.1038/s42003-023-04461-6.
#'   \item pdacR package, MIT license.
#' }
"pdac"

#' Processed PDAC TCGA dataset (reduced version)
#'
#' A processed subset of the TCGA pancreatic ductal adenocarcinoma (PAAD) dataset, prepared for flexible survival and molecular analyses. This reduced dataset contains clinical covariates, survival data, expression values of selected driver genes, and the top 3,000 most variable non-driver genes.
#'
#'
#' @details
#' This data frame includes:
#' \itemize{
#'   \item \strong{time}: Overall survival time in months.
#'   \item \strong{status}: Censoring indicator: 1 = event, 0 = censored.
#'   \item \strong{treatment}: Treatment indicator: 1 = radiation therapy, 0 = control
#'   \item \strong{age}: Age at initial pathologic diagnosis (numeric).
#'   \item \strong{sex}: Binary sex indicator: 1 = male, 0 = female.
#'   \item \strong{grade}: Tumor differentiation grade (ordinal: 1 = well, 2 = moderately, 3 = poorly, 4 = undifferentiated).
#'   \item \strong{tumor.cellularity}: Pathologist-reviewed tumor cellularity (numeric).
#'   \item \strong{tumor.purity}: Tumor purity class (binary: 1 = high, 0 = low).
#'   \item \strong{absolute.purity}: Absolute purity estimate (numeric).
#'   \item \strong{moffitt.cluster}: Moffitt transcriptional subtype (binary: 1 = basal-like, 0 = classical).
#'   \item \strong{meth.leukocyte.percent}: DNA methylation leukocyte percent estimate (numeric).
#'   \item \strong{meth.purity.mode}: DNA hypermethylation mode purity (numeric).
#'   \item \strong{stage}: Nodal stage indicator (binary: 1 = n1, 0 = n0).
#'   \item \strong{lymph.nodes}: Number of lymph nodes examined (numeric).
#'   \item \strong{Driver gene expression columns}: Expression values of known driver genes (KRAS, TP53, CDKN2A, SMAD4, BRCA1, BRCA2, etc.).
#'   \item \strong{Other gene expression columns}: Expression values for the top ~3,000 most variable non-driver genes (selected based on median absolute deviation).
#' }
#'
#' This dataset is derived from the TCGA PAAD cohort (accession: DbGaP phs000178, BroadGDAC: PAAD) and was processed as described in Raphael BJ et al., Cancer Cell, 2017 (PMID: 28810144). The dataset includes selected clinical covariates for potential prognostic and predictive modeling and is designed to facilitate reproducible survival and genomics analyses in pancreatic cancer.
#'
#' @source \doi{10.1016/j.ccell.2017.07.007}
#'
#' @references
#' Raphael BJ, et al. "Integrated genomic characterization of pancreatic ductal adenocarcinoma." Cancer Cell. 2017 Aug 14;32(2):185–203.e13.
"pdac_plain"
