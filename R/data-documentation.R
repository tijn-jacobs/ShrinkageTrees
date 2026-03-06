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


#' TCGA Ovarian Cancer Dataset
#'
#' Gene expression and clinical data for ovarian cancer patients from
#' The Cancer Genome Atlas (TCGA-OV). The dataset contains log2-normalised
#' gene expression profiles alongside overall survival outcomes, treatment
#' assignment, and clinical covariates.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{X}{A numeric matrix of dimensions n x 2000, where rows are patients
#'     and columns are genes. Values are log2(TPM + 1) normalised expression
#'     levels. Genes were selected as the top 2000 most variable genes across
#'     all TCGA-OV samples, ranked by median absolute deviation (MAD).}
#'   \item{clinical}{A data frame with n rows and 7 columns:
#'     \describe{
#'       \item{bcr_patient_barcode}{Character. TCGA patient barcode (12 characters),
#'         e.g. \code{"TCGA-04-1331"}.}
#'       \item{OS_time}{Numeric. Overall survival time in days. Defined as days
#'         to death for deceased patients and days to last follow-up for
#'         censored patients.}
#'       \item{OS_event}{Integer. Overall survival event indicator.
#'         1 = death observed, 0 = censored.}
#'       \item{age}{Integer. Age at initial pathologic diagnosis in years.}
#'       \item{figo_stage}{Integer. FIGO staging score coded as
#'         2 = Stage II, 3 = Stage III, 4 = Stage IV.}
#'       \item{tumor_grade}{Integer. Histologic tumor grade coded as
#'         2 = G2, 3 = G3, 4 = G4. GX (unknown grade) was coded as
#'         \code{NA} and excluded.}
#'       \item{treatment}{Integer. First-line platinum-based chemotherapy.
#'         1 = carboplatin, 0 = cisplatin. Patients who received both
#'         carboplatin and cisplatin were coded as 0 (cisplatin group).
#'         Patients with ambiguous or missing treatment records were excluded.}
#'     }
#'   }
#' }
#' @details
#' RNA-seq data were downloaded from the GDC portal using the
#' \code{TCGAbiolinks} package (STAR - Counts workflow). Expression values
#' were normalised to TPM and log2-transformed as log2(TPM + 1). Genes with
#' median TPM <= 1 across all samples were removed prior to MAD filtering.
#' Clinical data were obtained from the BCR Biotab clinical supplement.
#' Treatment assignment was derived from the drug table
#' (\code{clinical_drug_ov}), restricted to adjuvant (first-line) treatment
#' records. Samples were matched between expression and clinical data using
#' the 12-character TCGA patient barcode.
#'
#' @source \url{https://portal.gdc.cancer.gov/projects/TCGA-OV}
#'
#' @references
#' Cancer Genome Atlas Research Network (2011). Integrated genomic analyses
#' of ovarian carcinoma. \emph{Nature}, 474, 609--615.
#' \doi{10.1038/nature10166}
#'
#' Colaprico, A. et al. (2016). TCGAbiolinks: an R/Bioconductor package for
#' integrative analysis with GDC data. \emph{Nucleic Acids Research}, 44(8).
#' \doi{10.1093/nar/gkv1507}
#'
#' @examples
#' data(ovarian)
#'
#' # Expression matrix
#' dim(ovarian$X)
#'
#' # Survival outcome
#' head(ovarian$clinical[, c("OS_time", "OS_event", "treatment")])
#'
#' # KM plot by treatment
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   fit <- survfit(Surv(OS_time, OS_event) ~ treatment, data = ovarian$clinical)
#'   plot(fit, col = c("blue", "red"), xlab = "Time (days)", ylab = "Survival")
#'   legend("topright", c("Carboplatin", "Cisplatin"), col = c("blue", "red"), lty = 1)
#' }
"ovarian"