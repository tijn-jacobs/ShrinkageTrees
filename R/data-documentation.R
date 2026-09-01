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




#' Semi-synthesised TCGA Ovarian Cancer Dataset
#'
#' Gene expression and clinical covariates for ovarian cancer patients from
#' The Cancer Genome Atlas (TCGA-OV), combined with semi-synthetic
#' survival outcomes and treatment assignment. Real covariates (age,
#' FIGO stage, tumor grade, gene expression) are retained; survival
#' times, event indicator, and treatment assignment are simulated from a
#' known data-generating process so that the true treatment effect is
#' available for validation (see \code{\link{ovarian_truth}}).
#'
#' @format A data frame with 357 rows (patients) and 1004 columns:
#' \itemize{
#'   \item \strong{OS_time}: Numeric. Observed survival time in days (simulated).
#'   \item \strong{OS_event}: Integer. Event indicator (simulated).
#'     1 = event observed, 0 = right-censored.
#'   \item \strong{treatment}: Integer. Simulated treatment assignment.
#'     1 = carboplatin, 0 = cisplatin. Driven primarily by
#'     \code{year_of_diagnosis}, a confounder
#'     (cisplatin era pre-~2000, carboplatin after).
#'   \item \strong{age}: Integer. Age at initial pathologic diagnosis in years.
#'   \item \strong{figo_stage}: Integer. FIGO stage coded as 2 = Stage II,
#'     3 = Stage III, 4 = Stage IV.
#'   \item \strong{tumor_grade}: Integer. Histologic tumor grade coded as
#'     2 = G2, 3 = G3, 4 = G4. Rows with GX (unknown grade) were
#'     excluded.
#'   \item \strong{year_of_diagnosis}: Integer. Year of initial pathologic diagnosis
#'     (approx. 1992--2013). A confounder in the DGP: it drives treatment
#'     assignment and also has a direct prognostic effect.
#'   \item \strong{ENSG...}: Numeric. Log2(TPM + 1) normalised gene expression
#'     levels for 997 Ensembl genes (columns named by versioned Ensembl
#'     gene IDs, e.g. \code{ENSG00000270372.1}). Genes were selected as
#'     the most variable transcripts across TCGA-OV samples, ranked by
#'     median absolute deviation (MAD).
#' }
#'
#' @details
#' Every column except \code{OS_time}, \code{OS_event} and \code{treatment} is
#' a baseline covariate, so the design matrix is
#' \preformatted{
#'   drop <- c("OS_time", "OS_event", "treatment")
#'   X <- as.matrix(ovarian[, setdiff(names(ovarian), drop)])
#' }
#' Select by name rather than position. Earlier versions of this dataset
#' carried outcome-derived \code{left_time} and \code{right_time} columns
#' together with a duplicated \code{year_of_diagnosis}; positional selection
#' silently admitted them to the design matrix, leaking the outcome. They have
#' been removed. To demonstrate the interval-censored interface, build the
#' bounds from the outcome instead:
#' \preformatted{
#'   left  <- ovarian$OS_time
#'   right <- ifelse(ovarian$OS_event == 1, ovarian$OS_time, Inf)
#' }
#'
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
#' # Dimensions: patients x (6 clinical + 2000 gene columns)
#' dim(ovarian)
#'
#' # Survival outcome
#' head(ovarian[, c("OS_time", "OS_event", "treatment")])
#'
#' # KM plot by treatment
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   fit <- survfit(Surv(OS_time, OS_event) ~ treatment, data = ovarian)
#'   plot(fit, col = c("blue", "red"), xlab = "Time (days)", ylab = "Survival")
#'   legend("topright", c("Carboplatin", "Cisplatin"), col = c("blue", "red"), lty = 1)
#' }
"ovarian"


#' Ground-truth quantities for the semi-synthesised ovarian dataset
#'
#' The simulated quantities that correspond to the
#' \code{\link{ovarian}} dataset. Because the ovarian outcomes and
#' treatment assignment are generated from a known data-generating
#' process, the underlying potential outcomes, prognostic function,
#' conditional treatment effect, and propensity score are available for
#' validating estimators of treatment effects under right- and
#' interval-censored survival.
#'
#' @format A data frame with one row per patient in \code{\link{ovarian}}
#'   and the following columns:
#' \describe{
#'   \item{mu}{Numeric. Prognostic function \eqn{\mu(x)}: expected log
#'     survival under the control treatment.}
#'   \item{tau}{Numeric. Conditional average treatment effect
#'     \eqn{\tau(x)} on the log-survival scale. This is the CATE;
#'     \eqn{\exp(\tau)} is the acceleration factor.}
#'   \item{f}{Numeric. The true regression function
#'     \eqn{f(x) = \mu(x) + A\,\tau(x)}, that is
#'     \eqn{E[\log T \mid x, A]} under the treatment actually assigned.
#'     This is what a prediction model estimates when treatment is included
#'     as an ordinary covariate, so it is the target for prediction error.}
#'   \item{propensity}{Numeric. True propensity \eqn{e(x) = P(A = 1 \mid x)}
#'     used to simulate the observed assignment.}
#'   \item{log_time}{Numeric. True (uncensored) survival time on the log
#'     scale, under the treatment actually assigned.}
#'   \item{time}{Numeric. The same on the original scale, in months.}
#'   \item{censoring_time}{Numeric. The censoring time drawn for this
#'     patient, so it is visible why an observation was censored.}
#' }
#'
#' @details
#' The data-generating process is
#' \deqn{\log T = \mu(x) + A\,\tau(x) + \sigma W,}
#' with \eqn{W} standard normal. The residual scale \eqn{\sigma} is stored as
#' the attribute \code{"sigma"}, and the indices of the covariates that enter
#' \eqn{\mu} and \eqn{\tau} as \code{"active_prog"} and \code{"active_tau"}.
#' Attributes are dropped by row subsetting, so read them before subsetting.
#'
#' \code{f} is the target for a model that uses treatment as a covariate. A
#' model fitted without it estimates \eqn{\mu(x) + e(x)\,\tau(x)} instead, by
#' the law of total expectation. Score prediction error against whichever
#' matches the design matrix that was used.
#'
#' @seealso \code{\link{ovarian}} for the observed semi-synthesised data.
#'
#' @examples
#' data(ovarian)
#' data(ovarian_truth)
#' stopifnot(nrow(ovarian) == nrow(ovarian_truth))
#'
#' # True average treatment effect on the log-survival scale
#' mean(ovarian_truth$tau)
#'
#' # log T decomposes as f plus noise
#' stopifnot(all.equal(ovarian_truth$f,
#'                     ovarian_truth$mu + ovarian$treatment * ovarian_truth$tau))
#'
#' # Target for a model fitted without treatment
#' eta <- ovarian_truth$mu + ovarian_truth$propensity * ovarian_truth$tau
#'
#' # Best out-of-sample concordance the process allows: no model can beat this
#' sigma <- attr(ovarian_truth, "sigma")
#' rho   <- sd(eta) / sqrt(var(eta) + sigma^2)
#' 0.5 + asin(rho) / pi
"ovarian_truth"
