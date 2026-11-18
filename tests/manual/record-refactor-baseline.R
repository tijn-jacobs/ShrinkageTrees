# Record a baseline from the code BEFORE the 2.1.1 clean-up, so that
# tests/testthat/test-refactor-equivalence.R can compare the current code
# against it.
#
# Run from the package root, on the same machine and R version that will run
# the tests:
#
#   git worktree add ../ShrinkageTrees-baseline 8d4bc0d
#   Rscript tests/manual/record-refactor-baseline.R ../ShrinkageTrees-baseline
#   Rscript -e 'devtools::test(filter = "refactor")'
#   git worktree remove --force ../ShrinkageTrees-baseline
#
# 8d4bc0d ("update to 2.1.1") is the last commit before the clean-up. Any
# other commit, branch or tag works as the first argument of `git worktree`.
#
# The old package is installed into a temporary library. The package in the
# working tree is not touched and does not need to be installed.
#
# To see a per-scenario table instead of test output, add `--compare`:
#
#   Rscript tests/manual/record-refactor-baseline.R --compare
#
# This loads the CURRENT package with pkgload and prints, for every scenario,
# the largest absolute difference from the recorded baseline.

args      <- commandArgs(trailingOnly = TRUE)
fixture   <- file.path("tests", "testthat", "fixtures", "refactor-baseline.rds")
scenarios <- file.path("tests", "testthat", "helper-refactor-scenarios.R")

if (!file.exists(scenarios))
  stop("Run this script from the package root.")

if (length(args) == 1 && args[1] == "--compare") {

  # ── Compare the current code with the recorded baseline ───────────────────
  if (!file.exists(fixture))
    stop("No baseline found at ", fixture, ". Record one first.")
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Package 'pkgload' is needed for --compare.")

  pkgload::load_all(".", quiet = TRUE)
  source(scenarios)

  baseline <- readRDS(fixture)
  cat("Baseline: ShrinkageTrees ", baseline$package_version, " (",
      baseline$source, "), recorded ", format(baseline$recorded), "\n\n",
      sep = "")

  max_diff <- function(new, old) {
    fields <- intersect(names(new), names(old))
    if (length(fields) == 0) return(NA_real_)
    max(vapply(fields, function(f) {
      if (!identical(dim(new[[f]]), dim(old[[f]])) ||
          length(new[[f]]) != length(old[[f]])) return(Inf)
      max(abs(new[[f]] - old[[f]]))
    }, numeric(1)))
  }

  runs <- refactor_scenarios()
  rows <- lapply(names(runs), function(name) {
    old <- baseline$values[[name]]
    new <- refactor_run(runs[[name]])
    data.frame(
      scenario      = name,
      max_abs_diff  = max_diff(new, old),
      only_in_old   = paste(setdiff(names(old), names(new)), collapse = ","),
      only_in_new   = paste(setdiff(names(new), names(old)), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  print(table, row.names = FALSE, right = FALSE)

  ok <- all(is.finite(table$max_abs_diff) & table$max_abs_diff < 1e-8) &&
    all(table$only_in_old == "") && all(table$only_in_new == "")
  cat("\n", if (ok) "All scenarios match the baseline."
            else "Some scenarios differ from the baseline.", "\n", sep = "")
  quit(status = if (ok) 0 else 1)
}

# ── Record the baseline from an old source tree ─────────────────────────────
if (length(args) != 1 || !dir.exists(args[1]))
  stop("Usage: Rscript tests/manual/record-refactor-baseline.R <old-source-dir>")

old_source <- normalizePath(args[1])

lib <- tempfile("ShrinkageTrees-baseline-lib-")
dir.create(lib)
message("Installing the old package from ", old_source, " ...")
install.packages(old_source, lib = lib, repos = NULL, type = "source",
                 quiet = TRUE)
library(ShrinkageTrees, lib.loc = lib)

source(scenarios)

runs <- refactor_scenarios()
values <- lapply(names(runs), function(name) {
  message("  ", name)
  refactor_run(runs[[name]], baseline = TRUE)
})
names(values) <- names(runs)

dir.create(dirname(fixture), showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(
    values          = values,
    package_version = as.character(utils::packageVersion("ShrinkageTrees")),
    source          = old_source,
    r_version       = R.version.string,
    platform        = R.version$platform,
    recorded        = Sys.time()
  ),
  fixture
)
message("Baseline for ", length(values), " scenarios written to ", fixture)
