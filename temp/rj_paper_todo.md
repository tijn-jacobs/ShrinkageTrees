# R Journal paper — ShrinkageTrees: punch list to finish

## Context

The paper `paper/ShrinkageTrees.Rmd` is a late-stage draft for the R Journal
introducing the `ShrinkageTrees` package. The remaining work is now
primarily: run-the-simulation / run-the-benchmark items that need an R
session, one author decision (Bayesian bootstrap), and a few small
pre-submission chores that need info from Tijn (editor name, stale-file
confirmation).

Everything below is scoped to `paper/` unless noted.

---

## A. Explicit red `\textbf{TODO}` markers still in the Rmd

6. **[ShrinkageTrees.Rmd:451-452]** Replace the `\ldots\%` censoring
   percentages in the interval-censored sim description with real numbers
   (and fix grammar: "Approximately …% are interval-censored, …% are …").
7. **[ShrinkageTrees.Rmd:507]** Replace `...` RMSE / coverage cells in the
   `sim-ic-table` chunk with real results, and run the simulation for real
   (current chunk is `eval=FALSE`). Decide on number of replications and
   record it in the caption.
8. **[ShrinkageTrees.Rmd:526]** Update the concluding paragraph of the
   interval-censored simulation with actual findings; add an
   estimated-vs-true survival curve figure if space allows.
9. **[ShrinkageTrees.Rmd:587-597]** Run the benchmarks and fill in the
   Computational-speed subsection: (i) within-package table (BART vs DART vs
   HorseTrees vs CausalHorseForest), (ii) external comparison against
   `BART::abart()`, (iii) scaling plot over n, p, m (`figures/scaling.pdf`).
   Also fill in CPU model and R version placeholders (`Intel …, R version …`).

---

## B. Numbers that are stale / contradict the actual runs

16. **Placeholder censoring percentages in sim section (line 449)** –
    `\ldots\%` for interval-censored / right-censored / exactly-observed
    shares. Requires running the interval-censored simulation first
    (coupled to A.6).

---

## C. Figures: presence, regeneration, cleanup

20. `figures/ExampleTree_standalone.pdf` is a standalone pre-made diagram.
    Confirm it still renders at the current fig width; if it's from an older
    revision, regenerate or re-save.
22. Remaining figure to add (from §A): `scaling.pdf`; optional
    `sim_ic_curves.pdf` (only if item A.8 calls for it).

---

## F. Bayesian bootstrap — decision point

37. `temp/bayesian_bootstrap_implementation.md` drafts a ~15-line post-hoc
    `bayesian_bootstrap_ate()` function that reweights stored CATE samples
    with Dirichlet weights to produce proper PATE credible intervals.
    The Discussion currently frames this as future work (line 641).
    Decision to make:
    - **(a) Keep as future work** — no package or paper changes.
    - **(b) Ship it** — add the function to the package, rerun the ovarian
      causal block, report both MATE and PATE in §4.2, and delete the
      "future work" paragraph. This strengthens the paper (turns a
      limitation into a feature) at the cost of one package release and a
      rerun.
      Recommend (b) if time allows; the implementation is small and purely
      post-hoc.

---

## G. Pre-submission checks

40. **Motivation letter.** `motivation-letter/motivation-letter.md` is a
    good draft; fill in the editor's name (currently "Professor Tanaka" —
    confirm that's the correct current RJ editor) and regenerate the PDF.
42. **Session info.** The R Journal usually asks for a session info block
    at the end; either embed one or include it in a reproducibility
    appendix.
44. **Stale files in repo.** `ShrinkageTrees-improved.R`,
    `ShrinkageTrees-improved.tex`, `ShrinkageTrees_improved.html`,
    `ShrinkageTrees_oldversion.Rmd`, `initial_checks.log` look like
    leftovers from earlier drafts. Confirm and delete before submission.
    (`.Rhistory` and `.DS_Store` are now covered by `paper/.gitignore`.)

---

## H. Suggested order of attack

1. **Run the interval-censored simulation for real (A.6–A.8 plus B.16):**
   write a separate `scripts/simulation_interval_censored.R`, run with
   e.g. 100 replications, cache results to `outputs/sim_ic.rds`, wire into
   §4.3, and fill in the censoring percentages at line 449.
2. **Run benchmarks (A.9):** separate script
   `scripts/benchmarks.R`, cache to `outputs/benchmarks.rds`, produce
   `figures/scaling.pdf`.
3. **Bayesian bootstrap decision (F.37):** if (b), implement in package,
   rerun §4.2, update Discussion; if (a), leave as-is.
4. **Final read-through and pre-submission (G.40, G.42, G.44).**

## I. Notes from Tijn of other stuff

3. Fix the lay-out of the schematic overview of the architecture (done in
   an earlier pass; re-check once more on the rendered PDF).
