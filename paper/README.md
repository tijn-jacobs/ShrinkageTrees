# ShrinkageTrees — R Journal paper

Source for the R Journal submission describing the
[ShrinkageTrees](https://github.com/tijn-jacobs/ShrinkageTrees) package.

## Reproducing the paper

1. Install the package from source (one directory up):

   ```r
   devtools::install(".")
   ```

2. Generate all figures and numeric outputs:

   ```sh
   Rscript paper/scripts/generate_figures.R
   ```

   This fits `SurvivalBART()`, `HorseTrees()`, and `CausalHorseForest()` on
   the bundled `ovarian` dataset, writes summary/print/C-index text files
   to `paper/outputs/`, and saves figures to `paper/figures/`.

3. Knit the paper:

   ```r
   rmarkdown::render("paper/ShrinkageTrees.Rmd")
   ```

   The `show_output()` helper in the Rmd inlines the text files from step 2,
   so all numbers quoted in the summary/print blocks come directly from the
   fit and cannot drift.

## Directory layout

- `ShrinkageTrees.Rmd` — paper source.
- `RJreferences.bib` — bibliography.
- `scripts/generate_figures.R` — entry point that produces everything in
  `outputs/` and `figures/`.
- `figures/` — PDF figures, including hand-drawn standalone TikZ diagrams
  (`architecture.tex`, `regularisation_landscape.tex`, `tau_learner.tex`,
  `ExampleTree_standalone.pdf`).
- `outputs/` — text files (`print_*.txt`, `summary_*.txt`, `cindex*.txt`)
  inlined into the Rmd.
- `motivation-letter/` — editor-facing cover letter.

## Compiling a standalone figure

Each `figures/*.tex` file is a self-contained `standalone` TikZ document.
Compile in place, for example:

```sh
cd paper/figures
pdflatex architecture.tex
```
