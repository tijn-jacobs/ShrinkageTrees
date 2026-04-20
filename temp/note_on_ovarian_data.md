# TCGA-OV Analysis Dataset

## Background

The Cancer Genome Atlas (TCGA) Ovarian Serous Cystadenocarcinoma
(**TCGA-OV**) project profiled high-grade serous ovarian cancer across
molecular, clinical, and outcome dimensions. This dataset combines
gene expression, clinical covariates, treatment information, and
survival outcomes into a single patient-level analysis table for
downstream modeling (e.g. Cox regression, penalized survival models,
clustering).

All data were pulled from the GDC Data Portal via the Bioconductor
package `TCGAbiolinks`:

- **Expression**: RNA-seq STAR counts (`data.category = "Transcriptome Profiling"`),
  restricted to primary tumor samples, normalized to log2-CPM
  (without `edgeR`, using a stabilized pseudocount).
- **Clinical**: indexed clinical fields via `GDCquery_clinic()`.
- **Treatment**: BCR XML clinical files, drug-level records
  (`GDCprepare_clinic(..., clinical.info = "drug")`), harmonized by
  regex matching of free-text drug names into standardized categories.

## Cohort construction

Starting from 427 TCGA-OV patients with primary tumor RNA-seq, we
retained patients with:

- A usable survival outcome (`vital_status` ∈ {Alive, Dead} and a
  non-missing, positive follow-up time).
- A known treatment regimen (patients with blank drug records or
  trial-name-only entries were dropped).

The resulting analysis table contains one row per patient. Expression
values are restricted to the **top 1,000 genes by median absolute
deviation (MAD)** computed on log2-CPM values, after filtering low-
expression genes (≥10 counts in ≥20% of samples).

## Variables

### Identifier

- `submitter_id` — TCGA patient barcode (e.g. `TCGA-04-1331`).

### Survival outcome

- `time` — follow-up time in days. Equals `days_to_death` for patients
  who died; `days_to_last_follow_up` for those still alive.
- `delta` — event indicator: `1` = death observed, `0` = censored.

### Clinical covariates

- `age_at_index` — age in years at study entry.
- `figo_main` — FIGO stage collapsed to main stages (I / II / III / IV).
  Sub-stages (e.g. IIIA, IIIC) are merged into their main stage.
- `tumor_grade` — histologic grade (G1–G4, GB, GX).
- `race` — self-reported race.
- `year_of_diagnosis` — year of initial diagnosis. Important as a
  proxy for treatment era (supportive care and surgical practice
  changed substantially over TCGA's enrollment window).
- `prior_treatment` — whether the patient had received prior treatment
  before study entry (Yes / No).

### Treatment variables

Derived from the BCR XML drug table. Each patient may contribute
multiple drug records; these are aggregated into patient-level flags.

- `regimen` — coarse regimen category:
  - `platinum + taxane` — standard first-line HGSOC therapy.
  - `platinum only` — platinum without a taxane.
  - `taxane only` — taxane without a platinum.
  - `other` — non-platinum, non-taxane regimens (rare).
- `platinum_type` — platinum agent(s) received:
  `carbo only`, `cis only`, `carbo + cis`, `other platinum`, `no platinum`.
- `treatment` — simplified binary contrast:
  - `carbo only` — received carboplatin, no cisplatin.
  - `carbo and/or cis` — received cisplatin (with or without carboplatin).
- `any_carbo`, `any_cis`, `any_pacli`, `any_doce`, `any_doxo`, `any_gem`,
  `any_topo`, `any_bev` — logical flags for exposure to individual drugs:
  carboplatin, cisplatin, paclitaxel, docetaxel, doxorubicin, gemcitabine,
  topotecan, bevacizumab.
- `n_drug_records` — number of drug records on file for the patient
  (a rough proxy for lines of therapy / treatment intensity).

### Gene expression

- 1,000 columns named by Ensembl gene ID (e.g. `ENSG00000134184.13`),
  each containing the patient's **log2-CPM** expression for that gene.
- Genes were selected as the top 1,000 by MAD across samples after
  low-expression filtering, giving a compact set of the most variable
  transcripts in the cohort.

## Caveats

- **Treatment is nearly constant.** ~95% of patients received
  platinum + taxane (mostly carboplatin + paclitaxel), so `regimen`
  has limited discriminative power for survival modeling.
- **Cisplatin vs. carboplatin** is strongly confounded by year of
  diagnosis; cisplatin-treated patients were enrolled earlier. Any
  comparison of these subgroups should adjust for `year_of_diagnosis`.
- **Residual disease** is a key prognostic variable in HGSOC but is
  ~97% missing in the GDC clinical data and was therefore excluded.
  If needed, it can be recovered from the Bell et al. 2011 / Liu
  et al. 2018 curated clinical resources.
- **Ensembl IDs with version suffixes** are used as gene identifiers;
  a mapping to HGNC symbols is available via `rowData()` of the
  original SummarizedExperiment.

## File

- `tcga_ov_flat_cc.rds` — analysis-ready data frame, one row per
  patient, ~1,020 columns.

Load with:

    flat_cc <- readRDS("tcga_ov_flat_cc.rds")

## Survival outcome summaries

### By treatment group

| Treatment        |   n | Events | Censored | Event rate | Median time (days) | IQR      | Max time |
| ---------------- | --: | -----: | -------: | ---------: | -----------------: | -------- | -------: |
| carbo only       | 258 |    155 |      103 |       0.60 |                923 | 516–1497 |     4424 |
| carbo and/or cis | 125 |     77 |       48 |       0.62 |               1442 | 914–2218 |     4665 |

_Note: "Median time" here is the median observation time (ignoring censoring), not the median survival time._

### Kaplan–Meier median survival

| Treatment        |   n | Events | Median OS (days) | 95% CI    |
| ---------------- | --: | -----: | ---------------: | --------- |
| carbo only       | 258 |    155 |             1329 | 1158–1484 |
| carbo and/or cis | 125 |     77 |             1784 | 1484–2148 |

### Survival at clinical landmarks

| Treatment        | 1 year           | 3 years          | 5 years          |
| ---------------- | ---------------- | ---------------- | ---------------- |
| carbo only       | 0.93 (0.90–0.97) | 0.60 (0.54–0.68) | 0.29 (0.23–0.36) |
| carbo and/or cis | 0.97 (0.95–1.00) | 0.79 (0.71–0.87) | 0.49 (0.40–0.60) |

Patients in the cisplatin-containing group show substantially higher
survival at every landmark. This is consistent with the earlier
observation that these patients were enrolled earlier (cisplatin era)
and therefore have longer accumulated follow-up — a confound, not
necessarily a drug effect.

## AFT model (Weibull)

Fitted on `time ~ treatment + age_at_index + figo_main + tumor_grade + year_of_diagnosis`.
Coefficients are on the log-time scale; `exp(coef)` is a **time
ratio** (> 1 = longer survival). n = 381 (2 patients dropped for
missing covariates).

| Term                         |                      Coef |    SE | Time ratio | 95% CI    |      p |
| ---------------------------- | ------------------------: | ----: | ---------: | --------- | -----: |
| treatment (carbo and/or cis) |                     0.351 | 0.093 |       1.42 | 1.18–1.71 | 0.0002 |
| age_at_index (per year)      |                    -0.008 | 0.004 |       0.99 | 0.98–1.00 |  0.046 |
| figo_main II (vs I)          |                     0.408 | 0.271 |       1.50 | 0.88–2.56 |   0.13 |
| figo_main III (vs I)         |                     0.090 | 0.115 |       1.09 | 0.87–1.37 |   0.43 |
| figo_main IV (vs I)          |                        NA |    NA |         NA | NA        |     NA |
| year_of_diagnosis (per year) |                     0.031 | 0.011 |       1.03 | 1.01–1.06 |  0.006 |
| tumor_grade                  | unstable — see note below |

**Model fit:** Weibull scale = 0.647. Likelihood-ratio χ² = 33.3 on 11 df, p = 0.00047.

### Interpretation

- **Treatment** — receiving cisplatin (alone or with carboplatin) is
  associated with a **42% longer survival time** vs. carboplatin alone
  (time ratio 1.42, 95% CI 1.18–1.71). This is almost certainly
  driven by era confounding rather than a true drug effect: cisplatin
  patients were enrolled earlier and have had more opportunity to be
  observed as long survivors. The naive comparison should not be read
  as evidence that cisplatin is superior to carboplatin.
- **Age** — each additional year of age at diagnosis is associated
  with ~0.8% shorter survival time (time ratio 0.99, p = 0.046).
- **Year of diagnosis** — each later year of enrollment is associated
  with ~3% longer observed survival (time ratio 1.03, p = 0.006). This
  is in part a real improvement in outcomes over the TCGA enrollment
  window, and in part a reflection that later-enrolled patients have
  shorter possible follow-up (right-censoring pattern).
- **FIGO stage** — none of the non-reference levels reach significance,
  and stage IV collapses to `NA` because it's collinear with the other
  predictors once year/age/treatment are included (only ~1 Stage I
  patient in the cohort, so the reference category is essentially
  empty — the model is over-parameterized).
- **Tumor grade** — all grade coefficients have enormous standard
  errors (~2200) and effectively zero time ratios. This is a
  **numerical degeneracy**, not a real effect: the reference level
  (G1) probably has very few patients, so the model can't separate
  grade from the intercept. Drop or collapse grade before refitting.

### Suggested refit

```r
# Collapse sparse categories to stabilize the model
flat_cc <- flat_cc |>
  mutate(
    figo_late = factor(
      ifelse(figo_main %in% c("III","IV"), "III/IV", "I/II"),
      levels = c("I/II","III/IV")
    ),
    grade_hi  = factor(
      ifelse(tumor_grade %in% c("G3","G4"), "G3/G4", "G1/G2"),
      levels = c("G1/G2","G3/G4")
    )
  )

aft2 <- survreg(
  Surv(time, delta) ~ treatment + age_at_index + figo_late +
                      grade_hi + year_of_diagnosis,
  data = flat_cc,
  dist = "weibull"
)
summary(aft2)
```

This dichotomizes the two sparsely-populated categoricals into
clinically sensible high/low contrasts, which should eliminate the
`NA` and the 2000-sized standard errors.
