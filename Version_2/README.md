# NIAID-1418 Spectral Signature Analysis

Impact of **laser configuration (5L vs 6L)** and **fixation** on the spectral signatures of fluorochromes, measured on antibody-capture bead experiments run on a spectral flow cytometer.

This repository holds the analysis code, the parameterized Quarto report, and the rendered HTML reports for each colour dataset.

---

## Questions addressed

The aims come from the RTB Flow Core slide deck. Each report answers them in order.

| Aim | Question | Reported as |
|---|---|---|
| **Q1a** | Do spectral peaks shift (5L vs 6L, and treatment vs Unfixed)? | Six-category peak-outcome classification, Sony Spectral Viewer validation (63c), per-fluorochrome plots |
| **Q1b** | Is similarity between fluorochromes constant between 5L and 6L? | Pairwise cosine similarity, threshold counts, clustered heatmaps |
| **Q2a** | What is the overall impact of each fixation treatment, and is any worse? | RMSE vs Unfixed, omnibus and pairwise treatment tests |
| **Q2b** | Are particular fluorochromes or fluorochrome families more susceptible? | Susceptibility ranking, family-level test |
| **Q2c** | Does the fixation impact differ between 5L and 6L? | Paired 5L vs 6L fixation effects |

**Treatments:** Unfixed (baseline), 1% PFA, 4% PFA, Perm. **Configurations:** 5L and 6L (the 6L adds the 320 nm dUV laser block).

---

## Datasets and reports

The experiments were run as three colour panels (**63c**, **83c**, **90c**). Each is analysed on its own, and all three are also pooled (**ALL**). Every dataset has its own report.

| Dataset | Fluorochromes | Report | Results folder |
|---|---|---|---|
| 63c | Full | `V2_NIAID-1418_63c_full.html` | `Reports/` |
| 83c | Full (16 flagged dyes kept in) | `V2_NIAID-1418_83c_full.html` | `Reports/` |
| **83c (filtered)** | **16 flagged dyes removed** | `V2_NIAID-1418_83c_filtered.html` | `Reports/` |
| 90c | Full | `V2_NIAID-1418_90c_full.html` | `Reports/` |
| ALL (63c + 83c filtered + 90c) | Full | `V3_NIAID-1416_63_90_83_filtered_slim_down.html` | `Reports/` |

> **Recommended:** use the **83c (filtered)** version. See below.

---

## Full vs filtered, and which to use

### Why there is a filtered version

The PI reviewed the first results and asked us to keep an eye on 16 fluorochromes in the **83c** set. These were flagged as likely bad (signature errors), and in testing they behaved as clear outliers. Because the dataset is a set of single-stain bead signatures, a handful of faulty signatures can distort dataset-wide summaries: RMSE distributions, similarity statistics, family-level susceptibility, and the pooled **ALL** analysis.

To handle this transparently, **every dataset was processed two ways**:

- **Full:** all fluorochromes retained, nothing removed.
- **Filtered:** the 16 flagged fluorochromes removed *before* the analysis. Removal happens at the single point where the common fluorochrome set is defined, so every statistic (Q1a to Q2c) is computed natively on the reduced panel, not patched afterwards.

### Flagged fluorochromes (removed in filtered runs)

`BV605`, `SB645`, `StarBright Violet 670`, `BV650`, `SB702`, `StarBright Violet 710`, `BV711`, `BV750`, `VioBright-515`, `Vio515`, `BB515`, `NB510`, `KIRAVIA Blue`, `AF488`, `FITC`, `NB530`

Names are matched against each set's data. Filtering is only active where a flagged name is present; names absent from a given set are skipped (`allow_unmatched <- TRUE` in the driver).

### Recommendation

**Use the filtered version for interpretation and reporting.** The flagged dyes were identified independently of the results as having signature errors, so including them mixes measurement artefacts with genuine laser and fixation effects. The filtered reports give a cleaner answer to the actual scientific questions.

**Keep the full version for transparency.** It shows exactly what the data look like with nothing excluded, and lets readers see the effect of the exclusion. Where the full and filtered conclusions agree, that is supporting evidence that the findings do not depend on the excluded dyes. Where they differ, the filtered result is the one to trust, and the difference is attributable to the 16 flagged dyes.

---

## Repository layout (TO DO.....)

```

```

`<set>` is one of `63c`, `83c`, `90c`, `ALL`.

---

## Reproducing the analysis

### 1. Run the analysis (`driver_full.r`)

The driver runs `spectral-analysis_V4.R` once per set and mode. It does not edit V4 on disk. For each run it writes a temporary copy with the input workbooks, output folder and family table swapped in, and for filtered runs it also switches on the exclusion list.

Edit the config block at the top:

```r
sets  <- c("63c", "83c", "90c", "ALL")
modes <- c("full", "filtered")
```

then run:

```r
source("driver_full.r")
```

Outputs go to `spectral_analysis_results_<set>/` (full) and `spectral_analysis_results_<set>_filtered/` (filtered). A run summary table is printed at the end, listing any failed runs.

**Inputs**

- `recreated_xlsx/Compiled_matrixes_{5L,6L}_treatments_<set>.xlsx`: spectral matrices per laser configuration and treatment
- `Q2b_fluorochrome_family_annotation_CURATED_3datasets.csv`: curated family table (needed for Q2b)

### 2. Render the reports (`run_qmd.R`)

The report is a single parameterized Quarto document. The parameters are:

| Parameter | Meaning |
|---|---|
| `set` | `63c`, `83c`, `90c` or `ALL` |
| `outdir` | Results folder override (use `spectral_analysis_results_<set>_filtered` for filtered reports) |
| `link_base` | Base GitHub URL for the "View plot" links in the Q1a table |
| `plot_links` | `relevant` (plot that fits the row's category) or `all` |
| `debug_links` | Print Q1a rows that received no plot link |

Example:

```r
library(quarto)

# Full
quarto_render("v2_report_PARAMETERIZED_by_set.qmd",
              output_file = "V2_NIAID-1418_63c_full.html",
              execute_params = list(set = "63c"))

# Filtered
quarto_render("v2_report_PARAMETERIZED_by_set.qmd",
              output_file = "V2_NIAID-1418_83c_filtered.html",
              execute_params = list(set = "83c",
                                    outdir = "spectral_analysis_results_83c_filtered"))
```

`run_qmd.R` renders the full set of reports. The **ALL** reports are large; if the render fails or the HTML is too big, use the lighter `v3_report_PARAMETERIZED_by_set_ALL_safe.qmd`.

**R packages:** `quarto`, `dplyr`, `tidyr`, `ggplot2`, `readr`, `readxl`, `DT`

### Notes

- The Sony Spectral Viewer validation in Q1a uses `Sony-channels_and_peak-analysis-63c_20260921.xlsx` and runs for **63c only**. Other sets show a note instead.
- Q1a plots are hosted as one flat folder per set (`spectral_analysis_results_<set>/Q1a_plots/`) with the original file names, and the report's Q1a table links to them through `link_base`.
- The reports are self-contained HTML (`embed-resources: true`), so they can be opened directly in a browser.

---

## Metrics used

- **RMSE** measures the overall magnitude of spectral change versus Unfixed.
- **Cosine similarity** measures similarity in spectral shape (closely related to the Similarity Index used in full-spectrum cytometry).
- **Maximum detector difference** is the largest single-channel change.
- **Total intensity change** is the overall gain or loss of signal across channels.
- **Peak shift** is movement of the detector containing the maximum.

Full definitions and the six Q1a reporting categories are in the report appendix (A1, A2).

---

## Interactive exploration

An in-house Shiny application is available for interactive inspection of individual fluorochromes, treatments, laser configurations and spectral comparisons. It is an exploratory companion; the conclusions in the reports come from the scripted workflow in this repository.

---

## References

- Jensen HA, Wnek R (2021). Analytical performance of a 25-marker spectral cytometry immune monitoring assay in peripheral blood. *Cytometry A* 99:180–193. https://doi.org/10.1002/cyto.a.24290
- Rajwa B, Roederer M (2025). Lost in translation: Harmonizing terminology and defining mathematical tools for panel optimization. *Cytometry A* 107(12):793–816. https://doi.org/10.1002/cyto.a.24972
- Mage PL, Konecny AJ, Mair F (2026). Measurement and prediction of unmixing-dependent spreading due to collinearity in spectral flow cytometry. *Cytometry A* 109(6):436–455. https://doi.org/10.1002/cyto.a.70044
