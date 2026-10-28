# NIAID-1418 | Spectral Flow Cytometry: 5L vs 6L

## Overview

This repository contains manuscript-focused HTML reports comparing fluorochrome spectral signatures between **five-laser (5L)** and **six-laser (6L)** instrument configurations across three datasets: **63c (full)**, **83c (filtered)**, and **90c (full)**. The main report displays figures side by side. Each fixation treatment is analyzed independently; the primary comparison is **instrument configuration**, not treatment.

## Original collaborator request (verbatim)

> Thinking on ways to structure the data into figures that will build the manuscript, we would like to first approach the data looking at the impact of the 2 configurations separately on the fluorochromes. When just comparing between 5L and 6L, are there changes in the peaks? How different is the spectral signature of a fluorochrome in the 5L and 6L? Could you run this analysis? The focus would be the instrument configuration, not the treatment, with each treatment as a separate analysis, one independent of the other.

## Additional questions addressed

1. **Peak positions:** Does a fluorochrome's peak detector assignment change between 5L and 6L?
2. **Shared channels:** Do apparent peak shifts persist when considering only detectors active in both configurations?
3. **Spectral shape:** How similar are the 5L and 6L signatures on shared active channels?
4. **Difference magnitude:** Which fluorochromes show the largest configuration-associated changes?
5. **Representative spectra:** What do the overlaid signatures and marked peak positions look like for individual fluorochromes?
6. **Across datasets:** How do these results appear side by side for 63c, filtered 83c, and 90c, without pooling treatment effects?

## Datasets

| Dataset | Included fluorochromes | Comparison |
| --- | --- | --- |
| 63c | Full set | 5L vs 6L, separately by treatment |
| 83c | Filtered; 16 previously flagged fluorochromes excluded | 5L vs 6L, separately by treatment |
| 90c | Full set | 5L vs 6L, separately by treatment |

The four independent treatment analyses are **Unfixed**, **1% PFA**, **4% PFA**, and **Perm**. No treatment is compared with Unfixed in this configuration-focused report.

## What the HTML report shows

- **Peak detector assignments** for 5L and 6L, with shared-channel assessment.
- **Cosine similarity** of spectral signatures on common active detector channels.
- **Normalized RMSE** to quantify differences in spectral shape.
- **Lowest-similarity fluorochromes** for focused review.
- **5L/6L spectral overlays**, including marked maxima and zoomed views where available.
- **Side-by-side figures** for 63c, filtered 83c, and 90c, organized for manuscript figure selection.



## Project

**NIAID-1418** — spectral flow cytometry instrument-configuration analysis for manuscript development.
