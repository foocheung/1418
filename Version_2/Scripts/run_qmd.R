library(quarto)

qmd <- "v2_report_PARAMETERIZED_by_set.qmd"   # your actual filename

# --- Full (4) ---
quarto_render(qmd, output_file = "V2_NIAID-1418_63c_full.html",
              execute_params = list(set = "63c"))
quarto_render(qmd, output_file = "V2_NIAID-1418_83c_full.html",
              execute_params = list(set = "83c"))
quarto_render(qmd, output_file = "V2_NIAID-1418_90c_full.html",
              execute_params = list(set = "90c"))


# --- Filtered (2) — must pass outdir with the _filtered suffix ---
quarto_render(qmd, output_file = "V2_NIAID-1418_83c_filtered.html",
              execute_params = list(set = "83c", outdir = "spectral_analysis_results_83c_filtered"))

##MIGHT FAIL TOO BIG NEED THE OTHER QMD FILE TEST
qmd <- "v3_report_PARAMETERIZED_by_set_ALL_safe.qmd"

quarto_render(qmd, output_file = "V3_NIAID-1418_ALL_filtered.html",
              execute_params = list(set = "ALL", outdir = "spectral_analysis_results_ALL_filtered"))

quarto_render(qmd, output_file = "V2_NIAID-1418_ALL_full.html",
              execute_params = list(set = "ALL"))