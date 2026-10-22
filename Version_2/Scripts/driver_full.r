## run_V4_by_color_v3.R
## Runs spectral-analysis_V4.R once per colour set (63c, 83c, 90c, ALL),
## for BOTH exclusion modes:
##   mode "full"     -> spectral_analysis_results_<set>/
##   mode "filtered" -> spectral_analysis_results_<set>_filtered/
##
## V4 is not modified on disk. For each run the driver writes a temporary copy in
## which FILE_6L / FILE_5L / OUTDIR / FAMILY_FILE are swapped, and for filtered
## runs it also prepends EXCLUSION_MODE and EXCLUDE_LIST so V4's built-in
## exclusion logic activates. V4 removes the flagged dyes at the single point
## where common_fluors / partial_fluors are defined, so the entire pipeline and
## all statistics are computed natively on the reduced panel.
##
## Q2b needs the curated family table for all datasets. Put this file in the same
## folder as spectral-analysis_V4.R. Set xlsx_dir to the folder holding the
## recreated workbooks.

## ---------------------------------------------------------------------------
## Config
## ---------------------------------------------------------------------------
v4_script   <- "v4_full.r"
xlsx_dir    <- "recreated_xlsx"
family_file <- "Q2b_fluorochrome_family_annotation_CURATED_3datasets.csv"
##sets        <- c("63c", "83c", "90c", "ALL")   # drop any you don't want
#sets  <- c("63c", "90c")
sets  <- c("ALL")
modes <- c("filtered")
#modes       <- c("full", "filtered")           # "full", "filtered", or both
stop_on_error <- FALSE                          # FALSE = log failure, continue
##allow_unmatched <- TRUE

## Collaborator-flagged "possible bad" fluorochromes (dropped in filtered runs).
## NOTE corrected spacing on "StarBright Violet 710".
EXCLUDE_LIST <- c(
  "BV605", "SB645", "StarBright Violet 670", "BV650", "SB702",
  "StarBright Violet 710", "BV711", "BV750", "VioBright-515", "Vio515",
  "BB515", "NB510", "KIRAVIA Blue", "AF488", "FITC", "NB530"
)

## If TRUE, a filtered run proceeds even when some EXCLUDE_LIST names are not
## found in that set's data (runs with the matched subset). If FALSE, V4 STOPS
## on any unmatched name so you can fix spelling first. Keep FALSE until the
## first filtered run confirms every name matches.
#allow_unmatched <- FALSE
allow_unmatched <- TRUE
## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------
src <- readLines(v4_script, warn = FALSE)

if (!file.exists(family_file))
  stop("Family table '", family_file, "' not found in ", getwd(),
       ". Run extend_family_table.R first.")

swap <- function(lines, var, value) {
  pat <- sprintf("^\\s*%s\\s*<-\\s*\"[^\"]*\"", var)
  i <- grep(pat, lines)
  if (length(i) != 1) stop("Expected exactly one assignment of ", var, " in ",
                           v4_script, " (found ", length(i), ")")
  lines[i] <- sprintf("%s <- \"%s\"", var, value)
  lines
}

## Build the header injected at the TOP of the temp copy for filtered runs.
## These globals are read by V4 (it only sets defaults if they don't already
## exist), so defining them first makes V4 run in filtered mode.
make_excl_header <- function(exclude_vec, allow_unmatched) {
  c(
    'EXCLUSION_MODE <- "filtered"',
    paste0('EXCLUDE_LIST <- c(',
           paste(sprintf('"%s"', exclude_vec), collapse = ", "), ')'),
    paste0('EXCL_ALLOW_UNMATCHED <- ', if (isTRUE(allow_unmatched)) 'TRUE' else 'FALSE')
  )
}

status <- expand.grid(Set = sets, Mode = modes, stringsAsFactors = FALSE)
status$Result <- NA_character_
status$Outdir <- NA_character_

## ---------------------------------------------------------------------------
## Main loop
## ---------------------------------------------------------------------------
for (r in seq_len(nrow(status))) {
  s    <- status$Set[r]
  mode <- status$Mode[r]
  f5   <- file.path(xlsx_dir, sprintf("Compiled_matrixes_5L_treatments_%s.xlsx", s))
  f6   <- file.path(xlsx_dir, sprintf("Compiled_matrixes_6L_treatments_%s.xlsx", s))
  out  <- if (mode == "filtered")
    sprintf("spectral_analysis_results_%s_filtered", s)
  else
    sprintf("spectral_analysis_results_%s", s)

  if (!file.exists(f5) || !file.exists(f6)) {
    message("[", s, "/", mode, "] workbook(s) not found: ", f5, " / ", f6)
    status$Result[r] <- "missing input"; status$Outdir[r] <- out; next
  }

  lines <- src
  lines <- swap(lines, "FILE_6L",     f6)
  lines <- swap(lines, "FILE_5L",     f5)
  lines <- swap(lines, "OUTDIR",      out)
  lines <- swap(lines, "FAMILY_FILE", family_file)

  if (mode == "filtered") {
    lines <- c(make_excl_header(EXCLUDE_LIST, allow_unmatched), lines)
  }

  tmp <- tempfile(pattern = paste0("V4_", s, "_", mode, "_"), fileext = ".R")
  writeLines(lines, tmp)

  message("\n################ ", s, " / ", mode, " -> ", out, " ################")
  res <- try(source(tmp, local = new.env(parent = globalenv()), echo = FALSE),
             silent = FALSE)
  unlink(tmp)

  if (inherits(res, "try-error")) {
    status$Result[r] <- "FAILED"
    message("[", s, "/", mode, "] FAILED: ",
            conditionMessage(attr(res, "condition")))
    if (stop_on_error) stop("Stopped at ", s, "/", mode)
  } else {
    status$Result[r] <- "ok"
  }
  status$Outdir[r] <- out
}

cat("\n==================== RUN SUMMARY ====================\n")
print(status, row.names = FALSE)
cat("\nFiltered runs that FAILED with an unmatched-name error mean a flagged\n")
cat("name did not match any dye in that set. Fix the spelling in EXCLUDE_LIST\n")
cat("(check the printed present-identity sample), or set allow_unmatched <- TRUE\n")
cat("to run those sets with the matched subset only.\n")
