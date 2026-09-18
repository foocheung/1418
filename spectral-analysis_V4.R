# ============================================================
# FULL 5L vs 6L SPECTRAL-SIGNATURE ANALYSIS  (merged version)
#
# Implements the NIAID-1418 QMD plan:
#   Q1a  peak behavior / peak shift, evaluated comprehensively:
#          A) 5L vs 6L within each treatment (configuration comparison)
#          B) each fixation treatment vs Unfixed within 5L and within 6L
#          C) whether the treatment-associated peak response is the same in 5L and 6L
#          - peak calls retain the second-highest channel so near-ties are visible
#          - expected/reference peak location remains a separate comparison
#   Q1b  fluorochrome-pair distinguishability and 5L vs 6L similarity changes (Pearson/cosine)
#   Q2a  treatment-vs-Unfixed fixation effects (metrics kept separate)
#   Q2b  fluorochrome susceptibility rankings (per metric, no combined score)
#   Q2c  configuration-dependent fixation effect (6L vs 5L)
#
# KEY DATA FACT (verified against the supplied workbooks):
#   The 5L file keeps the 35 deep-UV (320 nm) detector rows physically present
#   but ZERO-FILLED, because the UV laser is off in the 5L configuration.
#   Treating those zeros as real 5L signal would (a) inflate cosine similarity
#   via matched zeros, (b) distort RMSE / max-diff / intensity, and (c) drag
#   UV-peaking dyes to a spurious 5L peak. We therefore:
#     - detect all-zero ("inactive") channels from the data (not by hard-coding
#       the 320 prefix), so the logic also works if the real 5L file deletes
#       the rows or drops a different laser;
#     - run cross-config comparisons (Q1a/Q1b/Q2c) on the LIVE SHARED channels;
#     - flag dyes whose true 6L peak is in the UV block (peak lost at 5L) so
#       they are not misread as "peak shifts".
# ============================================================

# -------------------------
# 0. Packages
# -------------------------
required <- c("readxl", "dplyr", "tidyr", "ggplot2", "purrr", "stringr", "readr",
              "ggrepel", "forcats", "scales")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(readxl);  library(dplyr);   library(tidyr)
  library(ggplot2); library(purrr);   library(stringr); library(readr)
  library(ggrepel);  library(forcats); library(scales)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -------------------------
# 1. Files / settings
# -------------------------
FILE_6L <- "./Compiled_matrixes_6L_treatments.xlsx"
FILE_5L <- "Compiled_matrixes_5L_treatments.xlsx"

# fall back to a Downloads path if the local copies are absent
alt6 <- "./RE__NIAID-1418__New_Request_Submitted/Compiled_matrixes_6L_treatments.xlsx"
alt5 <- "./RE__NIAID-1418__New_Request_Submitted/Compiled_matrixes_5L_treatments.xlsx"
if (!file.exists(FILE_6L) && file.exists(alt6)) FILE_6L <- alt6
if (!file.exists(FILE_5L) && file.exists(alt5)) FILE_5L <- alt5

OUTDIR     <- "spectral_analysis_results2"

# --------------------------------------------------------------------------
# NORMALIZATION (changed in V3 -- see 00_scaling_audit.csv)
# --------------------------------------------------------------------------
# The supplied workbooks are NOT on a common scale. Every 5L sheet and the
# 6L Perm sheet are max-normalized per dye (column max == 1), but 6L Unfixed,
# 6L 1PFA and 6L 4PFA leave 88/145 dye columns on an arbitrary raw scale with
# maxima up to ~1.5e3. Running the magnitude metrics on "Raw values" therefore
# compares a normalized spectrum against an unnormalized one in 6L, which makes
# Perm look enormously different from Unfixed for those 88 dyes purely as a
# scaling artefact, and makes RMSE_5L vs RMSE_6L (Q2c) incommensurable.
#
# Q1a (argmax) and Q1b (cosine/Pearson) are invariant to a positive per-dye
# scalar and are unaffected. Q2a/Q2b/Q2c are NOT. V3 therefore renormalizes
# every spectrum within the channel set actually used for each comparison.
#   "Max normalize each fluorochrome"  - peak height = 1 (default; robust, but
#                                        anchored on one possibly noisy detector)
#   "Area normalize each fluorochrome" - sum of positive signal = 1 (less
#                                        sensitive to a single detector)
#   "Raw values"                       - ONLY valid if the audit reports a
#                                        consistent native scale
PREPROCESS <- "Max normalize each fluorochrome"

# Sensitivity: also recompute Q2a/Q2c under this second normalization and write
# a concordance table. Set to NA to skip.
PREPROCESS_SENSITIVITY <- "Area normalize each fluorochrome"

NEAR_TIE_RATIO <- 0.99       # second/first ratio at which a maximum is "near tied"
Q1B_REFERENCE_COSINE <- 0.90 # literature reference cutoff for Q1b

# Optional user-curated fluorochrome family table for Q2b. If absent, the script
# writes a heuristic annotation you can edit and feed back in via this path.
FAMILY_FILE <- "Q2b_fluorochrome_family_annotation_CURATED.csv"

for (d in c("", "Q1a_peak_shift", "Q1b_pairwise_similarity", "Q2_fixation", "figures"))
  dir.create(file.path(OUTDIR, d), showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 2. Helpers
# -------------------------
treatment_name <- function(sheet) {
  x <- sub("^\\s*[65]L\\s*", "", sheet)
  x <- sub("\\s*Compil.*$", "", x)
  dplyr::recode(trimws(x),
    "Unfix" = "Unfixed", "1PFA" = "1% PFA", "4PFA" = "4% PFA",
    "Perm" = "Perm", .default = trimws(x))
}

preprocess_matrix <- function(mat, method = "Raw values") {
  out <- mat
  if (method == "Max normalize each fluorochrome") {
    for (j in seq_len(ncol(out))) {
      m <- max(abs(out[, j]), na.rm = TRUE)
      if (is.finite(m) && m > 0) out[, j] <- out[, j] / m
    }
  } else if (method == "Area normalize each fluorochrome") {
    # Normalize on total positive signal rather than a single peak detector.
    for (j in seq_len(ncol(out))) {
      x <- out[, j]
      a <- sum(x[is.finite(x) & x > 0], na.rm = TRUE)
      if (is.finite(a) && a > 0) out[, j] <- out[, j] / a
    }
  } else if (method == "Z-score each fluorochrome") {
    for (j in seq_len(ncol(out))) {
      x <- out[, j]; s <- sd(x, na.rm = TRUE); mu <- mean(x, na.rm = TRUE)
      out[, j] <- if (is.finite(s) && s > 0) (x - mu) / s else 0
    }
  }
  out
}

parse_sheet <- function(path, sheet, config) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE,
                    .name_repair = "minimal", trim_ws = FALSE)

  col_d <- trimws(as.character(raw[[4]]))
  header_row <- which(col_d == "Channels")[1]
  if (is.na(header_row)) stop("Could not find 'Channels' header in sheet: ", sheet)

  raw_names  <- trimws(as.character(unlist(raw[header_row, 5:ncol(raw), drop = FALSE])))
  keep_cols  <- which(!is.na(raw_names) & raw_names != "") + 4
  raw_names  <- raw_names[!is.na(raw_names) & raw_names != ""]

  # Normalize dye names so the same dye matches across files despite data-entry
  # inconsistencies. In these workbooks the only real difference is CASE
  # (e.g. "NY730" vs "ny730"); whitespace collapse guards against stray/interior
  # spaces. Bracket removal is a harmless safeguard. The FIRST-seen original
  # spelling is kept as the human-readable display name.
  clean_name <- function(x) {
    x <- gsub("[\\[\\]]", "", x, perl = TRUE)  # safeguard: drop any [ or ]
    x <- gsub("\\s+", " ", x)                   # collapse whitespace runs
    trimws(x)
  }
  base_names <- clean_name(raw_names)            # display form (brackets removed)
  match_key  <- tolower(base_names)              # case-insensitive matching key

  # preserve genuine repeated dyes explicitly, using the normalized key
  occurrence <- ave(seq_along(match_key), match_key, FUN = seq_along)
  fluor_id <- ifelse(ave(seq_along(match_key), match_key, FUN = length) > 1,
                     paste0(match_key, " [rep ", occurrence, "]"), match_key)

  candidate_rows <- seq.int(header_row + 1, nrow(raw))
  channel_num <- suppressWarnings(as.numeric(as.character(raw[[4]][candidate_rows])))
  rows <- candidate_rows[!is.na(channel_num) & channel_num >= 1 & channel_num <= 320]

  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(keep_cols))
  for (j in seq_along(keep_cols))
    mat[, j] <- suppressWarnings(as.numeric(as.character(raw[[keep_cols[j]]][rows])))
  colnames(mat) <- fluor_id                      # matched on normalized key

  channel_B <- trimws(as.character(raw[[2]][rows]))
  channels <- tibble(
    Channel_A     = trimws(as.character(raw[[1]][rows])),
    Channel_B     = channel_B,
    Channel_C     = trimws(as.character(raw[[3]][rows])),
    ChannelNumber = suppressWarnings(as.numeric(as.character(raw[[4]][rows]))),
    ChannelKey    = channel_B,
    Laser_nm      = suppressWarnings(as.numeric(sub("-.*$", "", channel_B))),
    Detector_nm   = suppressWarnings(as.numeric(sub("^.*-", "", channel_B)))
  )

  active <- apply(mat, 1, function(z) { z <- z[is.finite(z)]; length(z) > 0 && any(abs(z) > 0) })

  list(config = config, treatment = treatment_name(sheet), sheet = sheet,
       matrix = mat, channels = channels, active = active,
       fluor_table = tibble(Raw_name = raw_names, Fluorochrome = base_names,
                            Match_key = match_key, Occurrence = occurrence,
                            Fluor_ID = fluor_id))
}

load_workbook <- function(path, config) {
  if (!file.exists(path)) stop("File not found: ", path)
  sheets <- excel_sheets(path)
  out <- lapply(sheets, function(s) parse_sheet(path, s, config))
  names(out) <- vapply(out, `[[`, character(1), "treatment")
  out
}

# restrict a condition to a set of channel keys and fluorochromes
restrict_profile <- function(p, keys, fluor_ids) {
  idx <- match(keys, p$channels$ChannelKey)
  if (anyNA(idx)) stop("Requested channels absent in ", p$config, " ", p$treatment)
  m <- preprocess_matrix(p$matrix[idx, , drop = FALSE], PREPROCESS)
  list(matrix = m[, fluor_ids, drop = FALSE], channels = p$channels[idx, , drop = FALSE])
}

cosine <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]
  if (!length(a)) return(NA_real_)
  den <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (!is.finite(den) || den == 0) return(NA_real_) else sum(a * b) / den
}
rmse <- function(a, b) { ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) NA_real_ else sqrt(mean((a[ok] - b[ok])^2)) }
safe_cor <- function(a, b, method) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || sd(a[ok]) == 0 || sd(b[ok]) == 0) NA_real_
  else cor(a[ok], b[ok], method = method)
}

peak_info <- function(values, channels) {
  ok <- is.finite(values)
  if (!any(ok)) return(tibble(Peak_ChannelKey = NA_character_, Peak_ChannelNumber = NA_real_,
                              Peak_Laser_nm = NA_real_, Peak_Detector_nm = NA_real_,
                              Peak_Value = NA_real_))
  ii <- which(ok)[which.max(values[ok])]
  tibble(Peak_ChannelKey = channels$ChannelKey[ii],
         Peak_ChannelNumber = channels$ChannelNumber[ii],
         Peak_Laser_nm = channels$Laser_nm[ii],
         Peak_Detector_nm = channels$Detector_nm[ii],
         Peak_Value = values[ii])
}

similarity_matrix <- function(x, method) {  # x: fluorochromes x channels
  if (method == "pearson")
    return(cor(t(x), use = "pairwise.complete.obs", method = method))

  # V3: previously non-finite values were replaced with 0 before the crossproduct,
  # which silently treats a missing detector as a measured zero and inflates
  # cosine via matched zeros -- exactly the failure mode the header warns about
  # for the 5L dUV block. Compute pairwise-complete instead, matching cosine().
  n <- nrow(x)
  s <- matrix(NA_real_, n, n, dimnames = list(rownames(x), rownames(x)))
  for (i in seq_len(n)) {
    s[i, i] <- 1
    if (i < n) for (j in seq.int(i + 1, n)) {
      v <- cosine(x[i, ], x[j, ])
      s[i, j] <- v; s[j, i] <- v
    }
  }
  s
}
pair_frame <- function(sim, treatment, metric, config) {
  idx <- which(upper.tri(sim), arr.ind = TRUE)
  tibble(Treatment = treatment, Metric = metric, Config = config,
         Fluor_A = rownames(sim)[idx[, 1]], Fluor_B = rownames(sim)[idx[, 2]],
         Similarity = sim[idx])
}

# -------------------------
# 3. Load BOTH workbooks
# -------------------------
data6 <- load_workbook(FILE_6L, "6L")
data5 <- load_workbook(FILE_5L, "5L")

common_treatments <- intersect(names(data6), names(data5))
if (!length(common_treatments)) stop("No matching treatments between 5L and 6L.")

# Fluorochromes present in every treatment/config, matched on the NORMALIZED
# name key (whitespace/case standardized in parse_sheet).
all_profiles  <- c(data6[common_treatments], data5[common_treatments])

# ---------------------------------------------------------------------------
# Coverage-based cleanup. Every sheet has the same COLUMN COUNT, but the 5L
# 1PFA sheet swaps a real dye (pe-fire700) for a stray one-letter entry ("v"),
# so the sheets do not share identical dye sets. We therefore:
#   - DROP entries that appear in only one sheet (stray/erroneous cells like
#     "v") — they can't be compared and aren't real measurements;
#   - KEEP dyes missing from just a sheet or two (like pe-fire700, absent only
#     from 5L 1PFA) and log them, so a real dye is not silently discarded by a
#     strict all-sheet intersection.
# We do NOT filter by name length: "PE" is a legitimate 2-character dye.
# ---------------------------------------------------------------------------
n_sheets <- length(all_profiles)
all_ids  <- unique(unlist(lapply(all_profiles, function(p) colnames(p$matrix))))
coverage <- vapply(all_ids, function(id)
  sum(vapply(all_profiles, function(p) id %in% colnames(p$matrix), logical(1))),
  integer(1))

stray_ids <- all_ids[coverage <= 1]                 # present in <=1 sheet -> drop
partial_ids <- all_ids[coverage > 1 & coverage < n_sheets]  # incomplete -> keep+log
full_ids  <- all_ids[coverage == n_sheets]          # present everywhere

if (length(stray_ids))
  message("Dropping stray non-dye entries (present in <=1 sheet): ",
          paste(stray_ids, collapse = ", "))
if (length(partial_ids)) {
  message("Dyes NOT present in all sheets (kept where available, excluded from ",
          "comparisons that need a missing sheet):")
  for (id in partial_ids) {
    where <- vapply(all_profiles, function(p) id %in% colnames(p$matrix), logical(1))
    message(sprintf("  %-22s in %d/%d sheets; missing: %s", id, sum(where), n_sheets,
                    paste(vapply(all_profiles[!where],
                                 function(p) paste(p$config, p$treatment), character(1)),
                          collapse = ", ")))
  }
}

# common_fluors = dyes present in all sheets (used for the strict cross-config
# comparisons). partial dyes are handled per-analysis where noted.
common_fluors <- full_ids
partial_fluors <- partial_ids   # available downstream if a per-treatment analysis wants them

# ---- exclusions audit: every dye NOT in the strict overlap set, with reason ----
sheet_label <- function(p) paste(p$config, p$treatment)
present_in  <- function(id) vapply(all_profiles, function(p) id %in% colnames(p$matrix), logical(1))
excluded_ids <- c(stray_ids, partial_ids)

exclusions <- if (length(excluded_ids)) {
  bind_rows(lapply(excluded_ids, function(id) {
    where   <- present_in(id)
    reason  <- if (id %in% stray_ids)
                 "Stray/non-dye entry (present in <=1 sheet); dropped entirely"
               else
                 "Real dye missing from one or more sheets; excluded from strict overlap set"
    tibble(Fluor_ID = id,
           N_sheets_present = sum(where),
           N_sheets_total   = n_sheets,
           Present_in = paste(vapply(all_profiles[where],  sheet_label, character(1)), collapse = "; "),
           Missing_from = paste(vapply(all_profiles[!where], sheet_label, character(1)), collapse = "; "),
           Reason = reason,
           Action = if (id %in% stray_ids) "Removed" else "Excluded from overlap analyses (kept in raw data)")
  }))
} else {
  tibble(Fluor_ID = character(), N_sheets_present = integer(),
         N_sheets_total = integer(), Present_in = character(),
         Missing_from = character(), Reason = character(), Action = character())
}
write_csv(exclusions, file.path(OUTDIR, "00_excluded_fluorochromes.csv"))

n_cols_each <- vapply(all_profiles, function(p) ncol(p$matrix), integer(1))
cat("Treatments:", paste(common_treatments, collapse = ", "), "\n")
cat(sprintf("Dyes: %d in all %d sheets | %d partial | %d stray dropped (of %d columns/sheet)\n",
            length(common_fluors), n_sheets, length(partial_fluors),
            length(stray_ids), min(n_cols_each)))
if (nrow(exclusions))
  cat(sprintf("Excluded %d dye(s) logged to 00_excluded_fluorochromes.csv\n",
              nrow(exclusions)))

# data audit
audit <- bind_rows(lapply(all_profiles, function(p) tibble(
  Config = p$config, Treatment = p$treatment, Sheet = p$sheet,
  Matrix_rows = nrow(p$matrix), Fluorochrome_columns = ncol(p$matrix),
  Active_channels = sum(p$active), Inactive_zero_channels = sum(!p$active))))
write_csv(audit, file.path(OUTDIR, "00_data_audit.csv"))

# active/live channel sets (valid across every treatment within a config)
active6 <- Reduce(intersect, lapply(data6[common_treatments], \(p) p$channels$ChannelKey[p$active]))
active5 <- Reduce(intersect, lapply(data5[common_treatments], \(p) p$channels$ChannelKey[p$active]))
shared_active <- intersect(active6, active5)

# UV block: live at 6L, inactive (zeroed) at 5L -> the detectors lost by dropping the laser
uv_lost_keys <- intersect(active6, setdiff(data5[[common_treatments[1]]]$channels$ChannelKey, active5))

channel_audit <- tibble(
  Set = c("6L active", "5L active", "Shared active (used for Q1/Q2c)",
          "Live at 6L, zeroed at 5L (UV lost)"),
  N_channels = c(length(active6), length(active5), length(shared_active), length(uv_lost_keys)))
write_csv(channel_audit, file.path(OUTDIR, "00_channel_audit.csv"))
print(channel_audit)

# ---------------------------------------------------------------------------
# 3b. SCALING AUDIT  (new in V3)
# ---------------------------------------------------------------------------
# Establish empirically whether the native values are on a common scale before
# any magnitude metric (RMSE, max detector difference, total intensity) is
# computed. This is not cosmetic: in the supplied workbooks the 5L sheets and
# the 6L Perm sheet are max-normalized per dye while 6L Unfixed / 1PFA / 4PFA
# are not, so a raw-value Q2a would report a scaling artefact as a fixation
# effect. The audit is written for every run so the assumption is always
# visible in the outputs rather than assumed.
scaling_audit <- bind_rows(lapply(all_profiles, function(p) {
  m <- p$matrix
  col_max <- apply(m, 2, function(z) { z <- z[is.finite(z)]; if (!length(z)) NA_real_ else max(z) })
  col_min <- apply(m, 2, function(z) { z <- z[is.finite(z)]; if (!length(z)) NA_real_ else min(z) })
  tibble(
    Config = p$config, Treatment = p$treatment, Sheet = p$sheet,
    N_dye_columns = ncol(m),
    N_max_normalized_to_1 = sum(abs(col_max - 1) < 1e-4, na.rm = TRUE),
    N_not_normalized = sum(abs(col_max - 1) >= 1e-4, na.rm = TRUE),
    Smallest_column_max = suppressWarnings(min(col_max, na.rm = TRUE)),
    Largest_column_max  = suppressWarnings(max(col_max, na.rm = TRUE)),
    Most_negative_value = suppressWarnings(min(col_min, na.rm = TRUE))
  )
}))
scaling_audit <- scaling_audit |>
  mutate(Native_scale = case_when(
    N_not_normalized == 0 ~ "Max-normalized per dye (max = 1)",
    N_max_normalized_to_1 == 0 ~ "Raw / unnormalized",
    TRUE ~ "MIXED - some dyes normalized, some not"
  ))
write_csv(scaling_audit, file.path(OUTDIR, "00_scaling_audit.csv"))

cat("\nScaling audit (native per-dye column maxima):\n")
print(scaling_audit |> select(Config, Treatment, N_dye_columns,
                              N_max_normalized_to_1, Largest_column_max, Native_scale))

scale_is_consistent <- dplyr::n_distinct(scaling_audit$Native_scale) == 1 &&
  scaling_audit$Native_scale[1] != "MIXED - some dyes normalized, some not"

# Magnitude metrics that depend on absolute signal level (total intensity and
# its percent change) are only interpretable if the native scale is consistent
# AND we are not renormalizing. Otherwise they are emitted as NA with a reason,
# rather than silently reported as biology.
INTENSITY_METRICS_VALID <- scale_is_consistent && PREPROCESS == "Raw values"
INTENSITY_INVALID_REASON <- if (INTENSITY_METRICS_VALID) NA_character_ else {
  if (!scale_is_consistent)
    "Source sheets are not on a common scale (see 00_scaling_audit.csv); absolute intensity is not comparable"
  else
    paste0("Spectra renormalized with '", PREPROCESS,
           "'; absolute intensity is removed by normalization")
}

if (!scale_is_consistent) {
  warning(
    "Source workbooks are NOT on a common scale. See 00_scaling_audit.csv. ",
    "Magnitude metrics (Q2a/Q2b/Q2c RMSE, max detector difference) are computed ",
    "on spectra renormalized with '", PREPROCESS, "'; total-intensity metrics are ",
    "suppressed. Q1a peak calls and Q1b cosine are scale-invariant and unaffected.",
    call. = FALSE, immediate. = TRUE
  )
}

write_csv(
  tibble(
    Setting = c("PREPROCESS", "PREPROCESS_SENSITIVITY", "Native scale consistent",
                "Intensity metrics valid", "Intensity suppression reason",
                "Near-tie ratio", "Q1b reference cosine"),
    Value = c(PREPROCESS, as.character(PREPROCESS_SENSITIVITY),
              as.character(scale_is_consistent), as.character(INTENSITY_METRICS_VALID),
              INTENSITY_INVALID_REASON %||% "", as.character(NEAR_TIE_RATIO),
              as.character(Q1B_REFERENCE_COSINE))
  ),
  file.path(OUTDIR, "00_analysis_settings.csv")
)

# ---------------------------------------------------------------------------
# 3c. UV/lost-laser block identified from the DATA, not hard-coded (new in V3)
# ---------------------------------------------------------------------------
# The header already states the intent not to hard-code "320", but the Q1a
# category logic tested Laser_nm == 320 literally. Derive it instead, so the
# same code is correct if a different laser is dropped in a future experiment.
ch_ref <- data5[[common_treatments[1]]]$channels
LOST_LASERS <- sort(unique(ch_ref$Laser_nm[ch_ref$ChannelKey %in% uv_lost_keys]))
LOST_LASERS <- LOST_LASERS[is.finite(LOST_LASERS)]
LOST_LASER_LABEL <- if (length(LOST_LASERS))
  paste0(paste(LOST_LASERS, collapse = "/"), "-nm") else "lost-laser"
cat("Laser block(s) present at 6L but inactive at 5L:", LOST_LASER_LABEL, "\n")

# -------------------------

# ---------------------------------------------------------------------------
# Plot helper — detector-nm scale INSIDE each laser block
# ---------------------------------------------------------------------------
add_laser_detector_axis <- function(p, plot_df, detector_tick_every = 2) {
  axis_meta <- plot_df |>
    dplyr::distinct(ChannelNumber, Laser_nm, Detector_nm) |>
    dplyr::filter(is.finite(ChannelNumber),
                  is.finite(Laser_nm),
                  is.finite(Detector_nm)) |>
    dplyr::arrange(ChannelNumber)

  blocks <- axis_meta |>
    dplyr::group_by(Laser_nm) |>
    dplyr::summarise(
      xmin = min(ChannelNumber),
      xmax = max(ChannelNumber),
      xmid = (xmin + xmax) / 2,
      detector_min = min(Detector_nm),
      detector_max = max(Detector_nm),
      .groups = "drop"
    ) |>
    dplyr::arrange(xmin)

  ticks <- axis_meta |>
    dplyr::group_by(Laser_nm) |>
    dplyr::mutate(
      i = dplyr::row_number(),
      nn = dplyr::n(),
      keep = i == 1L | i == nn |
        ((i - 1L) %% detector_tick_every == 0L)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(keep)

  p +
    ggplot2::geom_vline(
      data = blocks,
      ggplot2::aes(xintercept = xmin),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = .35,
      alpha = .45
    ) +
    ggplot2::scale_x_continuous(
      breaks = ticks$ChannelNumber,
      labels = ticks$Detector_nm,
      minor_breaks = NULL,
      sec.axis = ggplot2::dup_axis(
        breaks = blocks$xmid,
        labels = paste0(blocks$Laser_nm, " nm"),
        name = "Laser blocks"
      )
    ) +
    ggplot2::labs(
      x = "Detector wavelength within each laser block (nm)"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = .5, hjust = 1, size = 6.5
      ),
      axis.text.x.top = ggplot2::element_text(size = 8),
      axis.title.x.top = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(t = 14, r = 8, b = 8, l = 8)
    )
}


# 4. Q1a-A — peak behavior across configurations: native 5L vs native 6L
# -------------------------
# Q1a is intentionally broader than a single 5L-vs-6L test because the collaborator's
# wording "Peak shift?" and project slides do not clearly restrict the question to one
# contrast. Part A holds treatment constant and compares native 5L with native 6L.
# Part B (created below from the same detector-level data) holds configuration constant
# and compares each treatment with Unfixed. Part C compares those treatment-associated
# peak responses between 5L and 6L. Expected/reference peaks remain separate.
# Near-tied maxima are retained throughout so tiny swaps between almost-equal channels
# are not overinterpreted as movement of the whole spectral signature.

top2_peak_info <- function(values, channels) {
  ok <- which(is.finite(values))
  if (!length(ok)) return(tibble(
    Peak_ChannelKey = NA_character_, Peak_ChannelNumber = NA_real_,
    Peak_Laser_nm = NA_real_, Peak_Detector_nm = NA_real_, Peak_Value = NA_real_,
    Second_ChannelKey = NA_character_, Second_ChannelNumber = NA_real_,
    Second_Laser_nm = NA_real_, Second_Detector_nm = NA_real_, Second_Value = NA_real_,
    Peak_margin = NA_real_, Second_to_first_ratio = NA_real_))
  ord <- ok[order(values[ok], decreasing = TRUE, na.last = NA)]
  i1 <- ord[1]; i2 <- if (length(ord) >= 2) ord[2] else NA_integer_
  v1 <- values[i1]; v2 <- if (is.finite(i2)) values[i2] else NA_real_
  tibble(
    Peak_ChannelKey = channels$ChannelKey[i1],
    Peak_ChannelNumber = channels$ChannelNumber[i1],
    Peak_Laser_nm = channels$Laser_nm[i1],
    Peak_Detector_nm = channels$Detector_nm[i1], Peak_Value = v1,
    Second_ChannelKey = if (is.finite(i2)) channels$ChannelKey[i2] else NA_character_,
    Second_ChannelNumber = if (is.finite(i2)) channels$ChannelNumber[i2] else NA_real_,
    Second_Laser_nm = if (is.finite(i2)) channels$Laser_nm[i2] else NA_real_,
    Second_Detector_nm = if (is.finite(i2)) channels$Detector_nm[i2] else NA_real_,
    Second_Value = v2,
    Peak_margin = if (is.finite(v2)) v1 - v2 else NA_real_,
    Second_to_first_ratio = if (is.finite(v2) && is.finite(v1) && v1 != 0) v2 / v1 else NA_real_)
}

q1a <- list()
for (tr in common_treatments) {
  # IMPORTANT: peak is found on each configuration's own native ACTIVE channels.
  # We deliberately do NOT use shared_active here.
  ids <- intersect(colnames(data6[[tr]]$matrix), colnames(data5[[tr]]$matrix))
  ids <- intersect(ids, c(common_fluors, partial_fluors))
  p6 <- restrict_profile(data6[[tr]], active6, ids)
  p5 <- restrict_profile(data5[[tr]], active5, ids)
  for (f in ids) {
    pk6 <- top2_peak_info(p6$matrix[, f], p6$channels)
    pk5 <- top2_peak_info(p5$matrix[, f], p5$channels)
    q1a[[length(q1a) + 1]] <- bind_cols(
      tibble(Treatment = tr, Fluor_ID = f),
      pk6 |> rename_with(~ paste0(.x, "_6L")),
      pk5 |> rename_with(~ paste0(.x, "_5L"))
    ) |>
      mutate(
        Identified_peak_changed = Peak_ChannelKey_5L != Peak_ChannelKey_6L,
        Same_laser = Peak_Laser_nm_5L == Peak_Laser_nm_6L,
        Peak_change_type = case_when(
          is.na(Identified_peak_changed) ~ "Unavailable",
          !Identified_peak_changed ~ "Same identified maximum",
          Identified_peak_changed & Same_laser ~ "Different maximum within same laser",
          Identified_peak_changed & !Same_laser ~ "Different maximum across lasers",
          TRUE ~ "Unavailable"
        ),
        Detector_nm_difference_5L_minus_6L = Peak_Detector_nm_5L - Peak_Detector_nm_6L,
        Abs_detector_nm_difference = abs(Detector_nm_difference_5L_minus_6L),
        # Flag near-tied maxima rather than hiding them. NEAR_TIE_RATIO (default
        # 0.99) means the second channel is at least that fraction of the highest.
        Near_tied_maximum_6L = is.finite(Second_to_first_ratio_6L) & Second_to_first_ratio_6L >= NEAR_TIE_RATIO,
        Near_tied_maximum_5L = is.finite(Second_to_first_ratio_5L) & Second_to_first_ratio_5L >= NEAR_TIE_RATIO,
        Near_tied_either_config = Near_tied_maximum_6L | Near_tied_maximum_5L
      )
  }
}
q1a_peak_shift <- bind_rows(q1a) |>
  arrange(Treatment, desc(Identified_peak_changed), Peak_change_type, Fluor_ID)
write_csv(q1a_peak_shift,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_all_fluorochrome_native_peak_comparisons.csv"))
# Backward-compatible filename used by the QMD reader.
write_csv(q1a_peak_shift,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_all_fluorochrome_peak_shifts.csv"))

q1a_summary <- q1a_peak_shift |>
  group_by(Treatment) |>
  summarise(
    N_fluorochromes = n(),
    N_same_identified_maximum = sum(!Identified_peak_changed, na.rm = TRUE),
    N_changed_identified_maximum = sum(Identified_peak_changed, na.rm = TRUE),
    N_changed_within_same_laser = sum(Peak_change_type == "Different maximum within same laser", na.rm = TRUE),
    N_changed_across_lasers = sum(Peak_change_type == "Different maximum across lasers", na.rm = TRUE),
    N_near_tied_maximum = sum(Near_tied_either_config, na.rm = TRUE),
    Percent_changed = 100 * N_changed_identified_maximum / pmax(N_fluorochromes, 1),
    .groups = "drop")
write_csv(q1a_summary,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_peak_shift_summary_by_treatment.csv"))
print(q1a_summary)

q1a_changed <- q1a_peak_shift |> filter(Identified_peak_changed)
write_csv(q1a_changed,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_changed_identified_maxima.csv"))
q1a_near_ties <- q1a_peak_shift |> filter(Near_tied_either_config)
write_csv(q1a_near_ties,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_near_tied_maxima.csv"))

# Counts figure: directly answers whether the identified maximum changes.
q1a_counts_long <- q1a_peak_shift |>
  count(Treatment, Peak_change_type, name = "N")
p_q1a_counts <- ggplot(q1a_counts_long, aes(Treatment, N, fill = Peak_change_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(N > 0, N, "")),
            position = position_stack(vjust = 0.5), size = 3.3) +
  labs(title = "Q1a — Does the identified spectral maximum change between 5L and 6L?",
       subtitle = "Maximum is identified independently on each configuration's native active spectrum; changes are separated into within-laser and across-laser cases",
       x = NULL, y = "Number of fluorochromes", fill = "Peak comparison") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(OUTDIR, "figures", "Q1a_summary_counts.png"),
       p_q1a_counts, width = 10, height = 6, dpi = 250)

# Detail figure for every changed maximum. This makes missed cases visible.
if (nrow(q1a_changed)) {
  q1a_changed_plot <- q1a_changed |>
    mutate(Pair = paste0(Peak_ChannelKey_6L, "  →  ", Peak_ChannelKey_5L),
           Label = paste0(Fluor_ID, "\n", Pair,
                          ifelse(Near_tied_either_config, "  [near-tied maximum]", "")))
  p_q1a_changed <- ggplot(q1a_changed_plot,
                          aes(x = Treatment, y = reorder(Label, Peak_ChannelNumber_6L),
                              shape = Peak_change_type)) +
    geom_point(size = 3) +
    labs(title = "Q1a — Fluorochromes whose identified maximum differs between 6L and 5L",
         subtitle = "Label shows 6L maximum → 5L maximum; near-tied maxima are explicitly flagged",
         x = NULL, y = "Fluorochrome and identified maximum", shape = "Change type") +
    theme_bw(base_size = 10)
  ggsave(file.path(OUTDIR, "figures", "Q1a_changed_identified_maxima.png"),
         p_q1a_changed, width = 12, height = max(6, 0.28*nrow(q1a_changed_plot)+2), dpi = 250, limitsize = FALSE)
}

# helper: tidy native-axis spectrum for a dye in a given config/treatment
spectrum_df <- function(p_full, fluor, config_label, active_keys) {
  rp <- restrict_profile(p_full, active_keys, fluor)
  tibble(Config = config_label,
         ChannelNumber = rp$channels$ChannelNumber,
         ChannelKey = rp$channels$ChannelKey,
         Detector_nm = rp$channels$Detector_nm,
         Laser_nm = rp$channels$Laser_nm,
         Value = rp$matrix[, 1])
}

# Save an overlay for EVERY changed case, not just one or two examples.
if (nrow(q1a_changed)) {
  for (i in seq_len(nrow(q1a_changed))) {
    f <- q1a_changed$Fluor_ID[i]; tr <- q1a_changed$Treatment[i]
    df <- bind_rows(
      spectrum_df(data6[[tr]], f, "6L", active6),
      spectrum_df(data5[[tr]], f, "5L", active5))
    peaks <- df |> group_by(Config) |> slice_max(Value, n = 1, with_ties = FALSE) |> ungroup()
    # Label the actual physical peak identity on the plot. The x-axis remains
    # native channel order because detector wavelength alone is NOT unique across
    # laser blocks (e.g. 320-583 and 561-583 are different channels but both end
    # in detector 583 nm).
    peaks <- peaks |>
      mutate(
        Peak_label = paste0(
          Config, " MAX\n",
          ChannelKey, "\n",
          "laser ", Laser_nm, " nm | detector ", Detector_nm, " nm"
        )
      )

    # Mark laser-block boundaries and label each active laser along the x-axis.
    laser_blocks <- df |>
      group_by(Config, Laser_nm) |>
      summarise(
        xmin = min(ChannelNumber, na.rm = TRUE),
        xmax = max(ChannelNumber, na.rm = TRUE),
        xmid = (xmin + xmax) / 2,
        .groups = "drop"
      )

    laser_boundaries <- laser_blocks |>
      group_by(Config) |>
      arrange(xmin, .by_group = TRUE) |>
      mutate(boundary = xmin) |>
      filter(row_number() > 1) |>
      ungroup()

    # Use the union of laser blocks for labels; a block present only in 6L
    # (notably 320 nm UV) is still shown and therefore visibly absent from 5L.
    laser_labels <- laser_blocks |>
      group_by(Laser_nm) |>
      summarise(
        xmin = min(xmin), xmax = max(xmax),
        xmid = median(xmid),
        .groups = "drop"
      ) |>
      arrange(xmin)

    pp <- ggplot(df, aes(ChannelNumber, Value, color = Config)) +
      geom_line(linewidth = .75, na.rm = TRUE) +
      geom_point(data = peaks, size = 3) +
      geom_vline(data = peaks, aes(xintercept = ChannelNumber, color = Config),
                 linetype = "dashed", alpha = .55) +
      labs(
        title = paste0("Q1a — ", f, " / ", tr),
        subtitle = paste0(
          "Native spectra. Laser blocks are shown above; detector wavelengths are shown below within each block. Peak dots and dashed lines mark the identified maxima. ",
          "6L ", q1a_changed$Peak_ChannelKey_6L[i],
          " → 5L ", q1a_changed$Peak_ChannelKey_5L[i],
          ifelse(q1a_changed$Near_tied_either_config[i],
                 ". Near-tied maximum present.", ".")
        ),
        y = "Spectral value", color = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top")

    pp <- add_laser_detector_axis(pp, df, detector_tick_every = 2)

    safe <- gsub("[^A-Za-z0-9]+", "_", paste(tr, f, sep = "_"))
    ggsave(file.path(OUTDIR, "figures", paste0("Q1a_changed_", safe, ".png")),
           pp, width = 9, height = 5, dpi = 220)
  }
}

cat(sprintf("Q1a: %d native 5L-vs-6L identified-maximum changes; %d comparisons have a near-tied maximum.\n",
            nrow(q1a_changed), nrow(q1a_near_ties)))

# ---------------------------------------------------------------------------
# 4b. QC MATRICES — detector-level values and easy-to-scan peak summary
# ---------------------------------------------------------------------------
# These exports are intentionally arranged for Excel review. They preserve the
# actual value at every native active detector rather than reducing a spectrum
# immediately to one maximum. This makes near-tied maxima and treatment-related
# peak swaps easy to verify against the source workbook / Jarina's slides.
#
# Output 1 (LONG): one row = fluorochrome x configuration x treatment x channel.
# Output 2 (WIDE): one row = fluorochrome x native channel; eight condition
#                  columns are placed side-by-side (5L first, then 6L).
# Output 3 (PEAK SUMMARY): one row = fluorochrome; peak/second-peak information
#                          for all eight conditions is side-by-side.

QC_DIR <- file.path(OUTDIR, "QC_matrices")
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)

treatment_order_qc <- c("Unfixed", "1% PFA", "4% PFA", "Perm")
condition_order_qc <- c(
  "5L_Unfixed", "5L_1PFA", "5L_4PFA", "5L_Perm",
  "6L_Unfixed", "6L_1PFA", "6L_4PFA", "6L_Perm"
)

qc_condition_label <- function(config, treatment) {
  tr <- dplyr::recode(treatment,
                      "Unfixed" = "Unfixed", "1% PFA" = "1PFA",
                      "4% PFA" = "4PFA", "Perm" = "Perm",
                      .default = treatment)
  paste(config, tr, sep = "_")
}

# Build detector-level native ACTIVE data for every available dye/condition.
qc_long_list <- list()
for (cfg in c("5L", "6L")) {
  dat <- if (cfg == "5L") data5 else data6
  active_keys <- if (cfg == "5L") active5 else active6

  for (tr in intersect(treatment_order_qc, names(dat))) {
    p <- dat[[tr]]
    ids <- setdiff(colnames(p$matrix), stray_ids)
    rp <- restrict_profile(p, active_keys, ids)

    # Matrix -> long table while retaining physical detector identity.
    vals <- as.data.frame(rp$matrix, check.names = FALSE) |>
      mutate(.QC_row = row_number()) |>
      pivot_longer(-.QC_row, names_to = "Fluor_ID", values_to = "Spectral_Value")

    ch <- rp$channels |>
      mutate(.QC_row = row_number()) |>
      select(.QC_row, ChannelNumber, ChannelKey, Laser_nm, Detector_nm,
             Channel_A, Channel_B, Channel_C)

    qc_long_list[[length(qc_long_list) + 1]] <- vals |>
      left_join(ch, by = ".QC_row") |>
      mutate(Config = cfg,
             Treatment = tr,
             Condition = qc_condition_label(cfg, tr)) |>
      select(Fluor_ID, Config, Treatment, Condition,
             ChannelNumber, ChannelKey, Laser_nm, Detector_nm,
             Channel_A, Channel_B, Channel_C, Spectral_Value)
  }
}

qc_channel_values_long <- bind_rows(qc_long_list) |>
  mutate(
    Config = factor(Config, levels = c("5L", "6L")),
    Treatment = factor(Treatment, levels = treatment_order_qc),
    Condition = factor(Condition, levels = condition_order_qc)
  ) |>
  arrange(Fluor_ID, Config, Treatment, ChannelNumber) |>
  mutate(Config = as.character(Config),
         Treatment = as.character(Treatment),
         Condition = as.character(Condition))

write_csv(qc_channel_values_long,
          file.path(QC_DIR, "QC_01_ALL_detector_values_LONG.csv"))

# Wide Excel-friendly detector matrix. ChannelKey is the physical detector key,
# so values from 5L/6L line up where the detector exists in both configurations.
qc_channel_values_wide <- qc_channel_values_long |>
  select(Fluor_ID, ChannelKey, Laser_nm, Detector_nm, Condition, Spectral_Value) |>
  distinct() |>
  pivot_wider(names_from = Condition, values_from = Spectral_Value) |>
  arrange(Fluor_ID, Laser_nm, Detector_nm)

# Force a predictable left-to-right condition order even if a condition is absent.
for (nm in condition_order_qc) {
  if (!nm %in% names(qc_channel_values_wide)) qc_channel_values_wide[[nm]] <- NA_real_
}
qc_channel_values_wide <- qc_channel_values_wide |>
  select(Fluor_ID, ChannelKey, Laser_nm, Detector_nm, all_of(condition_order_qc))

write_csv(qc_channel_values_wide,
          file.path(QC_DIR, "QC_02_detector_values_WIDE_5L_then_6L.csv"))

# Peak/second-peak summary for each condition. In addition to the peak channel,
# retain the runner-up and its closeness to the maximum. This is especially
# useful when a tiny value difference makes the identified maximum swap channels.
qc_peak_long <- qc_channel_values_long |>
  group_by(Fluor_ID, Config, Treatment, Condition) |>
  group_modify(~ {
    z <- .x |> arrange(desc(Spectral_Value), ChannelNumber)
    z <- z[is.finite(z$Spectral_Value), , drop = FALSE]
    if (!nrow(z)) return(tibble(
      Peak_Channel = NA_character_, Peak_Laser_nm = NA_real_, Peak_Detector_nm = NA_real_,
      Peak_Value = NA_real_, Second_Channel = NA_character_, Second_Laser_nm = NA_real_,
      Second_Detector_nm = NA_real_, Second_Value = NA_real_, Peak_minus_second = NA_real_,
      Second_to_first_percent = NA_real_, Near_tied_99pct = NA
    ))
    a <- z[1, ]; b <- if (nrow(z) >= 2) z[2, ] else z[1, ]
    second_ok <- nrow(z) >= 2
    tibble(
      Peak_Channel = a$ChannelKey,
      Peak_Laser_nm = a$Laser_nm,
      Peak_Detector_nm = a$Detector_nm,
      Peak_Value = a$Spectral_Value,
      Second_Channel = if (second_ok) b$ChannelKey else NA_character_,
      Second_Laser_nm = if (second_ok) b$Laser_nm else NA_real_,
      Second_Detector_nm = if (second_ok) b$Detector_nm else NA_real_,
      Second_Value = if (second_ok) b$Spectral_Value else NA_real_,
      Peak_minus_second = if (second_ok) a$Spectral_Value - b$Spectral_Value else NA_real_,
      Second_to_first_percent = if (second_ok && is.finite(a$Spectral_Value) && a$Spectral_Value != 0)
                                  100 * b$Spectral_Value / a$Spectral_Value else NA_real_,
      Near_tied_99pct = if (second_ok && is.finite(a$Spectral_Value) && a$Spectral_Value != 0)
                          (b$Spectral_Value / a$Spectral_Value) >= NEAR_TIE_RATIO else NA
    )
  }) |>
  ungroup() |>
  mutate(Condition = factor(Condition, levels = condition_order_qc)) |>
  arrange(Fluor_ID, Condition)

write_csv(qc_peak_long |> mutate(Condition = as.character(Condition)),
          file.path(QC_DIR, "QC_03_peak_and_second_peak_LONG.csv"))

# Wide peak summary: one fluorochrome per row. Conditions are grouped in a
# consistent order so it can be scanned horizontally in Excel.
qc_peak_wide <- qc_peak_long |>
  select(Fluor_ID, Condition, Peak_Channel, Peak_Laser_nm, Peak_Detector_nm,
         Peak_Value, Second_Channel, Second_Laser_nm, Second_Detector_nm,
         Second_Value, Peak_minus_second, Second_to_first_percent, Near_tied_99pct) |>
  pivot_wider(
    names_from = Condition,
    values_from = c(Peak_Channel, Peak_Laser_nm, Peak_Detector_nm, Peak_Value,
                    Second_Channel, Second_Laser_nm, Second_Detector_nm, Second_Value,
                    Peak_minus_second, Second_to_first_percent, Near_tied_99pct),
    names_glue = "{Condition}__{.value}"
  ) |>
  arrange(Fluor_ID)

# Reorder wide columns by CONDITION first (5L Unfixed -> ... -> 6L Perm), then
# by peak fields. This is much easier to inspect than pivot_wider's default
# metric-first ordering.
peak_fields_qc <- c("Peak_Channel", "Peak_Laser_nm", "Peak_Detector_nm", "Peak_Value",
                    "Second_Channel", "Second_Laser_nm", "Second_Detector_nm", "Second_Value",
                    "Peak_minus_second", "Second_to_first_percent", "Near_tied_99pct")
ordered_peak_cols <- unlist(lapply(condition_order_qc, function(cc)
  paste0(cc, "__", peak_fields_qc)))
ordered_peak_cols <- ordered_peak_cols[ordered_peak_cols %in% names(qc_peak_wide)]
qc_peak_wide <- qc_peak_wide |> select(Fluor_ID, all_of(ordered_peak_cols))

write_csv(qc_peak_wide,
          file.path(QC_DIR, "QC_04_peak_summary_WIDE_one_row_per_fluor.csv"))

# Compact change flags: lets us immediately filter Excel to cases where fixation
# changes the maximum within 5L or 6L, and to 5L-vs-6L differences per treatment.
qc_peak_flags <- qc_peak_long |>
  select(Fluor_ID, Config, Treatment, Peak_Channel, Peak_Laser_nm,
         Peak_Detector_nm, Peak_Value, Second_Channel, Second_Value,
         Second_to_first_percent, Near_tied_99pct) |>
  pivot_wider(names_from = Config,
              values_from = c(Peak_Channel, Peak_Laser_nm, Peak_Detector_nm,
                              Peak_Value, Second_Channel, Second_Value,
                              Second_to_first_percent, Near_tied_99pct),
              names_sep = "_") |>
  mutate(
    # V3: keep NA explicitly NA (an unavailable comparison is not "no difference"),
    # but make the comparison itself total so downstream filters behave predictably.
    Peak_comparison_available = !is.na(Peak_Channel_5L) & !is.na(Peak_Channel_6L),
    Peak_differs_5L_vs_6L = if_else(Peak_comparison_available,
                                    Peak_Channel_5L != Peak_Channel_6L, NA),
    Same_laser_5L_vs_6L = if_else(Peak_comparison_available,
                                  Peak_Laser_nm_5L == Peak_Laser_nm_6L, NA),
    Detector_nm_diff_5L_minus_6L = Peak_Detector_nm_5L - Peak_Detector_nm_6L
  ) |>
  arrange(factor(Treatment, levels = treatment_order_qc),
          desc(Peak_differs_5L_vs_6L), Fluor_ID)

write_csv(qc_peak_flags,
          file.path(QC_DIR, "QC_05_peak_flags_5L_vs_6L_by_treatment.csv"))

# Fixation-specific peak changes relative to Unfixed, separately within 5L/6L.
qc_fixation_peak_flags <- qc_peak_long |>
  select(Fluor_ID, Config, Treatment, Peak_Channel, Peak_Laser_nm,
         Peak_Detector_nm, Peak_Value, Second_Channel, Second_Value,
         Second_to_first_percent, Near_tied_99pct) |>
  group_by(Fluor_ID, Config) |>
  mutate(
    Unfixed_Peak_Channel = Peak_Channel[Treatment == "Unfixed"][1],
    Unfixed_Peak_Laser_nm = Peak_Laser_nm[Treatment == "Unfixed"][1],
    Unfixed_Peak_Detector_nm = Peak_Detector_nm[Treatment == "Unfixed"][1]
  ) |>
  ungroup() |>
  filter(Treatment != "Unfixed") |>
  mutate(
    Treatment_comparison_available = !is.na(Peak_Channel) & !is.na(Unfixed_Peak_Channel),
    Peak_changed_vs_Unfixed = if_else(Treatment_comparison_available,
                                      Peak_Channel != Unfixed_Peak_Channel, NA),
    Peak_change_same_laser = Peak_changed_vs_Unfixed %in% TRUE &
                             !is.na(Peak_Laser_nm) & !is.na(Unfixed_Peak_Laser_nm) &
                             Peak_Laser_nm == Unfixed_Peak_Laser_nm,
    Detector_nm_shift_vs_Unfixed = Peak_Detector_nm - Unfixed_Peak_Detector_nm
  ) |>
  arrange(Config, factor(Treatment, levels = treatment_order_qc),
          desc(Peak_changed_vs_Unfixed), Fluor_ID)

write_csv(qc_fixation_peak_flags,
          file.path(QC_DIR, "QC_06_fixation_peak_changes_vs_Unfixed.csv"))

# ---------------------------------------------------------------------------
# Q1a-B / Q1a-C — treatment-associated peak behavior within each configuration
# ---------------------------------------------------------------------------
# These are promoted from QC to formal Q1a outputs so "Peak shift?" is answered in
# both directions: across configurations AND after treatment within a configuration.
q1a_treatment_peak_changes <- qc_fixation_peak_flags |>
  mutate(
    Peak_change_type_vs_Unfixed = case_when(
      !Peak_changed_vs_Unfixed ~ "Same identified maximum",
      Peak_changed_vs_Unfixed & Peak_change_same_laser ~ "Different maximum within same laser",
      Peak_changed_vs_Unfixed & !Peak_change_same_laser ~ "Different maximum across lasers",
      TRUE ~ "Unavailable"
    )
  ) |>
  arrange(Fluor_ID, Config, factor(Treatment, levels = treatment_order_qc))

write_csv(q1a_treatment_peak_changes,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_treatment_vs_Unfixed_peak_changes_within_config.csv"))

q1a_treatment_peak_summary <- q1a_treatment_peak_changes |>
  group_by(Config, Treatment) |>
  summarise(
    N_fluorochromes = n(),
    N_same_identified_maximum = sum(!Peak_changed_vs_Unfixed, na.rm = TRUE),
    N_changed_identified_maximum = sum(Peak_changed_vs_Unfixed, na.rm = TRUE),
    N_changed_within_same_laser = sum(Peak_change_same_laser, na.rm = TRUE),
    N_changed_across_lasers = sum(Peak_changed_vs_Unfixed & !Peak_change_same_laser, na.rm = TRUE),
    N_near_tied_treated = sum(Near_tied_99pct, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(q1a_treatment_peak_summary,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_treatment_vs_Unfixed_peak_summary.csv"))

# One row per fluorochrome/treatment: did treatment change the peak in 5L, 6L, both, or neither?
q1a_peak_response_5L_vs_6L <- q1a_treatment_peak_changes |>
  select(Fluor_ID, Config, Treatment, Peak_Channel, Unfixed_Peak_Channel,
         Peak_changed_vs_Unfixed, Peak_change_same_laser, Detector_nm_shift_vs_Unfixed,
         Near_tied_99pct) |>
  pivot_wider(
    names_from = Config,
    values_from = c(Peak_Channel, Unfixed_Peak_Channel, Peak_changed_vs_Unfixed,
                    Peak_change_same_laser, Detector_nm_shift_vs_Unfixed, Near_tied_99pct),
    names_sep = "_"
  ) |>
  mutate(
    Treatment_peak_response = case_when(
      Peak_changed_vs_Unfixed_5L & Peak_changed_vs_Unfixed_6L ~ "Peak changed in both 5L and 6L",
      Peak_changed_vs_Unfixed_5L & !Peak_changed_vs_Unfixed_6L ~ "Peak changed in 5L only",
      !Peak_changed_vs_Unfixed_5L & Peak_changed_vs_Unfixed_6L ~ "Peak changed in 6L only",
      !Peak_changed_vs_Unfixed_5L & !Peak_changed_vs_Unfixed_6L ~ "No peak change in either configuration",
      TRUE ~ "Unavailable"
    ),
    Same_treatment_peak_response_5L_6L = Peak_changed_vs_Unfixed_5L == Peak_changed_vs_Unfixed_6L
  ) |>
  arrange(factor(Treatment, levels = treatment_order_qc),
          desc(!Same_treatment_peak_response_5L_6L), Fluor_ID)

write_csv(q1a_peak_response_5L_vs_6L,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_treatment_peak_response_5L_vs_6L.csv"))

# Excel-friendly summary of the collaborator-slide style flags: ANY treatment-associated
# same-laser peak change within 5L and within 6L for each fluorochrome.
q1a_any_treatment_peak_flag <- q1a_treatment_peak_changes |>
  group_by(Fluor_ID, Config) |>
  summarise(
    Any_peak_change_after_treatment = any(Peak_changed_vs_Unfixed, na.rm = TRUE),
    Any_same_laser_peak_change_after_treatment = any(Peak_change_same_laser, na.rm = TRUE),
    Treatments_with_same_laser_peak_change = paste(Treatment[Peak_change_same_laser], collapse = "; "),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = Config,
              values_from = c(Any_peak_change_after_treatment,
                              Any_same_laser_peak_change_after_treatment,
                              Treatments_with_same_laser_peak_change),
              names_sep = "_") |>
  arrange(desc(Any_same_laser_peak_change_after_treatment_5L),
          desc(Any_same_laser_peak_change_after_treatment_6L), Fluor_ID)

write_csv(q1a_any_treatment_peak_flag,
          file.path(OUTDIR, "Q1a_peak_shift", "Q1a_ANY_treatment_peak_change_flags_5L_6L.csv"))

cat("Q1a comprehensive peak outputs written: across-config, treatment-vs-Unfixed, and 5L-vs-6L treatment response.\n")

cat("\nQC matrices written to:", normalizePath(QC_DIR), "\n")
cat("  QC_01: all detector values, long format\n")
cat("  QC_02: detector values, wide 5L -> 6L format\n")
cat("  QC_03: peak + second peak, long format\n")
cat("  QC_04: peak summary, one fluorochrome per row\n")
cat("  QC_05: 5L-vs-6L peak flags by treatment\n")
cat("  QC_06: fixation peak changes vs Unfixed within each configuration\n\n")


# ---------------------------------------------------------------------------
# 4c. Q1a COMBINED PEAK-OUTCOME TABLE + FOUR-SPECTRUM QC PLOTS
# ---------------------------------------------------------------------------
# One row = fluorochrome x fixation treatment (1% PFA / 4% PFA / Perm).
# This joins BOTH directions of the broad Q1a "Peak shift?" question:
#
#   horizontal: Unfixed -> Treatment within 5L and within 6L
#   vertical:   5L -> 6L for Unfixed and for the Treatment
#
# Therefore each row/plot can show whether the maximum is:
#   - stable everywhere,
#   - configuration-dependent only,
#   - changed by treatment in both configs (same or different response),
#   - changed by treatment in 5L only,
#   - changed by treatment in 6L only,
#   - or affected by both configuration and treatment.
#
# Near-tied maxima are retained as QC flags and do NOT by themselves redefine
# the categorical peak call.

Q1A_COMBINED_DIR <- file.path(OUTDIR, "Q1a_peak_shift", "combined_peak_outcomes")
Q1A_COMBINED_FIG_DIR <- file.path(OUTDIR, "figures", "Q1a_combined_peak_outcomes")
dir.create(Q1A_COMBINED_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(Q1A_COMBINED_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

treated_q1a <- setdiff(treatment_order_qc, "Unfixed")

# Separate formal sheet 1: across 5L vs 6L, including Unfixed and every treatment.
q1a_across_config_sheet <- qc_peak_flags |>
  mutate(
    Comparison = "5L vs 6L within condition",
    Treatment = as.character(Treatment)
  ) |>
  arrange(factor(Treatment, levels = treatment_order_qc),
          desc(Peak_differs_5L_vs_6L), Fluor_ID)

write_csv(
  q1a_across_config_sheet,
  file.path(Q1A_COMBINED_DIR, "Q1a_01_ACROSS_5L_vs_6L_by_condition.csv")
)

# Separate formal sheet 2: each treatment vs Unfixed, separately within 5L and 6L.
q1a_within_config_sheet <- q1a_treatment_peak_changes |>
  mutate(Comparison = "Treatment vs Unfixed within configuration") |>
  arrange(factor(Treatment, levels = treated_q1a), Config,
          desc(Peak_changed_vs_Unfixed), Fluor_ID)

write_csv(
  q1a_within_config_sheet,
  file.path(Q1A_COMBINED_DIR, "Q1a_02_TREATMENT_vs_Unfixed_within_5L_6L.csv")
)

# Build one-row-per-fluor/treatment treatment table.
tx_wide <- q1a_treatment_peak_changes |>
  filter(Treatment %in% treated_q1a) |>
  select(
    Fluor_ID, Treatment, Config,
    Unfixed_Peak_Channel, Unfixed_Peak_Laser_nm, Unfixed_Peak_Detector_nm,
    Peak_Channel, Peak_Laser_nm, Peak_Detector_nm,
    Peak_changed_vs_Unfixed, Peak_change_same_laser,
    Detector_nm_shift_vs_Unfixed, Near_tied_99pct,
    Second_Channel, Second_Value, Second_to_first_percent
  ) |>
  pivot_wider(
    names_from = Config,
    values_from = c(
      Unfixed_Peak_Channel, Unfixed_Peak_Laser_nm, Unfixed_Peak_Detector_nm,
      Peak_Channel, Peak_Laser_nm, Peak_Detector_nm,
      Peak_changed_vs_Unfixed, Peak_change_same_laser,
      Detector_nm_shift_vs_Unfixed, Near_tied_99pct,
      Second_Channel, Second_Value, Second_to_first_percent
    ),
    names_sep = "_"
  )

# Across-configuration flags for the same treated condition.
treated_across <- qc_peak_flags |>
  filter(Treatment %in% treated_q1a) |>
  transmute(
    Fluor_ID, Treatment,
    Treated_5L_vs_6L_changed = Peak_differs_5L_vs_6L,
    Treated_same_laser_5L_vs_6L = Same_laser_5L_vs_6L,
    Treated_detector_nm_diff_5L_minus_6L = Detector_nm_diff_5L_minus_6L,
    Treated_5L_Peak = Peak_Channel_5L,
    Treated_6L_Peak = Peak_Channel_6L
  )

# Across-configuration flag for Unfixed; repeat it beside every treatment.
unfixed_across <- qc_peak_flags |>
  filter(Treatment == "Unfixed") |>
  transmute(
    Fluor_ID,
    Unfixed_5L_vs_6L_changed = Peak_differs_5L_vs_6L,
    Unfixed_same_laser_5L_vs_6L = Same_laser_5L_vs_6L,
    Unfixed_detector_nm_diff_5L_minus_6L = Detector_nm_diff_5L_minus_6L,
    Unfixed_5L_Peak_across = Peak_Channel_5L,
    Unfixed_6L_Peak_across = Peak_Channel_6L
  )

q1a_combined <- tx_wide |>
  left_join(treated_across, by = c("Fluor_ID", "Treatment")) |>
  left_join(unfixed_across, by = "Fluor_ID") |>
  mutate(
    # Did treatment alter the maximum in each configuration?
    Treatment_change_5L = Peak_changed_vs_Unfixed_5L,
    Treatment_change_6L = Peak_changed_vs_Unfixed_6L,

    # If both configurations changed after treatment, did they land on the same
    # treated maximum? This separates "same response" from "different response".
    Both_configs_changed = Treatment_change_5L %in% TRUE & Treatment_change_6L %in% TRUE,
    Same_treated_peak_after_change =
      Both_configs_changed &
      !is.na(Peak_Channel_5L) & !is.na(Peak_Channel_6L) &
      Peak_Channel_5L == Peak_Channel_6L,

    Any_configuration_difference =
      Unfixed_5L_vs_6L_changed %in% TRUE | Treated_5L_vs_6L_changed %in% TRUE,

    # Describe WHERE each 5L-vs-6L comparison occurs independently of the
    # mutually exclusive reporting category.  "Comparable area" can therefore
    # coexist with a treatment-effect category (e.g. AF555 / Perm).
    # V3: the "lost laser" is derived from the channel audit (LOST_LASERS)
    # rather than hard-coded as 320, so the same logic holds if a different
    # laser is dropped in another experiment or instrument.
    Unfixed_comparison_area = case_when(
      is.na(Unfixed_Peak_Channel_5L) | is.na(Unfixed_Peak_Channel_6L) ~
        "Unavailable",
      Unfixed_Peak_Channel_5L == Unfixed_Peak_Channel_6L ~
        "Same peak",
      Unfixed_Peak_Laser_nm_6L %in% LOST_LASERS ~
        paste0(LOST_LASER_LABEL, " dUV area absent from 5L"),
      TRUE ~ "Comparable area"
    ),
    Treated_comparison_area = case_when(
      is.na(Peak_Channel_5L) | is.na(Peak_Channel_6L) ~
        "Unavailable",
      Peak_Channel_5L == Peak_Channel_6L ~
        "Same peak",
      Peak_Laser_nm_6L %in% LOST_LASERS ~
        paste0(LOST_LASER_LABEL, " dUV area absent from 5L"),
      TRUE ~ "Comparable area"
    ),

    Any_treatment_change =
      Treatment_change_5L %in% TRUE | Treatment_change_6L %in% TRUE,

    Any_peak_change =
      Any_configuration_difference | Any_treatment_change,

    Near_tied_any =
      Near_tied_99pct_5L %in% TRUE | Near_tied_99pct_6L %in% TRUE,

    # V3: a near-tied maximum does not change the category, but a reader should
    # not have to open a separate QC file to learn that a "peak shift" rests on
    # a <1% difference between two channels. Promote it to a first-class field.
    Peak_call_confidence = if_else(
      Near_tied_any,
      "Provisional - runner-up within near-tie margin",
      "Firm - clear single maximum"
    ),

    # Six-group reporting category used by the QMD.
    # Configuration-only categories require NO treatment-associated peak change.
    # "dUV Block Difference" is reserved for a 6L maximum in the 320-nm dUV block that
    # does not exist in 5L; other configuration-only differences are assigned to
    # "Comparable Area" when they occur in laser regions available to both.
    # Treatment categories take precedence whenever treatment changes the maximum.
    # The separate Unfixed and Treated 5L-vs-6L flags are retained for audit.
    Reporting_category = case_when(
      !Any_treatment_change & !Any_configuration_difference ~
        "Stable",
      !Any_treatment_change & Any_configuration_difference &
        (Unfixed_Peak_Laser_nm_6L %in% LOST_LASERS | Peak_Laser_nm_6L %in% LOST_LASERS) ~
        "Configuration Only – dUV Block Difference",
      !Any_treatment_change & Any_configuration_difference ~
        "Configuration Only – Comparable Area",
      Treatment_change_5L %in% TRUE & !(Treatment_change_6L %in% TRUE) ~
        "Treatment Effect – 5L Only",
      !(Treatment_change_5L %in% TRUE) & Treatment_change_6L %in% TRUE ~
        "Treatment Effect – 6L Only",
      Treatment_change_5L %in% TRUE & Treatment_change_6L %in% TRUE ~
        "Treatment Effect – Both Configurations",
      TRUE ~ "Unavailable / incomplete comparison"
    ),

    # Detailed mutually exclusive outcome category.
    Peak_outcome_category = case_when(
      !Any_configuration_difference & !Any_treatment_change ~
        "Stable: no configuration or treatment peak change",

      !Any_configuration_difference &
        Both_configs_changed & Same_treated_peak_after_change ~
        "Treatment both configs - same response",

      !Any_configuration_difference &
        Both_configs_changed & !Same_treated_peak_after_change ~
        "Treatment both configs - different response",

      !Any_configuration_difference &
        Treatment_change_5L %in% TRUE & !(Treatment_change_6L %in% TRUE) ~
        "Treatment 5L only",

      !Any_configuration_difference &
        !(Treatment_change_5L %in% TRUE) & Treatment_change_6L %in% TRUE ~
        "Treatment 6L only",

      Any_configuration_difference & !Any_treatment_change ~
        "Configuration only",

      Any_configuration_difference &
        Both_configs_changed & Same_treated_peak_after_change ~
        "Configuration + treatment both configs - same treated peak",

      Any_configuration_difference &
        Both_configs_changed & !Same_treated_peak_after_change ~
        "Configuration + treatment both configs - different treated peaks",

      Any_configuration_difference &
        Treatment_change_5L %in% TRUE & !(Treatment_change_6L %in% TRUE) ~
        "Configuration + treatment 5L only",

      Any_configuration_difference &
        !(Treatment_change_5L %in% TRUE) & Treatment_change_6L %in% TRUE ~
        "Configuration + treatment 6L only",

      TRUE ~ "Unavailable / incomplete comparison"
    )
  ) |>
  select(
    Fluor_ID, Treatment, Reporting_category, Peak_outcome_category,
    Any_peak_change, Any_configuration_difference, Any_treatment_change,
    Treatment_change_5L, Treatment_change_6L,
    Unfixed_5L_vs_6L_changed, Treated_5L_vs_6L_changed,
    Unfixed_comparison_area, Treated_comparison_area,
    Same_treated_peak_after_change, Near_tied_any, Peak_call_confidence,

    Unfixed_Peak_Channel_5L, Peak_Channel_5L,
    Unfixed_Peak_Channel_6L, Peak_Channel_6L,

    Unfixed_Peak_Laser_nm_5L, Peak_Laser_nm_5L,
    Unfixed_Peak_Detector_nm_5L, Peak_Detector_nm_5L,
    Detector_nm_shift_vs_Unfixed_5L, Peak_change_same_laser_5L,
    Near_tied_99pct_5L, Second_Channel_5L, Second_to_first_percent_5L,

    Unfixed_Peak_Laser_nm_6L, Peak_Laser_nm_6L,
    Unfixed_Peak_Detector_nm_6L, Peak_Detector_nm_6L,
    Detector_nm_shift_vs_Unfixed_6L, Peak_change_same_laser_6L,
    Near_tied_99pct_6L, Second_Channel_6L, Second_to_first_percent_6L,

    Unfixed_detector_nm_diff_5L_minus_6L,
    Treated_detector_nm_diff_5L_minus_6L,
    Unfixed_same_laser_5L_vs_6L, Treated_same_laser_5L_vs_6L
  ) |>
  arrange(
    factor(Treatment, levels = treated_q1a),
    desc(Any_peak_change),
    Peak_outcome_category,
    Fluor_ID
  )

write_csv(
  q1a_combined,
  file.path(Q1A_COMBINED_DIR, "Q1a_03_COMBINED_peak_outcomes_one_row_per_fluor_treatment.csv")
)

# A compact combined sheet for quick Excel scanning.
q1a_combined_compact <- q1a_combined |>
  transmute(
    Fluor_ID, Treatment, Reporting_category, Peak_outcome_category,
    `5L Unfixed peak` = Unfixed_Peak_Channel_5L,
    `5L Treated peak` = Peak_Channel_5L,
    `5L Treatment changed` = Treatment_change_5L,
    `6L Unfixed peak` = Unfixed_Peak_Channel_6L,
    `6L Treated peak` = Peak_Channel_6L,
    `6L Treatment changed` = Treatment_change_6L,
    `Unfixed 5L vs 6L changed` = Unfixed_5L_vs_6L_changed,
    `Treated 5L vs 6L changed` = Treated_5L_vs_6L_changed,
    `Same treated peak after both changed` = Same_treated_peak_after_change,
    `Near tied any treated maximum` = Near_tied_any
  )

write_csv(
  q1a_combined_compact,
  file.path(Q1A_COMBINED_DIR, "Q1a_04_COMBINED_COMPACT_for_Excel.csv")
)

# Main Q1a reporting table: one row per fluorochrome x treatment, assigned to
# exactly one of the six reader-facing categories. Detailed analysis flags remain
# in q1a_combined and are included here where they help audit the classification.
q1a_reporting_table <- q1a_combined |>
  mutate(
    Plot_file = if_else(
      Any_peak_change %in% TRUE,
      paste0(
        "Q1a_combined_",
        gsub("[^A-Za-z0-9]+", "_", paste(Treatment, Fluor_ID, sep = "_")),
        ".png"
      ),
      NA_character_
    ),
    Plot_relative_path = if_else(
      Any_peak_change %in% TRUE,
      file.path("figures", "Q1a_combined_peak_outcomes", Plot_file),
      NA_character_
    )
  ) |>
  transmute(
    Fluorochrome = Fluor_ID, Treatment, Reporting_category,
    `5L Unfixed peak` = Unfixed_Peak_Channel_5L,
    `5L Treated peak` = Peak_Channel_5L,
    `6L Unfixed peak` = Unfixed_Peak_Channel_6L,
    `6L Treated peak` = Peak_Channel_6L,
    # Keep the two configuration comparisons separate. A single combined
    # "Configuration difference" flag was ambiguous because it could be TRUE
    # only after treatment even when the Unfixed 5L and 6L peaks were identical.
    `Unfixed 5L vs 6L difference` = Unfixed_5L_vs_6L_changed,
    `Unfixed comparison area` = Unfixed_comparison_area,
    `Treated 5L vs 6L difference` = Treated_5L_vs_6L_changed,
    `Treated comparison area` = Treated_comparison_area,
    `Treatment change 5L` = Treatment_change_5L,
    `Treatment change 6L` = Treatment_change_6L,
    `Near tie` = Near_tied_any,
    `Peak call confidence` = Peak_call_confidence,
    Plot_relative_path
  ) |>
  arrange(
    factor(Reporting_category, levels = c(
      "Configuration Only – dUV Block Difference",
      "Configuration Only – Comparable Area",
      "Treatment Effect – 5L Only",
      "Treatment Effect – 6L Only",
      "Treatment Effect – Both Configurations",
      "Stable",
      "Unavailable / incomplete comparison"
    )),
    factor(Treatment, levels = treated_q1a), Fluorochrome
  )

write_csv(
  q1a_reporting_table,
  file.path(Q1A_COMBINED_DIR, "Q1a_00_MAIN_six_category_reporting_table.csv")
)

# Six-category summary used by the QMD.
q1a_reporting_category_summary <- q1a_combined |>
  count(Treatment, Reporting_category, name = "N") |>
  group_by(Treatment) |>
  mutate(Percent = 100 * N / sum(N)) |>
  ungroup() |>
  arrange(factor(Treatment, levels = treated_q1a), desc(N))

write_csv(
  q1a_reporting_category_summary,
  file.path(Q1A_COMBINED_DIR, "Q1a_00_six_category_counts.csv")
)

# Category summary: verifies that every fluor x treatment lands in one detailed bucket.
q1a_category_summary <- q1a_combined |>
  count(Treatment, Peak_outcome_category, name = "N") |>
  group_by(Treatment) |>
  mutate(Percent = 100 * N / sum(N)) |>
  ungroup() |>
  arrange(factor(Treatment, levels = treated_q1a), desc(N))

write_csv(
  q1a_category_summary,
  file.path(Q1A_COMBINED_DIR, "Q1a_05_peak_outcome_category_counts.csv")
)

# Also make one combined CSV containing only non-stable cases.
q1a_nonstable <- q1a_combined |>
  filter(Any_peak_change %in% TRUE)

write_csv(
  q1a_nonstable,
  file.path(Q1A_COMBINED_DIR, "Q1a_06_COMBINED_NONSTABLE_cases_only.csv")
)

# Summary figure of ALL outcome categories.
p_q1a_categories <- ggplot(
  q1a_reporting_category_summary,
  aes(x = Treatment, y = N, fill = Reporting_category)
) +
  geom_col(width = 0.72) +
  labs(
    title = "Q1a — Combined peak outcomes across configuration and treatment",
    subtitle = "Each fluorochrome × treatment is assigned to one of the six reporting categories",
    x = NULL, y = "Number of fluorochromes", fill = "Reporting category"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(OUTDIR, "figures", "Q1a_COMBINED_peak_outcome_categories.png"),
  p_q1a_categories, width = 13, height = 8, dpi = 300
)

# Helper for the four-spectrum QC plot.
q1a_four_spectrum_df <- function(fluor, treatment) {
  bind_rows(
    spectrum_df(data5[["Unfixed"]], fluor, "5L Unfixed", active5),
    spectrum_df(data5[[treatment]], fluor, paste0("5L ", treatment), active5),
    spectrum_df(data6[["Unfixed"]], fluor, "6L Unfixed", active6),
    spectrum_df(data6[[treatment]], fluor, paste0("6L ", treatment), active6)
  )
}

# Channel metadata are instrument metadata and do not depend on treatment.
# Define this BEFORE the combined plotting loop so a clean-session run is valid.
meta_treatment <- common_treatments[1]

# One FOUR-spectrum plot for every non-stable fluorochrome x treatment.
# The plot simultaneously shows:
#   5L Unfixed -> 5L Treatment
#   6L Unfixed -> 6L Treatment
# and allows the 5L/6L peaks to be compared vertically.
if (nrow(q1a_nonstable)) {
  for (i in seq_len(nrow(q1a_nonstable))) {
    rr <- q1a_nonstable[i, ]
    f <- rr$Fluor_ID
    tr <- rr$Treatment

    # Skip incomplete dye/treatment combinations cleanly rather than stopping the run.
    needed <- c(
      f %in% colnames(data5[["Unfixed"]]$matrix),
      f %in% colnames(data6[["Unfixed"]]$matrix),
      f %in% colnames(data5[[tr]]$matrix),
      f %in% colnames(data6[[tr]]$matrix)
    )
    if (!all(needed)) next

    df4 <- q1a_four_spectrum_df(f, tr) |>
      mutate(
        Config = ifelse(grepl("^5L", Config), "5L", "6L"),
        State = ifelse(grepl("Unfixed$", Config), "Unfixed", tr)
      )

    # spectrum_df names its condition in Config; preserve that before making facets.
    df4 <- q1a_four_spectrum_df(f, tr) |>
      rename(Condition = Config) |>
      mutate(
        Config = ifelse(grepl("^5L", Condition), "5L", "6L"),
        State = ifelse(grepl("Unfixed$", Condition), "Unfixed", tr),
        State = factor(State, levels = c("Unfixed", tr))
      )

    peaks4 <- df4 |>
      group_by(Config, State, Condition) |>
      slice_max(Value, n = 1, with_ties = FALSE) |>
      ungroup()

    subtitle_txt <- paste0(
      "Category: ", rr$Reporting_category,
      " | 5L treatment change: ", ifelse(rr$Treatment_change_5L, "YES", "NO"),
      " | 6L treatment change: ", ifelse(rr$Treatment_change_6L, "YES", "NO"),
      " | Unfixed 5L↔6L: ", ifelse(rr$Unfixed_5L_vs_6L_changed, "DIFFERENT", "SAME"),
      " | Treated 5L↔6L: ", ifelse(rr$Treated_5L_vs_6L_changed, "DIFFERENT", "SAME"),
      ifelse(rr$Near_tied_any, " | near-tied treated maximum present", "")
    )

    # Make configuration effects visually explicit.
    # Both facets use the SAME native x-axis. The inactive 320-nm dUV block is
    # therefore visible as an intentionally empty/shaded region in 5L rather
    # than disappearing when free_x rescales the panel.
    channel_meta6 <- data6[[meta_treatment]]$channels
    channel_meta5 <- data5[[meta_treatment]]$channels

    axis_min4 <- min(c(channel_meta5$ChannelNumber, channel_meta6$ChannelNumber), na.rm = TRUE)
    axis_max4 <- max(c(channel_meta5$ChannelNumber, channel_meta6$ChannelNumber), na.rm = TRUE)

    # Use 6L metadata to define the physical laser-block spans because 6L has
    # all six active lasers. These coordinates are then repeated in both facets.
    laser_blocks_base4 <- channel_meta6 |>
      filter(Laser_nm %in% c(320, 355, 405, 488, 561, 637)) |>
      group_by(Laser_nm) |>
      summarise(
        xmin = min(ChannelNumber, na.rm = TRUE) - 0.5,
        xmax = max(ChannelNumber, na.rm = TRUE) + 0.5,
        xmid = (xmin + xmax) / 2,
        .groups = "drop"
      ) |>
      mutate(
        Laser_label = case_when(
          Laser_nm == 320 ~ "320 dUV",
          Laser_nm == 355 ~ "355 UV",
          TRUE ~ paste0(Laser_nm, " nm")
        )
      )

    laser_blocks4 <- tidyr::crossing(
      Config = c("5L", "6L"),
      laser_blocks_base4
    )

    uv320_band4 <- laser_blocks_base4 |>
      filter(Laser_nm == 320) |>
      transmute(
        Config = "5L", xmin, xmax, xmid,
        ymin = -Inf, ymax = Inf,
        Missing_label = "320 dUV\\nABSENT IN 5L"
      )

    # Short peak summary in the subtitle. This explicitly separates excitation
    # laser identity from detector wavelength, which is essential for cases such
    # as AF350 (355-432 in 5L versus 320-432 in 6L: same detector, different laser).
    pk5_u <- peaks4 |> filter(Config == "5L", State == "Unfixed") |> slice(1)
    pk5_t <- peaks4 |> filter(Config == "5L", State == tr) |> slice(1)
    pk6_u <- peaks4 |> filter(Config == "6L", State == "Unfixed") |> slice(1)
    pk6_t <- peaks4 |> filter(Config == "6L", State == tr) |> slice(1)

    peak_line4 <- paste0(
      "5L: ", pk5_u$ChannelKey, " → ", pk5_t$ChannelKey,
      "   |   6L: ", pk6_u$ChannelKey, " → ", pk6_t$ChannelKey
    )

    # Do not use large peak-label boxes: the dots/dashed lines identify maxima,
    # while the subtitle gives exact LASER-DETECTOR identities without obscuring spectra.
    pp4 <- ggplot(df4, aes(ChannelNumber, Value, color = State)) +
      geom_rect(
        data = uv320_band4,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE, fill = "grey85", alpha = .75
      ) +
      geom_vline(
        data = laser_blocks4,
        aes(xintercept = xmin),
        inherit.aes = FALSE, linetype = "dotted", linewidth = .35, alpha = .5
      ) +
      geom_line(linewidth = .75, na.rm = TRUE) +
      geom_point(data = peaks4, size = 2.8) +
      geom_vline(
        data = peaks4,
        aes(xintercept = ChannelNumber, color = State),
        linetype = "dashed", alpha = .5
      ) +
      geom_text(
        data = laser_blocks4,
        aes(x = xmid, y = Inf, label = Laser_label),
        inherit.aes = FALSE, vjust = 1.25, size = 2.8, fontface = "bold"
      ) +
      geom_text(
        data = uv320_band4,
        aes(x = xmid, y = 0.52, label = Missing_label),
        inherit.aes = FALSE, size = 3.0, fontface = "bold", angle = 90
      ) +
      facet_wrap(~ Config, ncol = 1, scales = "fixed") +
      scale_x_continuous(limits = c(axis_min4 - 0.5, axis_max4 + 0.5), expand = c(0, 0)) +
      labs(
        title = paste0("Q1a combined peak QC — ", f, " / ", tr),
        subtitle = paste0(
          subtitle_txt, "\\n", peak_line4,
          "\\nSame x-axis in both panels; grey 320-dUV block = laser absent in 5L. Peak IDs are LASER-DETECTOR."
        ),
        x = "Detector channels grouped by excitation laser (common native axis)",
        y = "Spectral value",
        color = NULL
      ) +
      coord_cartesian(clip = "off") +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "top",
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
      )

    safe <- gsub("[^A-Za-z0-9]+", "_", paste(tr, f, sep = "_"))
    ggsave(
      file.path(Q1A_COMBINED_FIG_DIR, paste0("Q1a_combined_", safe, ".png")),
      pp4, width = 11, height = 7, dpi = 240
    )
  }
}

cat("\nQ1a combined peak-outcome outputs written.\n")
cat("  Separate across-config sheet: Q1a_01_ACROSS_5L_vs_6L_by_condition.csv\n")
cat("  Separate treatment-vs-Unfixed sheet: Q1a_02_TREATMENT_vs_Unfixed_within_5L_6L.csv\n")
cat("  Full combined sheet: Q1a_03_COMBINED_peak_outcomes_one_row_per_fluor_treatment.csv\n")
cat("  Compact combined sheet: Q1a_04_COMBINED_COMPACT_for_Excel.csv\n")
cat("  Category counts: Q1a_05_peak_outcome_category_counts.csv\n")
cat("  Main six-category table: Q1a_00_MAIN_six_category_reporting_table.csv\n")
cat("  Six-category counts: Q1a_00_six_category_counts.csv\n")
cat("  Non-stable combined cases: Q1a_06_COMBINED_NONSTABLE_cases_only.csv\n")
cat("  Four-spectrum QC plots:", Q1A_COMBINED_FIG_DIR, "\n\n")



# ---------------------------------------------------------------------------
# Q1a exact laser/detector map for QC
# ---------------------------------------------------------------------------
# These CSVs show exactly which detector wavelengths belong to each laser block.
# Use the first matched treatment only as the source of channel metadata.
# Channel metadata are instrument metadata and do not depend on treatment.
# active5/active6 contain ChannelKey strings (e.g. "355-432"), NOT channel numbers.
q1a_channel_map <- bind_rows(
  data5[[meta_treatment]]$channels |>
    dplyr::filter(ChannelKey %in% active5) |>
    dplyr::transmute(Config = "5L", ChannelNumber, Laser_nm, Detector_nm, ChannelKey),
  data6[[meta_treatment]]$channels |>
    dplyr::filter(ChannelKey %in% active6) |>
    dplyr::transmute(Config = "6L", ChannelNumber, Laser_nm, Detector_nm, ChannelKey)
) |>
  dplyr::distinct() |>
  dplyr::arrange(Config, ChannelNumber)

write_csv(q1a_channel_map,
          file.path(OUTDIR, "Q1a_peak_shift",
                    "Q1a_channel_map_laser_and_detector_nm.csv"))

q1a_laser_detector_ranges <- q1a_channel_map |>
  group_by(Config, Laser_nm) |>
  summarise(
    Detector_min_nm = min(Detector_nm),
    Detector_max_nm = max(Detector_nm),
    N_detector_channels = n(),
    Detector_nm_values = paste(Detector_nm, collapse = "; "),
    .groups = "drop"
  )

write_csv(q1a_laser_detector_ranges,
          file.path(OUTDIR, "Q1a_peak_shift",
                    "Q1a_laser_detector_ranges.csv"))



# ---------------------------------------------------------------------------
# Q1a PLOTS — separate folders for the two comparison questions
# ---------------------------------------------------------------------------
# Folder 01 = CONFIGURATION comparison, treatment held constant:
#             5L vs 6L separately for Unfixed, 1% PFA, 4% PFA, Perm
#
# Folder 02 = FIXATION comparison, configuration held constant:
#             Unfixed vs 1% PFA, Unfixed vs 4% PFA, Unfixed vs Perm
#             separately for 5L and 6L
#
# All plots use the same clean laser-block / detector-nm axis.

Q1A_PLOT_ROOT <- file.path(OUTDIR, "Q1a_peak_shift", "plots")
Q1A_CONFIG_DIR <- file.path(Q1A_PLOT_ROOT, "01_5L_vs_6L_within_treatment")
Q1A_FIX_DIR <- file.path(Q1A_PLOT_ROOT, "02_treatment_vs_Unfixed")

dir.create(Q1A_CONFIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(Q1A_FIX_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 01. 5L vs 6L within each treatment
# -------------------------
# Re-save the already identified changed configuration cases into a dedicated
# folder so they are not mixed with fixation-vs-Unfixed plots.
if (nrow(q1a_changed) > 0) {
  for (i in seq_len(nrow(q1a_changed))) {
    tr <- q1a_changed$Treatment[i]
    f  <- q1a_changed$Fluor_ID[i]

    df <- bind_rows(
      spectrum_df(data5[[tr]], f, "5L", active5),
      spectrum_df(data6[[tr]], f, "6L", active6)
    ) |>
      filter(is.finite(ChannelNumber), is.finite(Laser_nm), is.finite(Detector_nm))

    peaks <- df |>
      group_by(Config) |>
      slice_max(Value, n = 1, with_ties = FALSE) |>
      ungroup()

    pp <- ggplot(df, aes(ChannelNumber, Value, color = Config)) +
      geom_line(linewidth = .75, na.rm = TRUE) +
      geom_point(data = peaks, size = 3) +
      geom_vline(data = peaks, aes(xintercept = ChannelNumber, color = Config),
                 linetype = "dashed", alpha = .55) +
      labs(
        title = paste0("Q1a — ", f, " / ", tr),
        subtitle = "Configuration comparison: 5L vs 6L with treatment held constant.",
        y = "Spectral value", color = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top")

    pp <- add_laser_detector_axis(pp, df, detector_tick_every = 2)

    safe <- gsub("[^A-Za-z0-9]+", "_", paste(tr, f, sep = "_"))
    ggsave(file.path(Q1A_CONFIG_DIR, paste0("Q1a_5L_vs_6L_", safe, ".png")),
           pp, width = 9, height = 5, dpi = 220)
  }
}

# -------------------------
# 02. Treatment vs Unfixed within each configuration
# -------------------------
# Generate a plot for every fluor/treatment where the identified maximum
# differs from Unfixed. This directly visualizes the fixation comparisons
# represented in QC_06_fixation_peak_changes_vs_Unfixed.csv.

fix_plot_cases <- qc_fixation_peak_flags |>
  filter(Peak_changed_vs_Unfixed %in% TRUE) |>
  distinct(Config, Treatment, Fluor_ID) |>
  mutate(
    Plot_file = paste0(
      "Q1a_treatment_vs_Unfixed_",
      gsub("[^A-Za-z0-9]+", "_", paste(Config, Treatment, Fluor_ID, sep = "_")),
      ".png"
    ),
    Plot_relative_path = file.path(
      "Q1a_peak_shift", "plots", "02_treatment_vs_Unfixed",
      Config, Plot_file
    )
  )

if (nrow(fix_plot_cases) > 0) {
  for (i in seq_len(nrow(fix_plot_cases))) {
    cfg <- fix_plot_cases$Config[i]
    tr  <- fix_plot_cases$Treatment[i]
    f   <- fix_plot_cases$Fluor_ID[i]

    dat <- if (cfg == "5L") data5 else data6
    active <- if (cfg == "5L") active5 else active6

    df <- bind_rows(
      spectrum_df(dat[["Unfixed"]], f, "Unfixed", active) |>
        rename(Condition = Config),
      spectrum_df(dat[[tr]], f, tr, active) |>
        rename(Condition = Config)
    ) |>
      filter(is.finite(ChannelNumber), is.finite(Laser_nm), is.finite(Detector_nm))

    peaks <- df |>
      group_by(Condition) |>
      slice_max(Value, n = 1, with_ties = FALSE) |>
      ungroup()

    pp <- ggplot(df, aes(ChannelNumber, Value, color = Condition)) +
      geom_line(linewidth = .75, na.rm = TRUE) +
      geom_point(data = peaks, size = 3) +
      geom_vline(data = peaks, aes(xintercept = ChannelNumber, color = Condition),
                 linetype = "dashed", alpha = .55) +
      labs(
        title = paste0("Q1a — ", f, " / ", cfg, " / ", tr, " vs Unfixed"),
        subtitle = paste0(
          "Fixation comparison within ", cfg,
          ": Unfixed → ", tr, "."
        ),
        y = "Spectral value", color = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top")

    pp <- add_laser_detector_axis(pp, df, detector_tick_every = 2)

    cfg_dir <- file.path(Q1A_FIX_DIR, cfg)
    dir.create(cfg_dir, recursive = TRUE, showWarnings = FALSE)

    plot_file <- fix_plot_cases$Plot_file[i]
    ggsave(file.path(cfg_dir, plot_file),
           pp, width = 9, height = 5, dpi = 220)
  }
}


# Add plot locations back to the fixation QC table so a reviewer can go
# directly from a changed row to its corresponding figure.
qc_fixation_peak_flags_with_plots <- qc_fixation_peak_flags |>
  left_join(
    fix_plot_cases |>
      select(Config, Treatment, Fluor_ID, Plot_file, Plot_relative_path),
    by = c("Config", "Treatment", "Fluor_ID")
  )

write_csv(
  qc_fixation_peak_flags_with_plots,
  file.path(QC_DIR, "QC_06_fixation_peak_changes_vs_Unfixed.csv")
)

write_csv(
  fix_plot_cases |>
    arrange(Config, factor(Treatment, levels = treatment_order_qc), Fluor_ID),
  file.path(OUTDIR, "Q1a_peak_shift",
            "Q1a_CHANGED_treatment_cases_with_plot_locations.csv")
)


# ---------------------------------------------------------------------------
# FINAL 01 / 02 TABLES — add plot path only to rows that have a changed-case plot
# ---------------------------------------------------------------------------

# 01: 5L vs 6L within treatment
config_plot_lookup <- q1a_changed |>
  distinct(Treatment, Fluor_ID) |>
  mutate(
    Plot_file = paste0(
      "Q1a_5L_vs_6L_",
      gsub("[^A-Za-z0-9]+", "_", paste(Treatment, Fluor_ID, sep = "_")),
      ".png"
    ),
    Plot_relative_path = file.path(
      "Q1a_peak_shift", "plots", "01_5L_vs_6L_within_treatment",
      Plot_file
    )
  )

q1a_across_config_sheet_with_plots <- q1a_across_config_sheet |>
  left_join(
    config_plot_lookup |>
      select(Treatment, Fluor_ID, Plot_file, Plot_relative_path),
    by = c("Treatment", "Fluor_ID")
  ) |>
  # Path must only appear for rows classified as changed.
  mutate(
    Plot_file = if_else(Peak_differs_5L_vs_6L %in% TRUE, Plot_file, NA_character_),
    Plot_relative_path = if_else(
      Peak_differs_5L_vs_6L %in% TRUE, Plot_relative_path, NA_character_
    )
  )

write_csv(
  q1a_across_config_sheet_with_plots,
  file.path(Q1A_COMBINED_DIR, "Q1a_01_ACROSS_5L_vs_6L_by_condition.csv")
)

# 02: treatment vs Unfixed within configuration
q1a_within_config_sheet_with_plots <- q1a_within_config_sheet |>
  left_join(
    fix_plot_cases |>
      select(Config, Treatment, Fluor_ID, Plot_file, Plot_relative_path),
    by = c("Config", "Treatment", "Fluor_ID")
  ) |>
  # Path must only appear for rows classified as changed.
  mutate(
    Plot_file = if_else(Peak_changed_vs_Unfixed %in% TRUE, Plot_file, NA_character_),
    Plot_relative_path = if_else(
      Peak_changed_vs_Unfixed %in% TRUE, Plot_relative_path, NA_character_
    )
  )

write_csv(
  q1a_within_config_sheet_with_plots,
  file.path(Q1A_COMBINED_DIR, "Q1a_02_TREATMENT_vs_Unfixed_within_5L_6L.csv")
)

cat("Updated combined 01 and 02 CSVs with plot paths for changed rows only.\n")

cat("\nQ1a plot folders:\n")
cat("  01 configuration: ", Q1A_CONFIG_DIR, "\n", sep = "")
cat("  02 fixation:      ", Q1A_FIX_DIR, "\n", sep = "")
cat("      ├─ 5L\n")
cat("      └─ 6L\n\n")


# -------------------------
# 5. Q1b — spectral distinguishability of fluorochrome pairs and its dependence on 5L vs 6L
#
# Q1b uses ALL unique pairs of different fluorochrome identities. Same-fluor
# replicate-vs-replicate pairs are excluded because they assess reproducibility,
# not distinguishability. Cosine similarity is the PRIMARY metric; Pearson is
# retained as a supporting/sensitivity metric. Spearman is intentionally omitted.
#
# Cross-configuration comparisons use only LIVE SHARED channels so the zero-filled
# 320-nm rows in 5L are not treated as real signal.
# -------------------------
metrics <- c("pearson", "cosine")
q1b_list <- list()
for (tr in common_treatments) {
  p6 <- restrict_profile(data6[[tr]], shared_active, common_fluors)
  p5 <- restrict_profile(data5[[tr]], shared_active, common_fluors)
  x6 <- t(p6$matrix); x5 <- t(p5$matrix)
  for (metric in metrics) {
    s6 <- similarity_matrix(x6, metric); s5 <- similarity_matrix(x5, metric)
    a <- pair_frame(s6, tr, metric, "6L") |> select(-Config) |> rename(Similarity_6L = Similarity)
    b <- pair_frame(s5, tr, metric, "5L") |> select(-Config) |> rename(Similarity_5L = Similarity)
    q1b_list[[length(q1b_list) + 1]] <-
      inner_join(a, b, by = c("Treatment", "Metric", "Fluor_A", "Fluor_B")) |>
      mutate(Delta_similarity_5L_minus_6L = Similarity_5L - Similarity_6L,
             Absolute_change_in_similarity = abs(Delta_similarity_5L_minus_6L),
             Mean_similarity_across_configs = (Similarity_5L + Similarity_6L) / 2,
             Minimum_similarity_across_configs = pmin(Similarity_5L, Similarity_6L))
  }
}
q1b_pairs <- bind_rows(q1b_list) |>
  # Q1b is about whether DIFFERENT fluorochromes can be distinguished.
  # Replicate columns of the same fluorochrome are useful for reproducibility,
  # but are not competing fluorochrome identities and are excluded here.
  mutate(
    Fluor_A_identity = str_remove(Fluor_A, " \\[rep [0-9]+\\]$"),
    Fluor_B_identity = str_remove(Fluor_B, " \\[rep [0-9]+\\]$")
  ) |>
  filter(Fluor_A_identity != Fluor_B_identity) |>
  select(-Fluor_A_identity, -Fluor_B_identity) |>
  arrange(Treatment, Metric, desc(Absolute_change_in_similarity))
write_csv(q1b_pairs,
          file.path(OUTDIR, "Q1b_pairwise_similarity", "Q1b_all_pairwise_similarity_5L_vs_6L.csv"))
write_csv(q1b_pairs,
          file.path(OUTDIR, "Q1b_pairwise_similarity", "Q1b_ALL_distinct_fluorochrome_pairs_5L_vs_6L.csv"))

# Q1b reporting focuses on COSINE similarity. Pearson remains in the complete
# all-pairs export as a supporting/sensitivity metric, but no Q1b plots are made.
#
# The literature reference cutoff is 0.90. Because a similarity cutoff is not an
# absolute biological boundary, several thresholds are summarized rather than
# treating 0.90 as proof that a pair cannot be unmixed.

q1b_cosine_pairs <- q1b_pairs |>
  filter(
    Metric == "cosine",
    is.finite(Similarity_5L),
    is.finite(Similarity_6L)
  ) |>
  mutate(
    Delta_cosine_5L_minus_6L = Similarity_5L - Similarity_6L,
    Abs_delta_cosine = abs(Delta_cosine_5L_minus_6L),
    `>=0.90 in 5L` = Similarity_5L >= 0.90,
    `>=0.90 in 6L` = Similarity_6L >= 0.90,
    `>=0.90 in both` = Similarity_5L >= 0.90 & Similarity_6L >= 0.90,
    `>=0.90 in either` = Similarity_5L >= 0.90 | Similarity_6L >= 0.90,
    Highest_cosine = pmax(Similarity_5L, Similarity_6L)
  ) |>
  select(
    Treatment, Fluor_A, Fluor_B,
    Cosine_5L = Similarity_5L,
    Cosine_6L = Similarity_6L,
    Delta_cosine_5L_minus_6L,
    Abs_delta_cosine,
    `>=0.90 in 5L`, `>=0.90 in 6L`,
    `>=0.90 in both`, `>=0.90 in either`,
    Highest_cosine
  ) |>
  arrange(Treatment, desc(Highest_cosine), desc(`>=0.90 in both`))

write_csv(
  q1b_cosine_pairs,
  file.path(
    OUTDIR, "Q1b_pairwise_similarity",
    "Q1b_ALL_unique_pairs_cosine_filterable.csv"
  )
)

q1b_thresholds <- c(0.90, 0.95, 0.98, 0.99)

q1b_threshold_summary <- bind_rows(lapply(q1b_thresholds, function(thr) {
  q1b_cosine_pairs |>
    group_by(Treatment) |>
    summarise(
      Threshold = thr,
      `Unique pairs >= threshold in 5L` = sum(Cosine_5L >= thr, na.rm = TRUE),
      `Unique pairs >= threshold in 6L` = sum(Cosine_6L >= thr, na.rm = TRUE),
      `Unique pairs >= threshold in both` =
        sum(Cosine_5L >= thr & Cosine_6L >= thr, na.rm = TRUE),
      `Unique pairs >= threshold in either` =
        sum(Cosine_5L >= thr | Cosine_6L >= thr, na.rm = TRUE),
      `Total unique pairs` = n(),
      .groups = "drop"
    )
})) |>
  arrange(Treatment, Threshold)

write_csv(
  q1b_threshold_summary,
  file.path(
    OUTDIR, "Q1b_pairwise_similarity",
    "Q1b_cosine_threshold_counts.csv"
  )
)

# Publication-reference subset (0.90) retained as a convenient CSV, while the
# report itself shows ALL unique pairs in a filterable table.
q1b_reference_090 <- q1b_cosine_pairs |>
  filter(`>=0.90 in either`) |>
  arrange(Treatment, desc(Highest_cosine))

write_csv(
  q1b_reference_090,
  file.path(
    OUTDIR, "Q1b_pairwise_similarity",
    "Q1b_pairs_cosine_ge_0.90_in_either_configuration.csv"
  )
)


# Q1b clustered cosine-similarity heatmaps
# Clustering is descriptive and uses the same literature-reference cosine cutoff
# (0.90) as the threshold summary. It is NOT a rule that only one fluorochrome
# may be selected from each cluster.

q1b_heatmap_dir <- file.path(OUTDIR, "Q1b_pairwise_similarity", "clustered_heatmaps")
dir.create(q1b_heatmap_dir, recursive = TRUE, showWarnings = FALSE)

make_q1b_clustered_heatmap <- function(dat, treatment, config = c("5L", "6L"),
                                       reference_cosine = Q1B_REFERENCE_COSINE) {
  config <- match.arg(config)
  value_col <- if (config == "5L") "Cosine_5L" else "Cosine_6L"

  d <- dat |>
    filter(Treatment == treatment) |>
    select(Fluor_A, Fluor_B, all_of(value_col))

  fluors_all <- sort(unique(c(d$Fluor_A, d$Fluor_B)))
  if (length(fluors_all) < 2) return(NULL)

  sim <- matrix(
    NA_real_,
    nrow = length(fluors_all),
    ncol = length(fluors_all),
    dimnames = list(fluors_all, fluors_all)
  )
  diag(sim) <- 1

  for (i in seq_len(nrow(d))) {
    a <- d$Fluor_A[i]
    b <- d$Fluor_B[i]
    v <- d[[value_col]][i]
    if (is.finite(v)) {
      sim[a, b] <- v
      sim[b, a] <- v
    }
  }

  # Keep the largest subset with complete finite pairwise similarities.
  # This avoids silently imputing missing similarities and still permits
  # genuine hierarchical clustering of the usable fluorochromes.
  keep <- rownames(sim)
  repeat {
    sub <- sim[keep, keep, drop = FALSE]
    bad_per_row <- rowSums(!is.finite(sub))
    if (all(bad_per_row == 0) || length(keep) <= 2) break
    drop_name <- names(which.max(bad_per_row))
    keep <- setdiff(keep, drop_name)
  }

  sim_use <- sim[keep, keep, drop = FALSE]
  excluded <- setdiff(fluors_all, keep)

  if (length(keep) < 2 || any(!is.finite(sim_use))) {
    return(list(
      index = tibble(
        Treatment = treatment,
        Configuration = config,
        N_fluorochromes_total = length(fluors_all),
        N_fluorochromes_clustered = 0L,
        N_fluorochromes_excluded = length(fluors_all),
        Excluded_fluorochromes = paste(fluors_all, collapse = "; "),
        Clustered = FALSE,
        Reference_cosine = reference_cosine,
        Plot_relative_path = NA_character_
      ),
      membership = tibble()
    ))
  }

  # IMPORTANT: pmax(0, matrix) drops matrix dimensions in R.
  # Preserve the square similarity matrix explicitly before as.dist().
  dist_matrix <- 1 - sim_use
  dist_matrix[dist_matrix < 0] <- 0
  diag(dist_matrix) <- 0

  stopifnot(
    is.matrix(dist_matrix),
    nrow(dist_matrix) == ncol(dist_matrix),
    nrow(dist_matrix) == length(keep),
    all(is.finite(dist_matrix))
  )

  dist_mat <- stats::as.dist(dist_matrix)
  hc <- stats::hclust(dist_mat, method = "average")
  ord <- hc$order
  sim_ord <- sim_use[ord, ord, drop = FALSE]

  # Reference grouping: distance <= 0.10 corresponds to cosine >= 0.90.
  # Hierarchical cut is descriptive; members of a cluster are not guaranteed
  # to have every pairwise cosine >= 0.90.
  cluster_id <- cutree(hc, h = 1 - reference_cosine)  # h = 1 - cosine cutoff

  membership <- tibble(
    Treatment = treatment,
    Configuration = config,
    Fluorochrome = names(cluster_id),
    Cluster = as.integer(cluster_id)
  ) |>
    group_by(Treatment, Configuration, Cluster) |>
    mutate(Cluster_size = n()) |>
    ungroup() |>
    arrange(Cluster, Fluorochrome)

  hm <- as.data.frame(as.table(sim_ord), stringsAsFactors = FALSE)
  names(hm) <- c("Fluor_row", "Fluor_col", "Cosine")
  hm$Fluor_row <- factor(hm$Fluor_row, levels = rev(rownames(sim_ord)))
  hm$Fluor_col <- factor(hm$Fluor_col, levels = colnames(sim_ord))

  p <- ggplot(hm, aes(x = Fluor_col, y = Fluor_row, fill = Cosine)) +
    geom_tile() +
    scale_fill_viridis_c(
      limits = c(0, 1),
      oob = scales::squish,
      name = "Cosine"
    ) +
    labs(
      title = paste0("Q1b clustered cosine similarity — ", config, " — ", treatment),
      subtitle = paste0(
        "Average-linkage clustering; distance = 1 - cosine. ",
        "Reference cluster cut: cosine ", sprintf("%.2f", reference_cosine),
        if (length(excluded) > 0)
          paste0(". Excluded for incomplete pairwise values: ", length(excluded))
        else ""
      ),
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    )

  outfile <- file.path(
    q1b_heatmap_dir,
    paste0(
      "Q1b_clustered_cosine_heatmap_",
      gsub("[^A-Za-z0-9_-]+", "_", treatment), "_", config, ".png"
    )
  )

  ggsave(outfile, p, width = 13, height = 11, dpi = 300, bg = "white")

  list(
    index = tibble(
      Treatment = treatment,
      Configuration = config,
      N_fluorochromes_total = length(fluors_all),
      N_fluorochromes_clustered = length(keep),
      N_fluorochromes_excluded = length(excluded),
      Excluded_fluorochromes =
        if (length(excluded) > 0) paste(excluded, collapse = "; ") else "",
      Clustered = TRUE,
      Reference_cosine = reference_cosine,
      Plot_relative_path = file.path(
        "Q1b_pairwise_similarity", "clustered_heatmaps", basename(outfile)
      )
    ),
    membership = membership
  )
}

q1b_heatmap_results <- lapply(sort(unique(q1b_cosine_pairs$Treatment)), function(trt) {
  list(
    make_q1b_clustered_heatmap(q1b_cosine_pairs, trt, "5L", 0.90),
    make_q1b_clustered_heatmap(q1b_cosine_pairs, trt, "6L", 0.90)
  )
}) |>
  unlist(recursive = FALSE)

q1b_heatmap_index <- bind_rows(lapply(q1b_heatmap_results, `[[`, "index"))
q1b_cluster_membership <- bind_rows(lapply(q1b_heatmap_results, `[[`, "membership"))

# Summarize how tightly similar each descriptive cluster actually is.
# IMPORTANT: calculate each cluster independently. Capture the current
# treatment/configuration/cluster as scalar values before filtering so dplyr
# does not accidentally compare columns to themselves inside filter().
cluster_keys <- q1b_cluster_membership |>
  distinct(Treatment, Configuration, Cluster) |>
  arrange(Treatment, Configuration, Cluster)

q1b_cluster_pairwise_summary <- bind_rows(lapply(seq_len(nrow(cluster_keys)), function(i) {

  trt_i <- cluster_keys$Treatment[i]
  cfg_i <- cluster_keys$Configuration[i]
  cl_i  <- cluster_keys$Cluster[i]

  members <- q1b_cluster_membership |>
    filter(
      Treatment == trt_i,
      Configuration == cfg_i,
      Cluster == cl_i
    ) |>
    pull(Fluorochrome) |>
    unique()

  n_fluors <- length(members)
  expected_pairs <- choose(n_fluors, 2)

  value_col <- if (cfg_i == "5L") "Cosine_5L" else "Cosine_6L"

  pair_dat <- q1b_cosine_pairs |>
    filter(
      Treatment == trt_i,
      Fluor_A %in% members,
      Fluor_B %in% members
    ) |>
    select(Fluor_A, Fluor_B, all_of(value_col)) |>
    distinct(Fluor_A, Fluor_B, .keep_all = TRUE)

  vals <- pair_dat[[value_col]]
  vals <- vals[is.finite(vals)]

  # Safety check: a cluster of N fluorochromes can never have more than
  # choose(N, 2) unique unordered within-cluster pairs.
  if (nrow(pair_dat) > expected_pairs) {
    stop(
      "Q1b cluster summary has too many pairs for ",
      trt_i, " / ", cfg_i, " / cluster ", cl_i,
      ": found ", nrow(pair_dat),
      " but expected at most ", expected_pairs, "."
    )
  }

  tibble(
    Treatment = trt_i,
    Configuration = cfg_i,
    Cluster = cl_i,
    N_fluorochromes = n_fluors,
    Expected_unique_pairs = expected_pairs,
    N_within_cluster_pairs = length(vals),
    Min_pairwise_cosine =
      if (length(vals) > 0) min(vals) else NA_real_,
    Median_pairwise_cosine =
      if (length(vals) > 0) median(vals) else NA_real_,
    Max_pairwise_cosine =
      if (length(vals) > 0) max(vals) else NA_real_,
    N_pairs_ge_0_90 =
      if (length(vals) > 0) sum(vals >= 0.90) else 0L,
    Percent_pairs_ge_0_90 =
      if (length(vals) > 0) 100 * sum(vals >= 0.90) / length(vals) else NA_real_,
    Fluorochromes = paste(sort(members), collapse = ", ")
  )
})) |>
  arrange(Treatment, Configuration, desc(N_fluorochromes), Cluster)

write_csv(
  q1b_heatmap_index,
  file.path(OUTDIR, "Q1b_pairwise_similarity", "Q1b_clustered_cosine_heatmap_index.csv")
)

write_csv(
  q1b_cluster_membership,
  file.path(OUTDIR, "Q1b_pairwise_similarity", "Q1b_cluster_membership_reference_0.90.csv")
)

write_csv(
  q1b_cluster_pairwise_summary,
  file.path(
    OUTDIR, "Q1b_pairwise_similarity",
    "Q1b_cluster_pairwise_cosine_summary_reference_0.90.csv"
  )
)

# 6. Q2a — fixation effect: each fluorochrome, treatment vs Unfixed
#          (uses each config's own ACTIVE channels)
# -------------------------
treated_conditions <- setdiff(common_treatments, "Unfixed")
if (!"Unfixed" %in% common_treatments) stop("Unfixed not found; Q2 requires it as reference.")

q2a_list <- list()
for (cfg in c("5L", "6L")) {
  dat <- if (cfg == "5L") data5 else data6
  active_keys <- if (cfg == "5L") active5 else active6
  ref <- restrict_profile(dat[["Unfixed"]], active_keys, common_fluors)
  for (tr in treated_conditions) {
    tx <- restrict_profile(dat[[tr]], active_keys, common_fluors)
    for (f in common_fluors) {
      # V3: renamed `t` -> `trt`; the old name masked base::t() inside the loop.
      u <- ref$matrix[, f]; trt <- tx$matrix[, f]
      pku <- peak_info(u, ref$channels); pkt <- peak_info(trt, tx$channels)
      total_u <- sum(u[is.finite(u)], na.rm = TRUE)
      total_t <- sum(trt[is.finite(trt)], na.rm = TRUE)
      dif <- abs(trt - u)
      q2a_list[[length(q2a_list) + 1]] <- tibble(
        Config = cfg, Treatment = tr, Fluor_ID = f,
        RMSE = rmse(trt, u),
        Pearson = safe_cor(trt, u, "pearson"),
        Cosine = cosine(trt, u),
        Max_abs_detector_difference = if (any(is.finite(dif)))
                                        max(dif, na.rm = TRUE) else NA_real_,
        # V3: absolute-intensity metrics are only emitted when the source scale
        # supports them (see 00_scaling_audit.csv). With mixed-scale workbooks or
        # renormalized spectra they would report a scaling artefact as biology.
        Total_intensity_unfixed = if (INTENSITY_METRICS_VALID) total_u else NA_real_,
        Total_intensity_treated = if (INTENSITY_METRICS_VALID) total_t else NA_real_,
        Total_intensity_change  = if (INTENSITY_METRICS_VALID) total_t - total_u else NA_real_,
        Total_intensity_change_percent = if (!INTENSITY_METRICS_VALID) NA_real_
          else ifelse(total_u == 0, NA_real_, 100 * (total_t - total_u) / total_u),
        Intensity_metric_status = if (INTENSITY_METRICS_VALID) "Valid" else INTENSITY_INVALID_REASON,
        Peak_unfixed = pku$Peak_ChannelKey, Peak_treated = pkt$Peak_ChannelKey,
        Peak_unfixed_nm = pku$Peak_Detector_nm, Peak_treated_nm = pkt$Peak_Detector_nm,
        Peak_shift_nm = pkt$Peak_Detector_nm - pku$Peak_Detector_nm,
        Abs_peak_shift_nm = abs(pkt$Peak_Detector_nm - pku$Peak_Detector_nm))
    }
  }
}
q2a_effects <- bind_rows(q2a_list)
write_csv(q2a_effects, file.path(OUTDIR, "Q2_fixation", "Q2a_all_fixation_effects.csv"))

q2a_treatment_summary <- q2a_effects |>
  group_by(Config, Treatment) |>
  summarise(N_fluorochromes = n(),
            Median_RMSE = median(RMSE, na.rm = TRUE), Mean_RMSE = mean(RMSE, na.rm = TRUE),
            Median_Pearson = median(Pearson, na.rm = TRUE),
            Median_Cosine = median(Cosine, na.rm = TRUE),
            Median_max_abs_detector_difference = median(Max_abs_detector_difference, na.rm = TRUE),
            Median_abs_intensity_change_percent = median(abs(Total_intensity_change_percent), na.rm = TRUE),
            Median_abs_peak_shift_nm = median(Abs_peak_shift_nm, na.rm = TRUE),
            .groups = "drop")
write_csv(q2a_treatment_summary, file.path(OUTDIR, "Q2_fixation", "Q2a_treatment_summary.csv"))

# ---------------------------------------------------------------------------
# Q2a INFERENCE — "is any treatment worse?" answered with a paired test (V3)
# ---------------------------------------------------------------------------
# The collaborator question "any worse?" is comparative, but the previous
# version answered it only by eyeballing overlapping boxplots. Every treatment
# is measured on the SAME set of fluorochromes, so fluorochrome is a natural
# blocking factor and a paired/blocked test is available without any extra data.
#
# This is inference across fluorochromes within one bead experiment. It answers
# "is treatment A systematically worse than treatment B across this dye panel?"
# It does NOT answer "is this reproducible across experiments" -- that still
# requires running the same workflow on the other bead runs. Stated explicitly
# in the output so the distinction is not lost downstream.

q2a_rmse_wide <- q2a_effects |>
  select(Config, Treatment, Fluor_ID, RMSE) |>
  filter(is.finite(RMSE)) |>
  pivot_wider(names_from = Treatment, values_from = RMSE) |>
  tidyr::drop_na()

q2a_friedman <- bind_rows(lapply(unique(q2a_rmse_wide$Config), function(cfg) {
  m <- q2a_rmse_wide |> filter(Config == cfg) |> select(-Config, -Fluor_ID) |> as.matrix()
  if (nrow(m) < 3 || ncol(m) < 3)
    return(tibble(Config = cfg, Test = "Friedman", N_fluorochromes = nrow(m),
                  Statistic = NA_real_, df = NA_real_, P_value = NA_real_,
                  Kendall_W = NA_real_))
  ft <- stats::friedman.test(m)
  # Kendall's W = Friedman chi-square / (n * (k - 1)); effect size on 0-1.
  W <- unname(ft$statistic) / (nrow(m) * (ncol(m) - 1))
  tibble(Config = cfg, Test = "Friedman (RMSE blocked by fluorochrome)",
         N_fluorochromes = nrow(m), Statistic = unname(ft$statistic),
         df = unname(ft$parameter), P_value = ft$p.value, Kendall_W = W)
}))
write_csv(q2a_friedman, file.path(OUTDIR, "Q2_fixation", "Q2a_omnibus_treatment_test.csv"))

# Pairwise follow-up: Wilcoxon signed-rank on paired per-fluorochrome RMSE,
# with the matched-pairs rank-biserial correlation as the effect size and a
# Benjamini-Hochberg adjustment within each configuration.
q2a_pairwise <- bind_rows(lapply(unique(q2a_rmse_wide$Config), function(cfg) {
  d <- q2a_rmse_wide |> filter(Config == cfg)
  trts <- setdiff(names(d), c("Config", "Fluor_ID"))
  cmb <- utils::combn(trts, 2, simplify = FALSE)
  bind_rows(lapply(cmb, function(pr) {
    a <- d[[pr[1]]]; b <- d[[pr[2]]]
    ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]
    if (length(a) < 5) return(tibble(
      Config = cfg, Treatment_A = pr[1], Treatment_B = pr[2],
      N_pairs = length(a), Median_RMSE_A = NA_real_, Median_RMSE_B = NA_real_,
      Median_paired_difference = NA_real_, Rank_biserial_r = NA_real_,
      P_value = NA_real_))
    wt <- suppressWarnings(stats::wilcox.test(a, b, paired = TRUE, exact = FALSE))
    dd <- a - b; dd <- dd[dd != 0]
    r  <- if (length(dd)) {
      rk <- rank(abs(dd))
      (sum(rk[dd > 0]) - sum(rk[dd < 0])) / sum(rk)
    } else NA_real_
    tibble(Config = cfg, Treatment_A = pr[1], Treatment_B = pr[2],
           N_pairs = length(a),
           Median_RMSE_A = median(a), Median_RMSE_B = median(b),
           Median_paired_difference = median(a - b),
           Rank_biserial_r = r, P_value = wt$p.value)
  }))
})) |>
  group_by(Config) |>
  mutate(P_adj_BH = p.adjust(P_value, method = "BH")) |>
  ungroup() |>
  mutate(
    Direction = case_when(
      !is.finite(P_adj_BH) ~ "Not testable",
      P_adj_BH >= 0.05 ~ "No detectable difference",
      Median_paired_difference > 0 ~ paste0(Treatment_A, " larger spectral change"),
      Median_paired_difference < 0 ~ paste0(Treatment_B, " larger spectral change"),
      TRUE ~ "No detectable difference"
    ),
    Scope = "Across fluorochromes within this bead experiment; not cross-experiment replication"
  ) |>
  arrange(Config, P_adj_BH)

write_csv(q2a_pairwise,
          file.path(OUTDIR, "Q2_fixation", "Q2a_pairwise_treatment_tests.csv"))

cat("\nQ2a — paired treatment comparisons (Wilcoxon signed-rank on per-fluorochrome RMSE):\n")
print(q2a_pairwise |> select(Config, Treatment_A, Treatment_B, Median_paired_difference,
                             Rank_biserial_r, P_adj_BH, Direction), n = Inf)

# ---- Q2a RMSE plot: label the 5 fluorochromes with the largest RMSE ----
# Top 5 are selected SEPARATELY within each Configuration x Treatment.
q2a_top5_rmse <- q2a_effects |>
  filter(is.finite(RMSE)) |>
  group_by(Config, Treatment) |>
  slice_max(order_by = RMSE, n = 5, with_ties = FALSE) |>
  ungroup()

write_csv(
  q2a_top5_rmse |>
    select(Config, Treatment, Fluor_ID, RMSE, Pearson, Cosine,
           Max_abs_detector_difference, Total_intensity_change_percent,
           Peak_shift_nm, Abs_peak_shift_nm) |>
    arrange(Config, Treatment, desc(RMSE)),
  file.path(OUTDIR, "Q2_fixation",
            "Q2a_top5_highest_RMSE_each_treatment_config.csv")
)

p_q2a_rmse <- ggplot(q2a_effects, aes(Treatment, RMSE)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15, alpha = .35, size = 1) +
  # Circle the five largest fixation effects in each treatment/configuration.
  geom_point(
    data = q2a_top5_rmse,
    shape = 21, size = 2.5, stroke = .7, fill = "white"
  ) +
  # Label those same five fluorochromes. ggrepel reduces label overlap.
  ggrepel::geom_text_repel(
    data = q2a_top5_rmse,
    aes(label = Fluor_ID),
    size = 2.7,
    box.padding = .35,
    point.padding = .25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 1418,
    show.legend = FALSE
  ) +
  facet_wrap(~ Config, scales = "free_y") +
  labs(
    title = "Q2a — What is the overall impact of each fixation treatment?",
    subtitle = "RMSE compares treated vs Unfixed within the same configuration; five highest-RMSE fluorochromes are labelled per treatment × configuration",
    x = NULL,
    y = "RMSE"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(OUTDIR, "figures", "Q2a_RMSE_by_treatment_and_config.png"),
  p_q2a_rmse,
  width = 12,
  height = 7,
  dpi = 300
)

cat("\nQ2a — top 5 highest-RMSE fluorochromes within each treatment x configuration:\n")
print(
  q2a_top5_rmse |>
    select(Config, Treatment, Fluor_ID, RMSE) |>
    arrange(Config, Treatment, desc(RMSE)),
  n = Inf
)

ggsave(file.path(OUTDIR, "figures", "Q2a_RMSE_heatmap.png"),
  ggplot(q2a_effects, aes(Treatment, Fluor_ID, fill = RMSE)) +
    geom_tile() + facet_wrap(~ Config) + scale_fill_viridis_c() +
    labs(title = "Q2a — Magnitude of fixation-related spectral change by fluorochrome", x = NULL, y = "Fluorochrome", fill = "RMSE") +
    theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 4)),
  width = 11, height = 18, dpi = 200)

# -------------------------
# 7. Q2b — susceptibility rankings (per metric, no combined score)
# -------------------------
q2b_rankings <- q2a_effects |>
  group_by(Config, Fluor_ID) |>
  summarise(Mean_RMSE = mean(RMSE, na.rm = TRUE), Median_RMSE = median(RMSE, na.rm = TRUE),
            Mean_Pearson = mean(Pearson, na.rm = TRUE),
            Mean_Cosine = mean(Cosine, na.rm = TRUE),
            Mean_max_abs_detector_difference = mean(Max_abs_detector_difference, na.rm = TRUE),
            Mean_abs_intensity_change_percent = mean(abs(Total_intensity_change_percent), na.rm = TRUE),
            Mean_abs_peak_shift_nm = mean(Abs_peak_shift_nm, na.rm = TRUE), .groups = "drop") |>
  group_by(Config) |>
  mutate(Rank_RMSE = rank(-Mean_RMSE, ties.method = "min"),
         Rank_Pearson_disruption = rank(Mean_Pearson, ties.method = "min"),
         Rank_Cosine_disruption = rank(Mean_Cosine, ties.method = "min"),
         Rank_Max_difference = rank(-Mean_max_abs_detector_difference, ties.method = "min"),
         Rank_Intensity_change = rank(-Mean_abs_intensity_change_percent, ties.method = "min"),
         Rank_Peak_shift = rank(-Mean_abs_peak_shift_nm, ties.method = "min")) |> ungroup()
write_csv(q2b_rankings,
          file.path(OUTDIR, "Q2_fixation", "Q2b_fluorochrome_susceptibility_rankings.csv"))

# ---------------------------------------------------------------------------
# Q2b — FLUOROCHROME TYPE / FAMILY (V4: manually curated table replaces the
# name-based regex heuristic entirely)
# ---------------------------------------------------------------------------
# The collaborator aim is "type of fluor more susceptible?". V3's regex
# heuristic (classify_fluorochrome()) got the family/tandem status wrong for
# ~90 of 144 dyes on this panel -- most importantly it missed several PE/APC
# Fire and Dazzle tandems (calling them bare protein), merged NovaFluor with
# BD Real and eFluor with cFluor into single bins, and left several dyes
# (Cy2/3/5, Texas Red, Live-or-Dye, a misspelled StarBright entry) Unassigned.
# V4 removes the heuristic from the analysis path: every dye on this panel is
# assigned by a human-reviewed lookup table (FAMILY_FILE), checked in
# alongside this script. There is no silent regex fallback -- if the table is
# missing, the script stops rather than quietly reverting to the old
# annotation.
#
# FAMILY_FILE columns (one row per Fluor_ID):
#   Fluor_ID, Fluorochrome, Family, Structural_class, Is_tandem,
#   Needs_review, Review_note, Annotation_source
# Needs_review == TRUE marks dyes with an unresolved naming ambiguity
# (probable duplicates such as byg584/cfluor-yg584, or an unverifiable
# catalogue name such as af561). Those rows are kept in the per-fluorochrome
# outputs but excluded from the family-level Kruskal-Wallis test, since a
# citable family comparison shouldn't rest on an unresolved identity.
# Any dye without an assignable family is "Unassigned" and is excluded from
# both, exactly as before.

if (!file.exists(FAMILY_FILE)) {
  stop("Q2b: curated fluorochrome family table not found at '", FAMILY_FILE,
       "'. This file is required -- there is no automatic heuristic fallback ",
       "in V4. Restore Q2b_fluorochrome_family_annotation_CURATED.csv next to ",
       "this script (see the column spec in the comment above this check).")
}

q2b_family <- read_csv(FAMILY_FILE, show_col_types = FALSE)

# Needs_review/Review_note are optional in FAMILY_FILE -- add them if absent,
# checked on the raw column names (not inside mutate(), where referencing a
# column that doesn't exist would error before any conditional runs).
if (!"Needs_review" %in% names(q2b_family)) q2b_family$Needs_review <- FALSE
if (!"Review_note"  %in% names(q2b_family)) q2b_family$Review_note  <- NA_character_

q2b_family <- q2b_family |>
  mutate(
    Is_tandem = as.logical(Is_tandem),
    Needs_review = as.logical(Needs_review),
    Annotation_source = "User-curated (manual)"
  )

missing_from_table <- setdiff(unique(q2a_effects$Fluor_ID), q2b_family$Fluor_ID)
if (length(missing_from_table) > 0) {
  warning("Q2b: ", length(missing_from_table),
          " Fluor_ID(s) in the data have no row in ", FAMILY_FILE,
          " and will be dropped from every Q2b output: ",
          paste(missing_from_table, collapse = ", "))
}

n_flagged <- sum(q2b_family$Needs_review, na.rm = TRUE)
if (n_flagged > 0) {
  message("Q2b: ", n_flagged,
          " dye(s) are flagged Needs_review in the curated table (kept in ",
          "per-fluorochrome outputs, excluded from the family-level test): ",
          paste(q2b_family$Fluor_ID[which(q2b_family$Needs_review)], collapse = ", "))
}

q2b_by_family <- q2a_effects |>
  inner_join(q2b_family |> select(Fluor_ID, Family, Structural_class, Is_tandem),
             by = "Fluor_ID") |>
  filter(Family != "Unassigned") |>
  group_by(Config, Treatment, Family) |>
  summarise(N_fluorochromes = n(),
            Median_RMSE = median(RMSE, na.rm = TRUE),
            IQR_RMSE = IQR(RMSE, na.rm = TRUE),
            Median_Cosine = median(Cosine, na.rm = TRUE),
            Median_abs_peak_shift_nm = median(Abs_peak_shift_nm, na.rm = TRUE),
            .groups = "drop") |>
  mutate(Annotation_source = "User-curated (manual)") |>
  arrange(Config, Treatment, desc(Median_RMSE))

write_csv(q2b_by_family,
          file.path(OUTDIR, "Q2_fixation", "Q2b_susceptibility_by_family.csv"))

# Group-level test: do families differ in fixation susceptibility? Kruskal-Wallis
# within each configuration x treatment, on per-fluorochrome RMSE. Families with
# fewer than 3 members are dropped -- a "family" of one dye is a dye, not a family.
q2b_family_test <- q2a_effects |>
  inner_join(q2b_family |> select(Fluor_ID, Family, Needs_review), by = "Fluor_ID") |>
  filter(Family != "Unassigned", is.finite(RMSE), !Needs_review) |>
  group_by(Config, Treatment, Family) |>
  filter(n() >= 3) |>
  ungroup() |>
  group_by(Config, Treatment) |>
  group_modify(~ {
    if (dplyr::n_distinct(.x$Family) < 2)
      return(tibble(N_families = dplyr::n_distinct(.x$Family), N_fluorochromes = nrow(.x),
                    Chi_squared = NA_real_, df = NA_real_, P_value = NA_real_,
                    Epsilon_squared = NA_real_))
    kt <- stats::kruskal.test(RMSE ~ factor(Family), data = .x)
    n <- nrow(.x); k <- dplyr::n_distinct(.x$Family)
    tibble(N_families = k, N_fluorochromes = n,
           Chi_squared = unname(kt$statistic), df = unname(kt$parameter),
           P_value = kt$p.value,
           Epsilon_squared = unname(kt$statistic) / ((n^2 - 1) / (n + 1)))
  }) |>
  ungroup() |>
  mutate(P_adj_BH = p.adjust(P_value, method = "BH"),
         Annotation_source = "User-curated (manual); Needs_review dyes excluded")

write_csv(q2b_family_test,
          file.path(OUTDIR, "Q2_fixation", "Q2b_family_level_test.csv"))

p_q2b_family <- q2a_effects |>
  inner_join(q2b_family |> select(Fluor_ID, Family), by = "Fluor_ID") |>
  filter(Family != "Unassigned", is.finite(RMSE)) |>
  group_by(Family) |> filter(n() >= 3) |> ungroup() |>
  mutate(Family = forcats::fct_reorder(Family, RMSE, .fun = median)) |>
  ggplot(aes(RMSE, Family)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = .18, alpha = .3, size = .9) +
  facet_grid(Config ~ Treatment, scales = "free_x") +
  labs(title = "Q2b — Is a fluorochrome type more susceptible to fixation?",
       subtitle = "Per-fluorochrome RMSE vs Unfixed, grouped by family (user-curated annotation)",
       x = "RMSE (treated vs Unfixed)", y = NULL) +
  theme_bw(base_size = 10)

ggsave(file.path(OUTDIR, "figures", "Q2b_RMSE_by_family.png"),
       p_q2b_family, width = 12, height = 7, dpi = 250)

# -------------------------
# 8. Q2c — configuration-dependent fixation effect (6L vs 5L)
# -------------------------
q2c <- q2a_effects |>
  select(Config, Treatment, Fluor_ID, RMSE, Pearson, Cosine,
         Max_abs_detector_difference, Total_intensity_change_percent,
         Peak_shift_nm, Abs_peak_shift_nm) |>
  pivot_wider(names_from = Config,
              values_from = c(RMSE, Pearson, Cosine,
                              Max_abs_detector_difference, Total_intensity_change_percent,
                              Peak_shift_nm, Abs_peak_shift_nm), names_sep = "_") |>
  mutate(
    Delta_RMSE_6L_minus_5L = RMSE_6L - RMSE_5L,
    Delta_Pearson_6L_minus_5L = Pearson_6L - Pearson_5L,
    Delta_Cosine_6L_minus_5L = Cosine_6L - Cosine_5L,
    Delta_Max_difference_6L_minus_5L = Max_abs_detector_difference_6L - Max_abs_detector_difference_5L,
    # NOTE: signed intensity delta — a sign flip between configs shows as a large value
    Delta_Intensity_change_percent_6L_minus_5L = Total_intensity_change_percent_6L - Total_intensity_change_percent_5L,
    Delta_Peak_shift_nm_6L_minus_5L = Peak_shift_nm_6L - Peak_shift_nm_5L)
write_csv(q2c, file.path(OUTDIR, "Q2_fixation", "Q2c_5L_vs_6L_fixation_effects.csv"))

q2c_ranked <- q2c |>
  mutate(Abs_Delta_RMSE = abs(Delta_RMSE_6L_minus_5L),
         Abs_Delta_Pearson = abs(Delta_Pearson_6L_minus_5L),
         Abs_Delta_Cosine = abs(Delta_Cosine_6L_minus_5L),
         Abs_Delta_Max_difference = abs(Delta_Max_difference_6L_minus_5L),
         Abs_Delta_Intensity = abs(Delta_Intensity_change_percent_6L_minus_5L),
         Abs_Delta_Peak_shift = abs(Delta_Peak_shift_nm_6L_minus_5L))
write_csv(q2c_ranked,
          file.path(OUTDIR, "Q2_fixation", "Q2c_configuration_dependence_ranked_source.csv"))

# ---- Q2c RMSE plot: label the 5 fluorochromes farthest from the diagonal ----
# Top 5 are selected SEPARATELY within each treatment using:
# |RMSE_6L - RMSE_5L|
q2c_top5_rmse_difference <- q2c_ranked |>
  filter(is.finite(RMSE_5L), is.finite(RMSE_6L), is.finite(Abs_Delta_RMSE)) |>
  group_by(Treatment) |>
  slice_max(order_by = Abs_Delta_RMSE, n = 5, with_ties = FALSE) |>
  ungroup()

write_csv(
  q2c_top5_rmse_difference |>
    select(Treatment, Fluor_ID, RMSE_5L, RMSE_6L,
           Delta_RMSE_6L_minus_5L, Abs_Delta_RMSE) |>
    arrange(Treatment, desc(Abs_Delta_RMSE)),
  file.path(OUTDIR, "Q2_fixation",
            "Q2c_top5_largest_RMSE_difference_each_treatment.csv")
)

p_q2c_rmse <- ggplot(q2c_ranked, aes(RMSE_5L, RMSE_6L)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(color = Abs_Delta_RMSE), alpha = .7) +
  geom_point(
    data = q2c_top5_rmse_difference,
    shape = 21, size = 2.8, stroke = .8, fill = "white", color = "black"
  ) +
  ggrepel::geom_text_repel(
    data = q2c_top5_rmse_difference,
    aes(label = Fluor_ID),
    size = 2.7,
    box.padding = .35,
    point.padding = .25,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 1418,
    show.legend = FALSE
  ) +
  facet_wrap(~ Treatment, scales = "free") +
  scale_color_viridis_c(name = "Absolute\nRMSE difference") +
  labs(
    title = "Q2c — Does the impact of fixation differ between 5L and 6L?",
    subtitle = "Most points can be assessed against the 1:1 line; five largest |RMSE 6L - RMSE 5L| differences are labelled per treatment",
    x = "Fixation RMSE (5L)",
    y = "Fixation RMSE (6L)"
  ) +
  theme_bw(base_size = 11)

ggsave(
  file.path(OUTDIR, "figures", "Q2c_RMSE_5L_vs_6L.png"),
  p_q2c_rmse,
  width = 10, height = 6, dpi = 300
)

cat("\nQ2c — top 5 fluorochromes farthest from the diagonal within each treatment:\n")
print(
  q2c_top5_rmse_difference |>
    select(Treatment, Fluor_ID, RMSE_5L, RMSE_6L,
           Delta_RMSE_6L_minus_5L, Abs_Delta_RMSE) |>
    arrange(Treatment, desc(Abs_Delta_RMSE)),
  n = Inf
)

# -------------------------
# 8b. Normalization sensitivity (V3)
# -------------------------
# Because the magnitude metrics now depend on a normalization choice, that
# choice must be shown not to drive the conclusions. Q2a is recomputed under the
# alternative normalization and the two rankings are compared. If the treatment
# ordering and the fluorochrome-level ranking are stable, the result is a
# property of the data rather than of the preprocessing.
if (!is.na(PREPROCESS_SENSITIVITY) && PREPROCESS_SENSITIVITY != PREPROCESS) {
  PREPROCESS_MAIN <- PREPROCESS
  PREPROCESS <- PREPROCESS_SENSITIVITY   # restrict_profile() reads this global

  alt_list <- list()
  for (cfg in c("5L", "6L")) {
    dat <- if (cfg == "5L") data5 else data6
    active_keys <- if (cfg == "5L") active5 else active6
    ref <- restrict_profile(dat[["Unfixed"]], active_keys, common_fluors)
    for (tr in treated_conditions) {
      tx <- restrict_profile(dat[[tr]], active_keys, common_fluors)
      for (f in common_fluors)
        alt_list[[length(alt_list) + 1]] <- tibble(
          Config = cfg, Treatment = tr, Fluor_ID = f,
          RMSE_alt = rmse(tx$matrix[, f], ref$matrix[, f]))
    }
  }
  PREPROCESS <- PREPROCESS_MAIN

  q2a_sensitivity <- q2a_effects |>
    select(Config, Treatment, Fluor_ID, RMSE) |>
    inner_join(bind_rows(alt_list), by = c("Config", "Treatment", "Fluor_ID"))

  q2a_sensitivity_summary <- q2a_sensitivity |>
    group_by(Config, Treatment) |>
    summarise(
      N_fluorochromes = n(),
      Spearman_rank_agreement = suppressWarnings(
        cor(RMSE, RMSE_alt, method = "spearman", use = "complete.obs")),
      Top10_overlap = length(intersect(
        Fluor_ID[order(-RMSE)][1:min(10, n())],
        Fluor_ID[order(-RMSE_alt)][1:min(10, n())])),
      .groups = "drop") |>
    mutate(Primary_normalization = PREPROCESS,
           Alternative_normalization = PREPROCESS_SENSITIVITY)

  write_csv(q2a_sensitivity,
            file.path(OUTDIR, "Q2_fixation", "Q2a_normalization_sensitivity_per_fluorochrome.csv"))
  write_csv(q2a_sensitivity_summary,
            file.path(OUTDIR, "Q2_fixation", "Q2a_normalization_sensitivity_summary.csv"))

  cat("\nNormalization sensitivity (primary vs alternative):\n")
  print(q2a_sensitivity_summary)
}

# -------------------------
# 9. Final summary
# -------------------------
cat("\n================ ANALYSIS COMPLETE ================\n")
cat("Results folder:", normalizePath(OUTDIR), "\n")
cat("Preprocessing:", PREPROCESS, "\n")
cat("Fluorochromes analyzed:", length(common_fluors), "\n")
cat("Treatments:", paste(common_treatments, collapse = ", "), "\n")
cat("Live shared channels (Q1b/Q2c):", length(shared_active),
    "| UV lost at 5L:", length(uv_lost_keys), "\n")
cat("Q1b pair rows:", nrow(q1b_pairs), "| Q2a effect rows:", nrow(q2a_effects), "\n")
cat("Native scale consistent across sheets:", scale_is_consistent, "\n")
cat("Absolute-intensity metrics emitted:", INTENSITY_METRICS_VALID, "\n")
cat("Q2b family annotation: USER-CURATED (manual table, ", FAMILY_FILE, "); ",
    n_flagged, " dye(s) flagged Needs_review and excluded from the family-level test\n", sep = "")
cat("===================================================\n")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
