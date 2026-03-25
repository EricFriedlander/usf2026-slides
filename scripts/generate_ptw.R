# generate_ptw.R
# Generates per-species PTW-aligned parquet files from combinedWide.parquet.
# Run once from the usf2026-slides/ directory before rendering the slides.
#
# Adapted from: https://github.com/aideneve15/Dadabay_Research/blob/main/scripts/importing/Align.R
# Fixes: each species now writes to its own output file (original had a copy-paste bug)
# Reference: uses time_point=30, sample=m4, enzyme=CYP+UGT; falls back to first row

library(arrow)
library(tidyverse)
library(ptw)
library(here)

# ── helpers ───────────────────────────────────────────────────────────────────

split_by_species <- function(df, species_col = "species") {
  split(df, df[[species_col]])
}

wide_to_long <- function(df) {
  df |>
    pivot_longer(
      cols = where(is.numeric),
      names_to = "time",
      values_to = "intensity"
    )
}

correct_baseline <- function(df, l = 1e9, maxit = 25, meta_in = TRUE) {
  if (meta_in) {
    meta <- df |> dplyr::select(where(~ !is.numeric(.)))
  }
  numeric_df <- df |> dplyr::select(where(is.numeric))
  numeric_corrected <- baseline.corr(numeric_df, lambda = l, p = 0.001,
                                     eps = 1e-8, maxit = maxit)
  if (meta_in) {
    return(bind_cols(meta, numeric_corrected))
  }
  return(numeric_corrected)
}

ptw_function <- function(df) {
  meta <- df |>
    dplyr::select(time_point, species, sample, enzyme, instrument)

  # Use a biologically meaningful reference: max incubation, m4, combined enzyme.
  # Fall back to first row if that sample isn't present for this species.
  reference <- df |>
    dplyr::filter(time_point == 30, sample == "m4", enzyme == "CYP+UGT") |>
    dplyr::slice(1)

  if (nrow(reference) == 0) {
    message("Reference sample not found for species '",
            unique(df$species), "'; using first row as reference.")
    reference <- df |> dplyr::slice(1)
  }

  samples <- df |>
    dplyr::anti_join(
      reference |> dplyr::select(time_point, species, sample, enzyme, instrument),
      by = c("time_point", "species", "sample", "enzyme", "instrument")
    )

  reference_meta <- reference |> dplyr::select(time_point, species, sample, enzyme, instrument)
  samples_meta   <- samples   |> dplyr::select(time_point, species, sample, enzyme, instrument)

  reference_num <- reference |> dplyr::select(where(is.numeric))
  samples_num   <- samples   |> dplyr::select(where(is.numeric))

  res <- ptw(reference_num, samples_num, warp.type = "individual")

  warped_df    <- as.data.frame(res$warped.sample)
  reference_df <- as.data.frame(res$reference)
  colnames(warped_df) <- colnames(reference_df)

  bind_rows(
    bind_cols(reference_meta, reference_df),
    bind_cols(samples_meta,   warped_df)
  )
}

# ── run ───────────────────────────────────────────────────────────────────────

allWide_un <- read_parquet(here("data", "prepared", "combinedWide.parquet"))

message("Applying baseline correction…")
allWide <- correct_baseline(allWide_un, meta_in = TRUE)

species_df <- split_by_species(allWide)

dir.create(here("data", "processed"), showWarnings = FALSE)

species_names <- c("green", "white", "early", "wyoming")
out_files <- c(
  green   = here("data", "processed", "greenPTW.parquet"),
  white   = here("data", "processed", "whitePTW.parquet"),
  early   = here("data", "processed", "earlyPTW.parquet"),
  wyoming = here("data", "processed", "wyomingPTW.parquet")
)

for (sp in species_names) {
  message("Aligning ", sp, "…")
  res <- ptw_function(species_df[[sp]])
  write_parquet(res, out_files[[sp]])
  message("  Saved → ", out_files[[sp]])
}

message("Done. All PTW parquet files written to data/processed/")
