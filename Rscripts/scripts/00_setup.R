# ============================================================================
# 00_setup.R
# Load once at the start of a session (or via source() at the top of any
# other script). Defines portable paths (via `here`), shared constants, and
# libraries used across the analysis scripts.
#
# IMPORTANT: run this only after opening the project via the .Rproj file, or
# after setting your working directory to the project root. `here::here()`
# then always resolves relative to that root, no matter which machine or
# which subfolder a script sits in - no more "H:/eclosion_monitor/..." paths.
# ============================================================================

# install.packages(c("dplyr","tidyr","ggplot2","readr","lubridate","circular","here"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(lubridate)
library(circular)
library(here)

# ---- Folder structure expected by the scripts in this project --------------
# data/P0/...       raw + merged P0 (parental generation) recordings
# data/F2_old/...   F2 rec_060-065, AmpSeq genotype calls
# data/F2_new/...   F2 rec_083-089, AmpSeq genotype calls
# output/figures/   all ggsave() output lands here
# output/tables/    all write.csv() output lands here

dir_data_P0     <- here("data", "P0")
dir_data_F2old  <- here("data", "F2_old")
dir_data_F2new  <- here("data", "F2_new")
dir_out_fig     <- here("output", "figures")
dir_out_tab     <- here("output", "tables")

dir.create(dir_out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_out_tab, recursive = TRUE, showWarnings = FALSE)

# ---- Shared constants (color scheme, ZT settings) --------------------------
col_day       <- "#E8C456"
col_night     <- "#4A5568"
bin_size      <- 1                 # 1-hour bins for all ZT histograms

lights_on_winter <- 7 + 20/60      # 07:20, used for P0_without_air F2_old (winter time)
lights_on_summer <- 8 + 20/60      # 08:20, used for P0 and F2_new (summer time)

all_genotypes <- c("Wildtype", "Heterozygous", "tim0")

# ---- Robust datetime parser --------------------------------------------------
# Some raw recording CSVs (confirmed for at least rec_088) use an alternate
# date format, M/D/YYYY H:MM instead of YYYY-MM-DD HH:MM:SS - likely from a
# row being deleted/edited in the CSV at some point, which changed the cell
# format for the remaining rows. A single-format as.POSIXct() call silently
# turns those rows into NA, which undercounts eclosed pupae without any
# warning. Use this everywhere a raw ecl_time column is parsed.
parse_ecl_time <- function(x) {
  parsed <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  still_na <- is.na(parsed) & !is.na(x) & x != ""
  parsed[still_na] <- as.POSIXct(x[still_na], format = "%m/%d/%Y %H:%M", tz = "UTC")
  parsed
}


# ---- Shared helper: convert ecl_time -> ZT / ZT_plot / daytime -------------
# `lights_on_h` lets you pass the correct offset for the dataset you're
# working with (lights_on_winter for P0/F2_old, lights_on_summer for F2_new).
add_ZT_columns <- function(data, ecl_time_col = "ecl_time", lights_on_h) {
  data %>%
    mutate(
      .ecl_time_local = force_tz(.data[[ecl_time_col]], "Europe/Berlin"),
      hour_of_day     = hour(.ecl_time_local) + minute(.ecl_time_local) / 60,
      ZT              = (hour_of_day - lights_on_h) %% 24,
      # Wraps into [-12, 12) consistently. (ZT > 12) ? ZT-24 : ZT would leave
      # the exact boundary ZT == 12 stranded in its own always-empty bin
      # instead of joining the -12 bin where ZT values just above 12 land -
      # this formula avoids that seam artifact.
      ZT_plot         = ((ZT + 12) %% 24) - 12,
      ZT_bin          = floor(ZT_plot / bin_size) * bin_size,
      daytime         = ifelse(ZT_bin >= 0 & ZT_bin < 12, "Day", "Night")
    ) %>%
    select(-.ecl_time_local)
}
