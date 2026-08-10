# ============================================================================
# detect_eclosion_template.R
# Detects eclosion time of Drosophila from rolling-median pixel-intensity
# traces (light/dark transition per well), for one recording at a time.
#
# Here the pipeline is a single parameterized template;
# you only change `rec_id` and the two folder paths below to reprocess a
# different recording, and manual corrections live in their own CSV under
# data/F2_new/manual_corrections/ instead of being hardcoded.
#
# ============================================================================

source(here::here("scripts", "00_setup.R"))
library(zoo)

# ---- Parameters: change these two lines to reprocess a different rec -------
rec_id  <- "rec_083"
roi_dir <- file.path(dir_data_F2new, rec_id, "ref")          # roi_coords.csv, well_coords.csv, WT_P0.txt, F2.txt
raw_dir <- file.path(dir_data_F2new, rec_id, "for_analysis") # results.csv (raw pixel intensity export)

corrections_path <- file.path(dir_data_F2new, "manual_corrections", paste0(rec_id, "_corrections.csv"))

# ---- Step 1: Load ROI coordinates and well layout ----------------------------
roi_coords <- read.csv(file.path(roi_dir, "roi_coords.csv")) %>%
  mutate(id = gsub(".*:", "", Label)) %>%
  select(id, X, Y)

well_coords <- read.csv(file.path(roi_dir, "well_coords.csv"))
wells <- c(outer(LETTERS[1:8], 1:12, paste0))
well_coords$plate <- rep(1:2, each = 96)

well_coords <- well_coords %>%
  mutate(
    well = rep(wells, times = 2),
    row  = substr(well, 1, 1),
    genotype = case_when(
      plate == 1                                    ~ "F2",
      plate == 2 & row %in% c("A", "B", "C")         ~ "WT_P0",
      plate == 2 & row %in% c("D", "E", "F", "G", "H") ~ "F2",
      TRUE ~ NA_character_
    ),
    plate_well = paste0(plate, "_", well)
  ) %>%
  select(well, X, Y, plate, genotype, plate_well)

# ---- Step 2: Match each ROI to its closest well ------------------------------
id_well <- expand.grid(
  id         = as.character(roi_coords$id),
  plate_well = as.character(well_coords$plate_well)
) %>%
  left_join(roi_coords %>% select(id, x_roi = X, y_roi = Y), by = "id") %>%
  left_join(well_coords %>% select(plate_well, x_well = X, y_well = Y, plate_well_num = plate, genotype),
            by = "plate_well") %>%
  mutate(D = sqrt((x_roi - x_well)^2 + (y_roi - y_well)^2)) %>%
  group_by(id) %>%
  filter(D == min(D)) %>%
  ungroup() %>%
  mutate(id = as.character(id), plate_well = as.character(plate_well))

# ---- Step 3: Load pixel intensity data and attach well/genotype -------------
d <- read.csv(file.path(raw_dir, "results.csv"))

datetime_raw <- sub(".*([0-9]{4}_[0-9]{2}_[0-9]{2}_[0-9]{2}_[0-9]{2}).*", "\\1", d$Label)
d$datetime <- as.POSIXct(datetime_raw, format = "%Y_%m_%d_%H_%M", tz = "Europe/Berlin")
d$id       <- as.character(sub(".*:", "", d$Label))

d <- merge(d, id_well, by = "id") %>%
  arrange(datetime, id) %>%
  mutate(plate = as.numeric(sub("_.*", "", plate_well)))

# ---- Step 4: Rolling median + eclosion detection helpers --------------------
comp_roll_med <- function(data, k = 5, width = 24 * 6) {
  data %>%
    group_by(id) %>%
    mutate(
      roll_median  = rollmedian(Median, k = k, fill = NA, align = "center"),
      med_before   = rollapply(roll_median, width = width, FUN = mean, align = "right", fill = NA, partial = FALSE),
      med_after    = rollapply(roll_median, width = width, FUN = mean, align = "left",  fill = NA, partial = FALSE),
      diff_medians = med_after - med_before
    ) %>%
    filter(!is.na(roll_median)) %>%
    ungroup()
}

get_ecl_time <- function(data, width = 24 * 6) {
  t_switch  <- data[which.min(data$diff_medians), ]$datetime
  data_pre  <- subset(data, datetime <  t_switch)
  data_post <- subset(data, datetime >= t_switch)
  data_pre  <- data_pre[max(1, nrow(data_pre) - width):nrow(data_pre), ]
  data_post <- data_post[1:min(nrow(data_post), width), ]
  if (nrow(data_pre) < width | nrow(data_post) < width) return(NA)
  if (median(data_pre$roll_median) > median(data_post$roll_median) &&
      wilcox.test(data_pre$roll_median, data_post$roll_median)$p.value < 0.01) {
    return(t_switch)
  } else {
    return(NA)
  }
}

# ---- Step 5: Detect eclosion times --------------------------------------------
d   <- comp_roll_med(d, k = 5, width = 24 * 6)
res <- d %>%
  group_by(plate_well) %>%
  group_modify(~ tibble(ecl_time = get_ecl_time(.x, width = 24 * 6))) %>%
  ungroup() %>%
  mutate(plate_well = as.character(plate_well))

cat("Automatically detected:", sum(!is.na(res$ecl_time)), "/", nrow(res), "\n")

# ---- Step 6: Visual QC (inspect before applying manual corrections) ---------
merge(d, res) %>%
  ggplot(aes(datetime, roll_median)) +
  geom_line() +
  geom_vline(aes(xintercept = ecl_time), linetype = 2, col = "red") +
  facet_wrap(~ plate_well) +
  theme_bw() +
  ggtitle(rec_id)

# ---- Step 7: Apply manual corrections from the CSV ----------------------------
res <- merge(res, id_well %>% select(plate_well, genotype), by = "plate_well")

if (file.exists(corrections_path)) {
  corrections <- read.csv(corrections_path, stringsAsFactors = FALSE)

  reset_wells   <- corrections$plate_well[corrections$action == "reset_to_NA"]
  keep_wells    <- corrections$plate_well[corrections$action == "keep_detected"]
  manual_wells  <- corrections$plate_well[corrections$action == "manual_override"]
  manual_times  <- corrections$ecl_time[corrections$action == "manual_override"]

  # rec_090 used the reverse logic in the original script (reset everything
  # except a short verified list) - if `keep_detected` rows are present,
  # treat this recording that way; otherwise reset only the listed wells.
  if (length(keep_wells) > 0) {
    res$ecl_time <- ifelse(res$plate_well %in% keep_wells, res$ecl_time, NA)
  } else {
    res$ecl_time[res$plate_well %in% reset_wells] <- NA
  }
  res$ecl_time <- as.POSIXct(res$ecl_time, tz = "Europe/Berlin", origin = "1970-01-01")

  if (length(manual_wells) > 0) {
    manual_ecl <- data.frame(plate_well = manual_wells,
                              ecl_time   = as.POSIXct(manual_times, tz = "Europe/Berlin"))
    res$ecl_time[match(manual_ecl$plate_well, res$plate_well)] <- manual_ecl$ecl_time
  }
} else {
  cat("No corrections file found for", rec_id, "- keeping automatic detections as-is.\n")
}

# ---- Step 8: Save the per-recording result -------------------------------------
res_final <- res %>%
  mutate(rec = rec_id) %>%
  select(rec, plate_well, id, ecl_time, genotype)

write.csv(res_final,
          file.path(dir_data_F2new, paste0("eclosiontime_", rec_id, "_results.csv")),
          quote = FALSE, row.names = FALSE)

cat("Saved:", file.path(dir_data_F2new, paste0("eclosiontime_", rec_id, "_results.csv")), "\n")
