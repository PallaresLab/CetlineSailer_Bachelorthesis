# ============================================================================
# 01_P0_positive_control.R
# P0 (parental generation): wildtype vs. tim01 mutant eclosion timing.
# Serves as the positive control confirming the established tim01 arrhythmic
# phenotype before moving on to the F2 recombinant population.
#
# Produces:
#   - eclosion rate per genotype (n eclosed / n total pupae)
#   - ZT histogram, wildtype vs. mutant
#   - circular descriptive stats (mean, SD, Rayleigh test) + % light/dark per genotype
#   - Watson's two-sample test (circular): wildtype vs. mutant
#
# Built directly from the three raw per-recording files (rec_067-069) in
# data/P0/, rather than from all_res_new_corrected.csv - that merged file
# turned out to have at least one row with a truncated, date-only ecl_time
# (no time-of-day), so it can't be trusted as the sole source. rec_069's
# timestamps also need the same +1h correction used in
# 06_appendix_P0_per_recording.R (its raw times ran 1h behind the other two
# recordings, e.g. a real 12:00 eclosion is stored as 11:00).
#
# Input files expected in data/P0/:
#   eclosiontime rec_067_results.csv   (wildtype, both plates)
#   eclosiontime rec_068_results.csv   (mutant,   both plates)
#   eclosiontime rec_069_results.csv   (plate 1 = wildtype, plate 2 = mutant)
# ============================================================================

source(here::here("scripts", "00_setup.R"))

# ---- Step 1: Load data from the three raw recordings, compute ZT -----------
files_P0 <- c(
  file.path(dir_data_P0, "eclosiontime rec_067_results.csv"),
  file.path(dir_data_P0, "eclosiontime rec_068_results.csv"),
  file.path(dir_data_P0, "eclosiontime rec_069_results.csv")
)

genotype_map_P0 <- data.frame(
  rec       = c("rec_067", "rec_067", "rec_068", "rec_068", "rec_069", "rec_069"),
  plate_num = c(1, 2, 1, 2, 1, 2),
  genotype  = c("wildtype", "wildtype", "tim0", "tim0", "wildtype", "tim0"),
  stringsAsFactors = FALSE
)

p0_raw <- bind_rows(lapply(files_P0, function(f) {
  df <- read.csv(f, na.strings = c("", "NA"))
  df$ecl_time <- parse_ecl_time(df$ecl_time)
  df$rec <- regmatches(f, regexpr("rec_\\d+", f))
  # raw timestamps in these recordings run 1h behind the real time
  # (a real 12:00 eclosion is stored as 11:00) - correct all of them.
  df$ecl_time <- df$ecl_time + hours(1)
  df
})) %>%
  mutate(plate_num = as.numeric(sub("_.*", "", plate_well))) %>%
  left_join(genotype_map_P0, by = c("rec", "plate_num"))

p0_data <- p0_raw %>%
  filter(!is.na(ecl_time), genotype %in% c("wildtype", "tim0")) %>%
  mutate(
    genotype_label = ifelse(genotype == "wildtype", "Wildtype", "Mutant"),
    genotype_label = factor(genotype_label, levels = c("Wildtype", "Mutant"))
  ) %>%
  add_ZT_columns(lights_on_h = lights_on_summer)

# Save the ZT-enriched dataset so 04_resampling_P0_vs_F2.R can reuse it
# directly, instead of re-deriving P0 from scratch with separate (and
# possibly out-of-sync) loading logic.
write.csv(p0_data, file.path(dir_out_tab, "P0_merged_with_ZT.csv"), row.names = FALSE)

# ---- eclosion rate per genotype -----------------------------------
# n_pupae_total must be supplied from the plate design (not derivable from
# eclosion times alone, since it also counts pupae that never eclosed).
n_pupae_total <- c(Wildtype = 281, Mutant = 283)

table_12 <- p0_data %>%
  count(genotype_label, name = "total_eclosed") %>%
  mutate(
    total_pupae = n_pupae_total[as.character(genotype_label)],
    percent     = round(total_eclosed / total_pupae * 100, 1)
  )
print(table_12)
write.csv(table_12, file.path(dir_out_tab, "P0_eclosion_rate.csv"), row.names = FALSE)

# ---- ZT histogram, wildtype vs. mutant ---------------------------
n_wt  <- sum(p0_data$genotype_label == "Wildtype")
n_mut <- sum(p0_data$genotype_label == "Mutant")

fig_11 <- p0_data %>%
  count(genotype_label, ZT_bin, daytime) %>%
  group_by(genotype_label) %>%
  complete(ZT_bin, fill = list(n = 0, daytime = NA)) %>%
  ungroup() %>%
  mutate(daytime = ifelse(ZT_bin >= 0 & ZT_bin < 12, "Day", "Night")) %>%
  ggplot(aes(x = ZT_bin + bin_size / 2, y = n, fill = daytime)) +
  geom_vline(xintercept = 0,  linetype = "dashed", color = "red",       linewidth = 0.7) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "steelblue", linewidth = 0.7) +
  annotate("text", x = -0.5, y = Inf, label = "Lights on",  vjust = 2, color = "red",       size = 3, hjust = 1) +
  annotate("text", x = 11.7, y = Inf, label = "Lights off", vjust = 2, color = "steelblue", size = 3, hjust = 1) +
  geom_col(width = bin_size, position = "identity", color = "white", linewidth = 0.2, alpha = 0.85) +
  scale_fill_manual(values = c("Day" = col_day, "Night" = col_night), na.value = "grey80") +
  scale_x_continuous(name = "Zeitgeber Time (ZT, h)", breaks = seq(-12, 12, by = 12),
                      limits = c(-12, 13.5), labels = function(x) paste0("ZT ", x)) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(name = "Number of eclosions", limits = c(0, NA),
                      expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ genotype_label, ncol = 2,
             labeller = as_labeller(c(
               "Wildtype" = paste0("Wildtype (n == ", n_wt, ")"),
               "Mutant"   = paste0("italic(tim)^{'01'}~(n == ", n_mut, ")")
             ), label_parsed)) +
  labs(fill = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 11, face = "bold"),
    legend.position   = "bottom",
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )
fig_11
ggsave(file.path(dir_out_fig, "P0_eclosion_timing.pdf"), fig_11, width = 10, height = 5)

# ---- circular descriptive stats + % light/dark -------------------
table_13 <- p0_data %>%
  group_by(genotype_label) %>%
  summarise(
    n            = n(),
    circ_mean_ZT = as.numeric(mean(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)),
    circ_sd_ZT   = as.numeric(sd(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)) * 24 / (2 * pi),
    rayleigh_p   = rayleigh.test(circular(ZT_plot, units = "hours", template = "clock24"))$p.value,
    pct_light_on = round(sum(daytime == "Day",   na.rm = TRUE) / n() * 100, 1),
    pct_dark     = round(sum(daytime == "Night", na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  )
print(table_13)
write.csv(table_13, file.path(dir_out_tab, "P0_circular_stats.csv"), row.names = FALSE)

# ---- Watson's two-sample test: wildtype vs. mutant --------------------------
wt_circ  <- circular(p0_data$ZT_plot[p0_data$genotype_label == "Wildtype"], units = "hours", template = "clock24")
mut_circ <- circular(p0_data$ZT_plot[p0_data$genotype_label == "Mutant"],   units = "hours", template = "clock24")

watson_p0 <- watson.two.test(wt_circ, mut_circ)
print(watson_p0)
