# ============================================================================
# 03_F2_new_genotype_comparison.R
# F2 population, second round / new genomic background (rec_083-089):
# merges AmpSeq genotype calls with eclosion times, runs the Mendelian
# check, ZT histograms, and circular statistics - INCLUDING the pairwise
# Watson's two-sample tests that were still missing for this dataset.
#
# Produces:
#   eclosion numbers/rates + Mendelian ratio per genotype
#   ZT histogram, wildtype vs. mutant
#   circular descriptive stats, wildtype + mutant
#   (NEW)- pairwise Watson's two-sample tests (all 3 genotypes),
#                    the comparison you were missing for F2_new
#
# Input files expected in data/F2_new/:
#   genotypes_filtered_new.txt          (AmpSeq calls, depth-filtered, tab-sep)
#   eclosiontime_rec0<83..89>_results.csv   (one file per recording)
# ============================================================================

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "utils_circ_stats.R"))

rec_ids <- 83:89

# ---- Step 1: Load filtered AmpSeq genotype calls -----------------------------
genotypes_filtered <- read.table(
  file.path(dir_data_F2new, "genotypes_filtered_new.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
) %>%
  mutate(plate_well = gsub("^0+(\\d)", "\\1", plate_well))

# ---- Step 2: Load & combine eclosion time data from individual rec files ----
# (parse_ecl_time() is defined in 00_setup.R - handles rec_088's alternate
# date format automatically)
eclosion <- bind_rows(lapply(rec_ids, function(recid) {
  read.csv(file.path(dir_data_F2new, paste0("eclosiontime_rec0", recid, "_results.csv")),
           na.strings = c("", "NA")) %>%
    mutate(
      rec      = paste0("rec_0", recid),
      ecl_time = parse_ecl_time(ecl_time)
    )
}))

# Keep the F2 wells only (WT_P0 control wells added per recording are handled
# separately, see 04_resampling_P0_vs_F2.R)
eclosion_F2 <- eclosion %>%
  filter(genotype == "F2") %>%
  select(-genotype)

# Each recording covers 2 plates; renumber plate_well so plate numbers are
# unique across all 7 recordings (rec_083 -> plates 1-2, rec_084 -> 3-4, ...)
plate_offset_map <- data.frame(
  rec    = paste0("rec_0", rec_ids),
  offset = (seq_along(rec_ids) - 1) * 2
)

eclosion_F2 <- eclosion_F2 %>%
  left_join(plate_offset_map, by = "rec") %>%
  mutate(
    plate_num_old = as.numeric(sub("^(\\d+)_.*$", "\\1", plate_well)),
    well_part     = sub("^\\d+_(.*)$", "\\1", plate_well),
    plate_num_new = plate_num_old + offset,
    plate_well    = paste0(plate_num_new, "_", well_part)
  ) %>%
  select(-plate_num_old, -well_part, -plate_num_new, -offset)

# ---- Step 3: Merge genotype with eclosion data (left join) ------------------
merged <- merge(eclosion_F2, genotypes_filtered, by = "plate_well", all.x = TRUE)

merged_clean <- merged %>%
  filter(!is.na(genotype)) %>%
  mutate(
    genotype_label = case_when(
      genotype == 0 ~ "Wildtype",
      genotype == 1 ~ "Heterozygous",
      genotype == 2 ~ "tim0",
      TRUE ~ NA_character_
    )
  ) %>%
  select(rec, plate_well, id, ecl_time, genotype, genotype_label)

write.csv(merged_clean, file.path(dir_out_tab, "F2_new_eclosion_genotypes_merged.csv"), row.names = FALSE)

# ---- Step 4: Prepare ZT variables (summer time!) -----------------------------
plot_data_F2new <- merged_clean %>%
  filter(!is.na(ecl_time)) %>%
  mutate(genotype = factor(genotype_label, levels = all_genotypes)) %>%
  add_ZT_columns(lights_on_h = lights_on_summer)

# Save the ZT-enriched dataset so 04_resampling_P0_vs_F2.R can reuse it
# directly instead of re-deriving F2_new from scratch.
write.csv(plot_data_F2new, file.path(dir_out_tab, "F2_new_merged_with_ZT.csv"), row.names = FALSE)

n_wt  <- sum(plot_data_F2new$genotype == "Wildtype")
n_het <- sum(plot_data_F2new$genotype == "Heterozygous")
n_mut <- sum(plot_data_F2new$genotype == "tim0")

# ---- eclosion numbers/rates + Mendelian ratio ---------------------
eclosion_rate_F2new <- merged_clean %>%
  group_by(genotype_label) %>%
  summarise(
    n_genotyped = n(),
    n_eclosed   = sum(!is.na(ecl_time)),
    eclosion_rate_pct = round(n_eclosed / n_genotyped * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(observed_pct = round(n_genotyped / sum(n_genotyped) * 100, 1))

observed <- setNames(rep(0, 3), all_genotypes)
observed[eclosion_rate_F2new$genotype_label] <- eclosion_rate_F2new$n_genotyped
chisq_result <- chisq.test(observed, p = c(0.25, 0.5, 0.25))
print(chisq_result)

table_18 <- eclosion_rate_F2new %>% mutate(genotype_label = factor(genotype_label, levels = all_genotypes))
print(table_18)
write.csv(table_18, file.path(dir_out_tab, "F2_new_eclosion_mendelian.csv"), row.names = FALSE)


# ---- Robustness check: is the genotype effect on eclosion rate consistent --
# across recordings, or driven by one or two outlier recordings? Simple
# descriptive check (not a formal test) - percent eclosed per genotype,
# broken down per recording.
eclosion_rate_per_rec_new <- merged %>%
  mutate(
    genotype_label = case_when(
      genotype == 0 ~ "Wildtype",
      genotype == 1 ~ "Heterozygous",
      genotype == 2 ~ "tim0",
      TRUE ~ NA_character_
    ),
    eclosed = !is.na(ecl_time)
  ) %>%
  filter(!is.na(genotype_label)) %>%
  group_by(rec, genotype_label) %>%
  summarise(pct_eclosed = round(mean(eclosed) * 100, 1), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = genotype_label, values_from = c(pct_eclosed, n))
print(eclosion_rate_per_rec_new)
write.csv(eclosion_rate_per_rec_new, file.path(dir_out_tab, "eclosion_rate_per_rec_F2_new.csv"), row.names = FALSE)

# ---- ZT histogram, wildtype vs. mutant only -----------------------
# (per Results 4.3.4: analysis for the new background focuses on WT vs.
# mutant, heterozygotes excluded from this particular comparison)
plot_data_wt_mut <- plot_data_F2new %>%
  filter(genotype %in% c("Wildtype", "tim0")) %>%
  mutate(genotype = droplevels(genotype))

fig_20 <- plot_data_wt_mut %>%
  count(genotype, ZT_bin, daytime) %>%
  group_by(genotype) %>%
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
  facet_wrap(~ genotype, ncol = 2,
             labeller = as_labeller(c(
               "Wildtype" = paste0("Wildtype (n == ", n_wt, ")"),
               "tim0"     = paste0("italic(tim)^{'01'}~(n == ", n_mut, ")")
             ), label_parsed)) +
  labs(fill = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 11, face = "bold"),
    legend.position   = "bottom",
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )
fig_20
ggsave(file.path(dir_out_fig, "F2_new_eclosion_timing_WT_vs_mutant.pdf"), fig_20, width = 10, height = 5)

# ---- circular descriptive stats, wildtype + mutant ----------------
table_19 <- plot_data_wt_mut %>%
  group_by(genotype) %>%
  summarise(
    n            = n(),
    circ_mean_ZT = as.numeric(mean(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)),
    circ_sd_ZT   = as.numeric(sd(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)) * 24 / (2 * pi),
    rayleigh_p   = rayleigh.test(circular(ZT_plot, units = "hours", template = "clock24"))$p.value,
    pct_light_on = round(sum(daytime == "Day",   na.rm = TRUE) / n() * 100, 1),
    pct_dark     = round(sum(daytime == "Night", na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  )
print(table_19)
write.csv(table_19, file.path(dir_out_tab, "F2_new_circular_stats_WT_vs_mutant.csv"), row.names = FALSE)

# ==============================================================================
# NEW: pairwise Watson's two-sample tests (circular), F2_new
# This is the comparison that was still missing across genotypes for the
# rec_083-089 dataset. Mirrors the approach used for F2_old
# (02_F2_old_genotype_comparison.R), so results are directly comparable.
# ==============================================================================
wt_circ_new  <- circular(plot_data_F2new$ZT_plot[plot_data_F2new$genotype == "Wildtype"],     units = "hours", template = "clock24")
het_circ_new <- circular(plot_data_F2new$ZT_plot[plot_data_F2new$genotype == "Heterozygous"], units = "hours", template = "clock24")
mut_circ_new <- circular(plot_data_F2new$ZT_plot[plot_data_F2new$genotype == "tim0"],         units = "hours", template = "clock24")

cat("\nWatson's two-sample test (F2_new): WT vs Mutant\n");  watson_wt_mut_new  <- watson.two.test(wt_circ_new, mut_circ_new);  print(watson_wt_mut_new)
cat("\nWatson's two-sample test (F2_new): WT vs Het\n");     watson_wt_het_new  <- watson.two.test(wt_circ_new, het_circ_new);  print(watson_wt_het_new)
cat("\nWatson's two-sample test (F2_new): Het vs Mutant\n"); watson_het_mut_new <- watson.two.test(het_circ_new, mut_circ_new); print(watson_het_mut_new)

table_19b <- data.frame(
  comparison = c("WT vs. mutant", "WT vs. heterozygous", "Heterozygous vs. mutant"),
  watson_U2  = c(watson_wt_mut_new$statistic, watson_wt_het_new$statistic, watson_het_mut_new$statistic),
  p_value    = c(watson_pvalue_text(watson_wt_mut_new), watson_pvalue_text(watson_wt_het_new), watson_pvalue_text(watson_het_mut_new))
)
print(table_19b)
write.csv(table_19b, file.path(dir_out_tab, "F2_new_watson_tests.csv"), row.names = FALSE)
# Note: watson.two.test() only computes a critical-value band, not an exact
# numeric p-value - watson_pvalue_text() (in utils_circ_stats.R) extracts
# that band from the test's print() output automatically.

# ==============================================================================
# Dominance/recessiveness check: is tim01 recessive or dominant?
# Read this off table_19b:
#   WT vs. Het n.s., Het vs. Mutant sig.  -> heterozygote resembles WT      -> recessive
#   WT vs. Het sig.,  Het vs. Mutant n.s. -> heterozygote resembles mutant  -> dominant
#   both significant (Het falls "in between")                              -> incomplete dominance / additive
#
# and give the visual + descriptive-stats counterpart
# to that call: all three genotypes side by side (unlike , which
# only shows WT vs. mutant), so you can actually see where the heterozygote
# sits relative to the two homozygous classes.
# ==============================================================================

n_het_new <- sum(plot_data_F2new$genotype == "Heterozygous")

fig_20b <- plot_data_F2new %>%
  count(genotype, ZT_bin, daytime) %>%
  group_by(genotype) %>%
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
  scale_y_continuous(name = "Number of eclosions", limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ genotype, ncol = 3,
             labeller = as_labeller(c(
               "Wildtype"     = paste0("Wildtype (n == ", n_wt, ")"),
               "Heterozygous" = paste0("Heterozygous (n == ", n_het_new, ")"),
               "tim0"         = paste0("italic(tim)^{'01'}~(n == ", n_mut, ")")
             ), label_parsed)) +
  labs(fill = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 11, face = "bold"),
    legend.position   = "bottom",
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )
fig_20b
ggsave(file.path(dir_out_fig, "F2_new_eclosion_timing_all_genotypes.pdf"), fig_20b, width = 14, height = 5)

table_19c <- plot_data_F2new %>%
  group_by(genotype) %>%
  summarise(
    n            = n(),
    circ_mean_ZT = as.numeric(mean(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)),
    circ_sd_ZT   = as.numeric(sd(circular(ZT_plot, units = "hours", template = "clock24"), na.rm = TRUE)) * 24 / (2 * pi),
    rayleigh_p   = rayleigh.test(circular(ZT_plot, units = "hours", template = "clock24"))$p.value,
    pct_light_on = round(sum(daytime == "Day",   na.rm = TRUE) / n() * 100, 1),
    pct_dark     = round(sum(daytime == "Night", na.rm = TRUE) / n() * 100, 1),
    .groups = "drop"
  )
print(table_19c)
write.csv(table_19c, file.path(dir_out_tab, "F2_new_circular_stats_all_genotypes.csv"), row.names = FALSE)
