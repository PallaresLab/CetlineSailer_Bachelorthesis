# ============================================================================
# 02_F2_old_genotype_comparison.R
# F2 population, first round (rec_060-065): merges AmpSeq genotype calls
# (0 = wildtype, 1 = heterozygous, 2 = tim0/tim0) with eclosion times and
# runs the Mendelian check, ZT histograms, and pairwise significance tests
# reported for this dataset.
#
# Produces:
#   eclosion numbers/rates + Mendelian ratio per genotype
#   ZT histogram, all three genotypes
#   circular descriptive stats per genotype
#   pairwise Watson's two-sample tests (circular)
#
# Input files expected in data/F2_old/:
#   genotypes_filtered.txt        (AmpSeq calls, tab-separated, no header)
#   all_recordings_combined.csv   (eclosion times, rec_060-065)
# ============================================================================

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "utils_circ_stats.R"))

# ---- Step 1: Load AmpSeq genotype calls -------------------------------------
genotypes <- read.table(
  file.path(dir_data_F2old, "genotypes_filtered.txt"),
  header = FALSE,
  col.names = c("plate_well_raw", "gene", "allele", "genotype", "depth_in", "depth_out")
) %>%
  mutate(plate_well = gsub("^p(\\d+)([A-H]\\d+)$", "\\1_\\2", plate_well_raw))

# ---- Step 2: Load eclosion time data -----------------------------------------
eclosion <- read.csv(
  file.path(dir_data_F2old, "all_recordings_combined.csv"),
  na.strings = c("", "NA")
) %>%
  mutate(ecl_time = parse_ecl_time(ecl_time))

# ---- Step 3: Merge genotype with eclosion data (left join) ------------------
merged <- merge(eclosion, genotypes %>% select(plate_well, genotype),
                 by = "plate_well", all.x = TRUE)

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

write.csv(merged_clean, file.path(dir_out_tab, "F2_old_eclosion_genotypes_merged.csv"), row.names = FALSE)

# ---- Step 4: Prepare ZT variables --------------------------------------------
plot_data_F2 <- merged_clean %>%
  filter(!is.na(ecl_time)) %>%
  mutate(genotype = factor(genotype_label, levels = all_genotypes)) %>%
  add_ZT_columns(lights_on_h = lights_on_winter)

# Save the ZT-enriched dataset so 04_resampling_P0_vs_F2.R can reuse it
# directly instead of re-deriving F2_old from scratch.
write.csv(plot_data_F2, file.path(dir_out_tab, "F2_old_merged_with_ZT.csv"), row.names = FALSE)

n_wt  <- sum(plot_data_F2$genotype == "Wildtype")
n_het <- sum(plot_data_F2$genotype == "Heterozygous")
n_mut <- sum(plot_data_F2$genotype == "tim0")

# ---- ZT histogram, all three genotypes ----------------------------
fig_16 <- plot_data_F2 %>%
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
  facet_wrap(~ genotype, ncol = 3,
             labeller = as_labeller(c(
               "Wildtype"     = paste0("Wildtype (n == ", n_wt, ")"),
               "Heterozygous" = paste0("Heterozygous (n == ", n_het, ")"),
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
fig_16
ggsave(file.path(dir_out_fig, "F2_old_eclosion_timing_all_genotypes.pdf"), fig_16, width = 14, height = 5)

# ---- circular descriptive stats per genotype ----------------------
table_16 <- plot_data_F2 %>%
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
print(table_16)
write.csv(table_16, file.path(dir_out_tab, "F2_old_circular_stats.csv"), row.names = FALSE)

# ---- pairwise Watson's two-sample tests (circular) ----------------
wt_circ  <- circular(plot_data_F2$ZT_plot[plot_data_F2$genotype == "Wildtype"],     units = "hours", template = "clock24")
het_circ <- circular(plot_data_F2$ZT_plot[plot_data_F2$genotype == "Heterozygous"], units = "hours", template = "clock24")
mut_circ <- circular(plot_data_F2$ZT_plot[plot_data_F2$genotype == "tim0"],         units = "hours", template = "clock24")

watson_wt_mut <- watson.two.test(wt_circ, mut_circ)
watson_wt_het <- watson.two.test(wt_circ, het_circ)
watson_het_mut <- watson.two.test(het_circ, mut_circ)

table_17 <- data.frame(
  comparison = c("WT vs. mutant", "WT vs. heterozygous", "Heterozygous vs. mutant"),
  watson_U2  = c(watson_wt_mut$statistic, watson_wt_het$statistic, watson_het_mut$statistic),
  p_value    = c(watson_pvalue_text(watson_wt_mut), watson_pvalue_text(watson_wt_het), watson_pvalue_text(watson_het_mut))
)
print(table_17)
write.csv(table_17, file.path(dir_out_tab, "F2_old_watson_tests.csv"), row.names = FALSE)

# ---- eclosion numbers/rates + Mendelian ratio ---------------------
eclosion_rate_F2 <- merged_clean %>%
  group_by(genotype_label) %>%
  summarise(
    n_genotyped = n(),
    n_eclosed   = sum(!is.na(ecl_time)),
    eclosion_rate_pct = round(n_eclosed / n_genotyped * 100, 1),
    .groups = "drop"
  ) %>%
  mutate(observed_pct = round(n_genotyped / sum(n_genotyped) * 100, 1))

observed <- setNames(rep(0, 3), all_genotypes)
observed[eclosion_rate_F2$genotype_label] <- eclosion_rate_F2$n_genotyped
chisq_result <- chisq.test(observed, p = c(0.25, 0.5, 0.25))
print(chisq_result)

table_14 <- eclosion_rate_F2 %>% mutate(genotype_label = factor(genotype_label, levels = all_genotypes))
print(table_14)
write.csv(table_14, file.path(dir_out_tab, "F2_old_eclosion_mendelian.csv"), row.names = FALSE)


# ---- Robustness check: is the genotype effect on eclosion rate consistent --
# across recordings, or driven by one or two outlier recordings? Simple
# descriptive check (not a formal test) - percent eclosed per genotype,
# broken down per recording.
eclosion_rate_per_rec <- merged %>%
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
print(eclosion_rate_per_rec)
write.csv(eclosion_rate_per_rec, file.path(dir_out_tab, "eclosion_rate_per_rec_F2_old.csv"), row.names = FALSE)
