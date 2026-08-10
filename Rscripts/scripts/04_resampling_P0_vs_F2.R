# ============================================================================
# 04_resampling_P0_vs_F2.R
# Circular resampling comparisons between the P0 wildtype reference and (a)
# the WT_P0 control flies added into each F2_new recording (technical
# control), and (b) the F2 wildtype/mutant data itself, pooled per
# generation (background effect). Uses resample_2() / summarise_circ_stats()
# / plot_resample_result() from utils_circ_stats.R.
#
# Produces:
#   resampled P0 wildtype vs. WT_P0 control, per recording
#               (rec_083-089, including rec_089 as shown in the thesis).
#               This is the only per-recording comparison that cleanly
#               isolates technical/recording variability, because WT_P0
#               controls are genetically identical to the p0 reference -
#               any difference found here can't be genetic.
#   pooled P0 vs. F2 overlap (wildtype + mutant, F2_old and
#               F2_new each get their own figure): answers whether the
#               genomic background itself shifts the phenotype. F2 wildtype/
#               mutant are NOT genetically identical to p0 (recombinant
#               background), so this question is deliberately answered with
#               the pooled comparison rather than per-recording, which
#               couldn't distinguish technical from biological variation
#               there anyway.
#
# IMPORTANT: run 01_P0_positive_control.R, 02_F2_old_genotype_comparison.R,
# and 03_F2_new_genotype_comparison.R BEFORE this script. Rather than
# re-deriving P0/F2_old/F2_new from raw files with its own separate loading
# logic (which had drifted out of sync with the fixes in those three
# scripts - wrong genotype label, wrong DST offset, single-format date
# parsing), this script now reads the already-corrected, ZT-enriched
# outputs those scripts write to output/tables/:
#   output/tables/P0_merged_with_ZT.csv       (from 01_)
#   output/tables/F2_old_merged_with_ZT.csv   (from 02_)
#   output/tables/F2_new_merged_with_ZT.csv   (from 03_)
#
# ZT_plot from those files is used directly as the resampling input
# (equivalent to the old hour_of_day_rel - circular() wraps it the same way
# regardless of which range it started in, since it applies modulo 2pi
# internally).
#
# Input file still read fresh here (not produced by 01-03):
#   data/F2_new/eclosiontime_rec0<83..89>_results.csv (for WT_P0 control wells)
# ============================================================================

source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "utils_circ_stats.R"))

rec_ids_new <- 83:89
rec_ids_old <- 60:65

# ---- Reference datasets: read the corrected, ZT-enriched outputs of 01-03 --
p0 <- read.csv(file.path(dir_out_tab, "P0_merged_with_ZT.csv")) %>%
  rename(hour_of_day_rel = ZT_plot)

F2_new <- read.csv(file.path(dir_out_tab, "F2_new_merged_with_ZT.csv")) %>%
  rename(hour_of_day_rel = ZT_plot)

F2_old <- read.csv(file.path(dir_out_tab, "F2_old_merged_with_ZT.csv")) %>%
  rename(hour_of_day_rel = ZT_plot)

# ---- WT_P0 controls added per F2_new recording -------------------------------
WT_P0_controls <- bind_rows(lapply(rec_ids_new, function(recid) {
  read.csv(file.path(dir_data_F2new, paste0("eclosiontime_rec0", recid, "_results.csv")),
           na.strings = c("", "NA")) %>%
    mutate(rec = paste0("rec_0", recid), ecl_time = parse_ecl_time(ecl_time))
})) %>%
  filter(genotype == "WT_P0") %>%
  add_ZT_columns(lights_on_h = lights_on_summer) %>%
  rename(hour_of_day_rel = ZT_plot)

WT_P0_only <- WT_P0_controls %>%
  filter(!is.na(hour_of_day_rel)) %>%
  mutate(group = "WT_P0_control") %>%
  select(hour_of_day_rel, group)

p0_WT_pool <- p0 %>%
  filter(genotype == "wildtype", !is.na(hour_of_day_rel)) %>%
  mutate(group = "p0_Wildtype") %>%
  select(hour_of_day_rel, group)

# ============================================================================
# WT_P0 controls vs. p0 wildtype - pooled and per recording
# ============================================================================
WT_P0_vs_p0_WT <- rbind(WT_P0_only, p0_WT_pool)
result_WT_P0_vs_p0_WT <- resample_2(data = WT_P0_vs_p0_WT, group = "group", period = 24, n_rep = 1000)
print(result_WT_P0_vs_p0_WT$obs)
print(result_WT_P0_vs_p0_WT$pvals)
plot_resample_result(result_WT_P0_vs_p0_WT, "WT_P0 controls vs. p0 Wildtype",
                      save_path = file.path(dir_out_fig, "resampling_WT_P0_vs_p0_WT.pdf"))

resample_per_rec <- function(recid, controls, pool, label) {
  rec_name <- paste0("rec_0", recid)
  this_rec <- controls %>%
    filter(rec == rec_name, !is.na(hour_of_day_rel)) %>%
    mutate(group = label) %>%
    select(hour_of_day_rel, group)
  if (nrow(this_rec) == 0) {
    cat("Skipping", rec_name, "- no data with valid hour_of_day_rel\n")
    return(NULL)
  }
  combined <- rbind(this_rec, pool)
  result   <- resample_2(data = combined, group = "group", period = 24, n_rep = 1000)
  plot_resample_result(result, paste0(rec_name, " - ", label, " vs. p0"))
  data.frame(
    rec = rec_name, n = nrow(this_rec),
    obs_dif_mean = round(result$obs$dif_mean_circ, 3), p_val_circ_mean = round(result$pvals["p_val_circ_mean"], 3),
    obs_dif_sd   = round(result$obs$dif_sd_circ, 3),   p_val_circ_sd   = round(result$pvals["p_val_circ_sd"], 3)
  )
}

pdf(file.path(dir_out_fig, "resampling_WT_P0_control_per_rec.pdf"), width = 10, height = 5)
summary_WT_P0_per_rec <- bind_rows(lapply(rec_ids_new, resample_per_rec,
                                           controls = WT_P0_controls, pool = p0_WT_pool, label = "WT_P0_control"))
dev.off()
print(summary_WT_P0_per_rec)  # includes rec_089
write.csv(summary_WT_P0_per_rec, file.path(dir_out_tab, "summary_WT_P0_per_rec.csv"), row.names = FALSE)

# ============================================================================
# Per-recording resampling for F2 wildtype and mutant, both generations,
# each against the pooled p0 reference. Note the same caveat as before:
# unlike the WT_P0 controls (genetically identical to p0), F2_old/F2_new are
# NOT genetically identical to p0 (recombinant background) - so a difference
# found here can be technical, recording-specific, AND/OR genetic/background-
# driven, and this comparison alone can't separate those explanations. Use
# together with the pooled comparison (below) for that question.
# ============================================================================

p0_MT_pool <- p0 %>%
  filter(genotype == "tim0", !is.na(hour_of_day_rel)) %>%
  mutate(group = "p0_tim0") %>%
  select(hour_of_day_rel, group)

# ---- F2_old wildtype vs. p0 wildtype, per recording (rec_060-065) ----------
F2_old_WT <- F2_old %>%
  filter(genotype == "Wildtype", !is.na(hour_of_day_rel)) %>%
  mutate(group = "F2_old_Wildtype")

pdf(file.path(dir_out_fig, "resampling_F2_old_WT_per_rec.pdf"), width = 10, height = 5)
summary_F2_old_WT_per_rec <- bind_rows(lapply(rec_ids_old, resample_per_rec,
                                               controls = F2_old_WT, pool = p0_WT_pool, label = "F2_old_Wildtype"))
dev.off()
print(summary_F2_old_WT_per_rec)  # includes rec_064, referenced in the thesis text
write.csv(summary_F2_old_WT_per_rec, file.path(dir_out_tab, "summary_F2_old_WT_per_rec.csv"), row.names = FALSE)

# ---- F2_old mutant vs. p0 mutant, per recording (rec_060-065) --------------
F2_old_MT <- F2_old %>%
  filter(genotype == "tim0", !is.na(hour_of_day_rel)) %>%
  mutate(group = "F2_old_tim0")

pdf(file.path(dir_out_fig, "resampling_F2_old_MT_per_rec.pdf"), width = 10, height = 5)
summary_F2_old_MT_per_rec <- bind_rows(lapply(rec_ids_old, resample_per_rec,
                                               controls = F2_old_MT, pool = p0_MT_pool, label = "F2_old_tim0"))
dev.off()
print(summary_F2_old_MT_per_rec)
write.csv(summary_F2_old_MT_per_rec, file.path(dir_out_tab, "summary_F2_old_MT_per_rec.csv"), row.names = FALSE)

# ---- F2_new wildtype vs. p0 wildtype, per recording (rec_083-089) ----------
F2_new_WT <- F2_new %>%
  filter(genotype == "Wildtype", !is.na(hour_of_day_rel)) %>%
  mutate(group = "F2_new_Wildtype")

pdf(file.path(dir_out_fig, "resampling_F2_new_WT_per_rec.pdf"), width = 10, height = 5)
summary_F2_new_WT_per_rec <- bind_rows(lapply(rec_ids_new, resample_per_rec,
                                               controls = F2_new_WT, pool = p0_WT_pool, label = "F2_new_Wildtype"))
dev.off()
print(summary_F2_new_WT_per_rec)
write.csv(summary_F2_new_WT_per_rec, file.path(dir_out_tab, "summary_F2_new_WT_per_rec.csv"), row.names = FALSE)

# ---- F2_new mutant vs. p0 mutant, per recording (rec_083-089) --------------
F2_new_MT <- F2_new %>%
  filter(genotype == "tim0", !is.na(hour_of_day_rel)) %>%
  mutate(group = "F2_new_tim0")

pdf(file.path(dir_out_fig, "resampling_F2_new_MT_per_rec.pdf"), width = 10, height = 5)
summary_F2_new_MT_per_rec <- bind_rows(lapply(rec_ids_new, resample_per_rec,
                                               controls = F2_new_MT, pool = p0_MT_pool, label = "F2_new_tim0"))
dev.off()
print(summary_F2_new_MT_per_rec)
write.csv(summary_F2_new_MT_per_rec, file.path(dir_out_tab, "summary_F2_new_MT_per_rec.csv"), row.names = FALSE)

# ============================================================================
# pooled P0 vs. F2 overlap - histogram overlay (A, B) + circular
# mean +/- circular SD (C), built separately for F2_new and F2_old so both
# generations get their own figure (matches your Results 4.3.3 wording:
# "the eclosion events of each genotype were plotted as an overlap...").
#
# NOTE: this reconstructs and fixes two bugs present in your original
# Resampling_new.R: (1) WT_F2/MT_F2 were referenced but never actually
# assigned (the rbind() results were accidentally stored back into
# WT_data/MT_data instead), and (2) summary_F2 was missing the mean_plot
# column that pC_F2 required (present only in the F2_old version,
# summary_old). Both are fixed below.
# ============================================================================

library(patchwork)

to_ZT12 <- function(x) ifelse(x > 12, x - 24, x)

build_overlap_figure <- function(f2_data, f2_genotype_wt, f2_genotype_mt,
                                  label_wt_f2, label_mt_f2, title_suffix, save_name) {

  WT_data <- rbind(
    p0 %>% filter(genotype == "wildtype", !is.na(hour_of_day_rel)) %>% mutate(group = "P0 tim+") %>% select(hour_of_day_rel, group),
    f2_data %>% filter(genotype == f2_genotype_wt, !is.na(hour_of_day_rel)) %>% mutate(group = label_wt_f2) %>% select(hour_of_day_rel, group)
  ) %>% mutate(ZT_plot = to_ZT12(hour_of_day_rel))

  MT_data <- rbind(
    p0 %>% filter(genotype == "tim0", !is.na(hour_of_day_rel)) %>% mutate(group = "P0 tim01") %>% select(hour_of_day_rel, group),
    f2_data %>% filter(genotype == f2_genotype_mt, !is.na(hour_of_day_rel)) %>% mutate(group = label_mt_f2) %>% select(hour_of_day_rel, group)
  ) %>% mutate(ZT_plot = to_ZT12(hour_of_day_rel))

  WT_hist <- WT_data %>% mutate(bin = floor(ZT_plot)) %>% count(group, bin) %>%
    complete(group, bin = -12:11, fill = list(n = 0)) %>% mutate(ZT_mid = bin + 0.5)
  MT_hist <- MT_data %>% mutate(bin = floor(ZT_plot)) %>% count(group, bin) %>%
    complete(group, bin = -12:11, fill = list(n = 0)) %>% mutate(ZT_mid = bin + 0.5)

  # Computed before pA/pB so their p-values can go straight into the panel
  # subtitles, instead of only living in a separate CSV.
  result_WT <- resample_2(WT_data %>% select(hour_of_day_rel, group), group = "group", period = 24, n_rep = 1000)
  result_MT <- resample_2(MT_data %>% select(hour_of_day_rel, group), group = "group", period = 24, n_rep = 1000)

  fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
  subtitle_WT <- paste0("Resampling P0 vs. ", label_wt_f2, ": mean p ", fmt_p(result_WT$pvals["p_val_circ_mean"]),
                         ", SD p ", fmt_p(result_WT$pvals["p_val_circ_sd"]))
  subtitle_MT <- paste0("Resampling P0 vs. ", label_mt_f2, ": mean p ", fmt_p(result_MT$pvals["p_val_circ_mean"]),
                         ", SD p ", fmt_p(result_MT$pvals["p_val_circ_sd"]))

  pA <- ggplot(WT_hist, aes(x = ZT_mid, y = n, colour = group, group = group)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-12, 12), breaks = c(-12, 0, 12), labels = c("ZT -12", "ZT 0", "ZT 12")) +
    scale_colour_manual(values = setNames(c("#9ecae1", "#08519c"), c("P0 tim+", label_wt_f2))) +
    labs(title = "A  Wildtype (tim+)", subtitle = subtitle_WT, x = "Zeitgeber time (ZT)", y = "Eclosion events") +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 8, color = "grey30"))

  pB <- ggplot(MT_hist, aes(x = ZT_mid, y = n, colour = group, group = group)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-12, 12), breaks = c(-12, 0, 12), labels = c("ZT -12", "ZT 0", "ZT 12")) +
    scale_colour_manual(values = setNames(c("#fdae6b", "#e6550d"), c("P0 tim01", label_mt_f2))) +
    labs(title = "B  tim01 mutant", subtitle = subtitle_MT, x = "Zeitgeber time (ZT)", y = "Eclosion events") +
    theme_classic() +
    theme(plot.subtitle = element_text(size = 8, color = "grey30"))

  # Built from the actual group_0/group_1 names returned by
  # summarise_circ_stats() rather than assuming P0 is always first - that
  # assumption was wrong because group_by() sorts alphabetically, and
  # "F2_..." sorts before "P0..." (F < P), which silently swapped the P0
  # and F2 values in this panel. Confirmed by comparing the rendered
  # figure against /18 - the labels didn't match the correct values.
  summary_df <- bind_rows(
    data.frame(group = result_WT$obs$group_0, mean_circ = result_WT$obs$mean_circ_0, sd_circ = result_WT$obs$sd_circ_0),
    data.frame(group = result_WT$obs$group_1, mean_circ = result_WT$obs$mean_circ_1, sd_circ = result_WT$obs$sd_circ_1),
    data.frame(group = result_MT$obs$group_0, mean_circ = result_MT$obs$mean_circ_0, sd_circ = result_MT$obs$sd_circ_0),
    data.frame(group = result_MT$obs$group_1, mean_circ = result_MT$obs$mean_circ_1, sd_circ = result_MT$obs$sd_circ_1)
  ) %>%
    mutate(
      group     = factor(group, levels = c("P0 tim+", label_wt_f2, "P0 tim01", label_mt_f2)),
      mean_plot = ifelse(mean_circ > 12, mean_circ - 24, mean_circ)
    )

  # Dynamic y-limits: a fixed +/-12h window clips the error bars when SD is
  # large (e.g. F2_new tim01 has SD ~9.6h, near-uniform distribution, since
  # its Rayleigh test wasn't significant - see table_19). Pad the actual
  # mean+/-SD range by 10% instead of assuming +/-12h always fits.
  y_range <- range(c(summary_df$mean_plot - summary_df$sd_circ, summary_df$mean_plot + summary_df$sd_circ))
  y_pad   <- diff(y_range) * 0.1
  y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  y_breaks <- pretty(y_limits, n = 5)

  pC <- ggplot(summary_df, aes(x = group, y = mean_plot)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_plot - sd_circ, ymax = mean_plot + sd_circ), width = 0.2) +
    geom_text(aes(y = mean_plot + sd_circ, label = sprintf("ZT %.2f \u00b1 %.2f", mean_plot, sd_circ)),
              vjust = -0.6, size = 2.8, color = "grey20") +
    scale_y_continuous(breaks = y_breaks, labels = function(x) paste0("ZT ", x)) +
    coord_cartesian(ylim = c(y_limits[1], y_limits[2] + diff(y_limits) * 0.08)) +
    labs(title = "C  Circular mean \u00b1 circular SD", x = "", y = "Eclosion time (ZT)") +
    theme_classic()

  final_fig <- (pA | pB) / pC +
    plot_annotation(title = paste0("Effect of genomic background on eclosion rhythm - ", title_suffix))
  # print(final_fig) removed - ggsave() below saves the PDF independently,
  # and print() to the RStudio plot pane was causing "grid.newpage(): write
  # failed" errors on this machine. Open the saved PDF directly to view it.
  ggsave(file.path(dir_out_fig, save_name), final_fig, width = 10, height = 12)

  pvals_df <- data.frame(
    comparison       = c(paste0("P0 vs. ", label_wt_f2), paste0("P0 vs. ", label_mt_f2)),
    obs_dif_mean     = c(result_WT$obs$dif_mean_circ, result_MT$obs$dif_mean_circ),
    p_val_circ_mean  = c(result_WT$pvals["p_val_circ_mean"], result_MT$pvals["p_val_circ_mean"]),
    obs_dif_sd       = c(result_WT$obs$dif_sd_circ, result_MT$obs$dif_sd_circ),
    p_val_circ_sd    = c(result_WT$pvals["p_val_circ_sd"], result_MT$pvals["p_val_circ_sd"])
  )
  print(pvals_df)

  list(figure = final_fig, summary = summary_df, pvals = pvals_df)
}

fig17_F2_new <- build_overlap_figure(
  f2_data = F2_new, f2_genotype_wt = "Wildtype", f2_genotype_mt = "tim0",
  label_wt_f2 = "F2_new tim+", label_mt_f2 = "F2_new tim01",
  title_suffix = "F2 new (rec_083-089)", save_name = "P0_vs_F2_new_overlap_summary.pdf"
)

fig17_F2_old <- build_overlap_figure(
  f2_data = F2_old, f2_genotype_wt = "Wildtype", f2_genotype_mt = "tim0",
  label_wt_f2 = "F2_old tim+", label_mt_f2 = "F2_old tim01",
  title_suffix = "F2 old (rec_060-065)", save_name = "P0_vs_F2_old_overlap_summary.pdf"
)

write.csv(fig17_F2_new$summary, file.path(dir_out_tab, "P0_vs_F2_new_circular_summary.csv"), row.names = FALSE)
write.csv(fig17_F2_old$summary, file.path(dir_out_tab, "P0_vs_F2_old_circular_summary.csv"), row.names = FALSE)

# Resampling p-values for the P0-vs-F2 comparisons shown in # (statistical backing for what the panel C error bars show visually)
write.csv(fig17_F2_new$pvals, file.path(dir_out_tab, "P0_vs_F2_new_resampling_pvals.csv"), row.names = FALSE)
write.csv(fig17_F2_old$pvals, file.path(dir_out_tab, "P0_vs_F2_old_resampling_pvals.csv"), row.names = FALSE)

# ============================================================================
# OPTIONAL / not confirmed in results.docx: p0_JI (Jun's winter-time P0
# reference, rec_031) vs. F2_new mutant. This comparison exists in your
# original Resampling_new.R but I could not confirm it appears anywhere in
# your results.docx - check with Luisa/Jun whether it belongs in the thesis
# before relying on it. Left here, clearly separated, rather than silently
# dropped or silently kept.
# ============================================================================
# p0_JI <- read.csv(file.path(dir_data_P0, "rec_031_result.csv")) %>%
#   mutate(
#     ecl_time_local  = force_tz(parse_ecl_time(ecl_time), "Europe/Berlin"),
#     hour_of_day     = hour(ecl_time_local) + minute(ecl_time_local) / 60,
#     hour_of_day_rel = hour_of_day - lights_on_winter
#   )
# p0_JI_MT <- p0_JI %>% filter(line == "BL-80930", !is.na(hour_of_day_rel)) %>% mutate(group = "p0_JI") %>% select(hour_of_day_rel, group)
# F2_new_MT <- F2_new %>% filter(genotype == "tim0", !is.na(hour_of_day_rel)) %>% mutate(group = "F2_new") %>% select(hour_of_day_rel, group)
# result_MT_JI <- resample_2(rbind(p0_JI_MT, F2_new_MT), group = "group", period = 24, n_rep = 1000)
# plot_resample_result(result_MT_JI, "Mutant: p0 (JI) vs. F2 new",
#                       save_path = file.path(dir_out_fig, "resampling_MT_JI.pdf"))
