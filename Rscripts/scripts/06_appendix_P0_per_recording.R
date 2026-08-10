# ============================================================================
# 06_appendix_P0_per_recording.R
# Appendix eclosion timing of the (air-supplied) P0 recordings,
# broken down per recording (rec_067, rec_068, rec_069).
#
# rec_069 needs a +1h daylight-saving-time correction (its timestamps were
# recorded 1 hour behind the other two recordings) - see memory note from
# the June/July debugging pass.
#
# Input files expected in data/P0/:
#   eclosiontime rec_067_results.csv   (wildtype, both plates)
#   eclosiontime rec_068_results.csv   (mutant,   both plates)
#   eclosiontime rec_069_results.csv   (plate 1 = wildtype, plate 2 = mutant)
# ============================================================================

source(here::here("scripts", "00_setup.R"))

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

all_res_P0 <- bind_rows(lapply(files_P0, function(f) {
  df <- read.csv(f)
  df$ecl_time <- parse_ecl_time(df$ecl_time)
  df$rec <- regmatches(f, regexpr("rec_\\d+", f))
  # DST correction: rec_069 timestamps are 1h behind the other two recordings
  if (grepl("rec_069", f)) df$ecl_time <- df$ecl_time + hours(1)
  df
})) %>%
  mutate(plate_num = as.numeric(sub("_.*", "", plate_well))) %>%
  left_join(genotype_map_P0, by = c("rec", "plate_num"))

# ---- ZT histogram, faceted by recording ----------------------------
plot_data_P0 <- all_res_P0 %>%
  filter(!is.na(ecl_time), !is.na(genotype)) %>%
  mutate(genotype = factor(genotype, levels = c("wildtype", "tim0"), labels = c("Wildtype", "tim0"))) %>%
  add_ZT_columns(lights_on_h = lights_on_summer)

fig_22 <- plot_data_P0 %>%
  count(rec, genotype, ZT_bin, daytime) %>%
  complete(rec, genotype, ZT_bin, fill = list(n = 0, daytime = NA)) %>%
  mutate(daytime = ifelse(ZT_bin >= 0 & ZT_bin < 12, "Day", "Night")) %>%
  ggplot(aes(x = ZT_bin + bin_size / 2, y = n, fill = daytime)) +
  geom_vline(xintercept = 0,  linetype = "dashed", color = "red",       linewidth = 0.5) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "steelblue", linewidth = 0.5) +
  annotate("text", x = -0.5, y = Inf, label = "Lights on",  vjust = 2, color = "red",       size = 3, hjust = 1) +
  annotate("text", x = 11.7, y = Inf, label = "Lights off", vjust = 2, color = "steelblue", size = 3, hjust = 1) +
  geom_col(width = bin_size, position = "identity", color = "white", linewidth = 0.1, alpha = 0.85) +
  scale_fill_manual(values = c("Day" = col_day, "Night" = col_night), na.value = "grey80") +
  scale_x_continuous(name = "Zeitgeber Time (ZT, h)", breaks = seq(-12, 12, by = 12),
                      limits = c(-12, 13.5), labels = function(x) paste0("ZT ", x)) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(name = "Number of eclosions", limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  facet_grid(genotype ~ rec,
             labeller = labeller(
               genotype = as_labeller(c("Wildtype" = "Wildtype", "tim0" = "italic(tim)^{'01'}"), label_parsed),
               rec = label_value
             )) +
  labs(title = "Eclosion time distribution per recording", fill = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 9, face = "bold"),
    legend.position   = "bottom",
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )
fig_22
ggsave(file.path(dir_out_fig, "P0_per_recording.pdf"), fig_22, width = 12, height = 6)
