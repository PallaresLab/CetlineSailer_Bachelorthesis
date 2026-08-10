# ============================================================================
# 05_appendix_P0_without_air.R
# Appendix eclosion timing of the three P0 recordings that were
# run without active air supply (rec_053, rec_056, rec_057), included as a
# comparison/control against the main air-supplied P0 recordings.
#
# Input files expected in data/P0/without_air/:
#   eclosiontime rec_053_results.csv   (wildtype, all 3 plates)
#   eclosiontime rec_056_results.csv   (mutant,   plates 1-2)
#   eclosiontime rec_057_results.csv   (mutant,   plate 3)
# ============================================================================

source(here::here("scripts", "00_setup.R"))

dir_data_P0_noair <- file.path(dir_data_P0, "without_air")

files_noair <- c(
  file.path(dir_data_P0_noair, "eclosiontime rec_053_results.csv"),
  file.path(dir_data_P0_noair, "eclosiontime rec_056_results.csv"),
  file.path(dir_data_P0_noair, "eclosiontime rec_057_results.csv")
)

genotype_map_noair <- data.frame(
  rec       = c("rec_053", "rec_053", "rec_053", "rec_056", "rec_056", "rec_057"),
  plate_num = c(1, 2, 3, 1, 2, 3),
  genotype  = c("wildtype", "wildtype", "wildtype", "tim0", "tim0", "tim0"),
  stringsAsFactors = FALSE
)

# ---- Robust datetime parser (some files use the alternate M/D/YYYY format) -
parse_ecl_time_ymd_or_mdy <- function(x) {
  parsed <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "CET")
  if (all(is.na(parsed))) parsed <- as.POSIXct(x, format = "%m/%d/%Y %H:%M", tz = "CET")
  parsed
}

all_res_noair <- bind_rows(lapply(files_noair, function(f) {
  df <- read.csv(f)
  df$ecl_time <- parse_ecl_time_ymd_or_mdy(as.character(df$ecl_time))
  df$rec <- regmatches(f, regexpr("rec_\\d+", f))
  df
})) %>%
  mutate(plate_num = as.numeric(sub("_.*", "", plate_well))) %>%
  left_join(genotype_map_noair, by = c("rec", "plate_num"))

# ---- Eclosion numbers per genotype -------------------------------------------
all_res_noair %>%
  group_by(genotype) %>%
  summarise(total = n(), with_time = sum(!is.na(ecl_time)), .groups = "drop") %>%
  print()

# ---- P0 without air supply, ZT histogram --------------------------
plot_data_noair <- all_res_noair %>%
  filter(!is.na(ecl_time), !is.na(genotype)) %>%
  mutate(genotype = factor(genotype, levels = c("wildtype", "tim0"), labels = c("Wildtype", "tim0"))) %>%
  add_ZT_columns(lights_on_h = lights_on_winter)

n_summary_noair <- plot_data_noair %>% count(genotype)

fig_21 <- plot_data_noair %>%
  count(genotype, ZT_bin, daytime) %>%
  complete(genotype, ZT_bin, fill = list(n = 0, daytime = NA)) %>%
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
  facet_wrap(~ genotype, ncol = 2,
             labeller = as_labeller(c(
               "Wildtype" = paste0("Wildtype (n == ", n_summary_noair$n[n_summary_noair$genotype == "Wildtype"], ")"),
               "tim0"     = paste0("italic(tim)^{'01'}~(n == ", n_summary_noair$n[n_summary_noair$genotype == "tim0"], ")")
             ), label_parsed)) +
  labs(fill = "") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(size = 12, face = "bold"),
    legend.position   = "bottom",
    panel.border      = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )
fig_21
ggsave(file.path(dir_out_fig, "P0_without_air_supply.pdf"), fig_21, width = 10, height = 5)

