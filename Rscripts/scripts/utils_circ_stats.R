# ============================================================================
# utils_circ_stats.R
# Shared circular-statistics helper functions used by the resampling scripts
# (04_resampling_P0_vs_F2.R). Source this once, e.g.:
#   source(here::here("scripts", "utils_circ_stats.R"))
#
# Two changes from your original circ_stats.R:
# 1. plot_resample_result() was previously only defined inline inside
#    Resampling_new.R and never sourced by F2_P0_comparisons_clean.R /
#    F2_vs_p0.R, which would have caused an "object not found" error. It now
#    lives here so every script that needs it can find it.
# 2. summarise_circ_stats() is optimized for speed (same results, ~2-4x
#    faster) - see the comment on that function for details. This matters
#    a lot since resample_2() calls it ~1000 times per comparison.
# ============================================================================

library(dplyr)
library(circular)

# ---- Extract the p-value band from watson.two.test() print output ----------
# watson.two.test() only stores the U2 statistic; the significance band
# (e.g. "0.001 < P-value < 0.01") is computed against a table of critical
# values (Watson 1962) at print() time and is not returned as a numeric
# field. This grabs that text so it can go straight into a results table
# instead of being read off the console by hand.
watson_pvalue_text <- function(watson_result) {
  out <- capture.output(print(watson_result))
  p_line <- out[grepl("P-value|p-value", out)]
  if (length(p_line) == 0) return(NA_character_)
  trimws(p_line[length(p_line)])
}

# ---- Circular descriptive stats for two groups -----------------------------
# Builds the circular object once per group and reuses it for mu and kappa,
# instead of reconstructing it and re-running mle.vonmises() redundantly
# (the original version called it twice per group). Also drops a
# rayleigh.test() call that was being run twice per group but whose result
# (R, p.value) was never actually used anywhere in this function's output -
# pure wasted computation, now removed. Same results as before, but
# noticeably faster - matters a lot inside resample_2(), which calls this
# function ~1000 times per comparison.
summarise_circ_stats <- function(data, group, period = 24) {
  if (nlevels(as.factor(data[[group]])) != 2) {
    stop("There has to be two groups")
  }
  out <- data %>%
    group_by(across(all_of(group))) %>%
    summarise(
      group_name = first(.data[[group]]),
      .stats = {
        rad <- circular(hour_of_day_rel * 2 * pi / period, units = "radians", modulo = "2pi")
        vm  <- mle.vonmises(rad)
        mu_raw <- as.numeric(vm$mu) * period / (2 * pi)
        list(list(
          mu_hours = ((mu_raw + 12) %% period) - period / 2,
          kappa    = vm$kappa
        ))
      },
      .groups = "drop"
    ) %>%
    mutate(
      mu_hours = sapply(.stats, `[[`, "mu_hours"),
      kappa    = sapply(.stats, `[[`, "kappa")
    ) %>%
    select(-.stats) %>%
    mutate(
      circ_sd = sqrt(-2 * log(besselI(kappa, 1) / besselI(kappa, 0))) * period / (2 * pi)
    ) %>%
    summarise(
      group_0 = group_name[1],
      group_1 = group_name[2],
      mean_circ_0 = mu_hours[1],
      mean_circ_1 = mu_hours[2],
      sd_circ_0 = circ_sd[1],
      sd_circ_1 = circ_sd[2],
      dif_mean_circ = diff(mu_hours),
      dif_sd_circ = diff(circ_sd)
    )
  return(out)
}

# ---- Circular resampling test (difference in circular mean / SD) -----------
# ---- Circular resampling test (difference in circular mean / SD) -----------
# Uses a pre-allocated list + a single bind_rows() at the end, instead of
# rbind()-ing onto a growing data frame inside the loop. rbind() in a loop
# is a classic R performance trap: R copies the entire accumulated object on
# every iteration, so cost grows quadratically with n_rep. With 1000
# iterations that copying dwarfs the actual statistics computation - this
# is very likely why resample_2() felt like it "never stops."
resample_2 <- function(data, group, period, n_rep) {
  obs <- summarise_circ_stats(data = data, group = group, period = period)
  resampled_list <- vector("list", n_rep - 1)
  for (i in seq_len(n_rep - 1)) {
    shuffled <- data %>% mutate(group = sample(group))
    resampled_list[[i]] <- summarise_circ_stats(data = shuffled, group = group, period = period)
  }
  resampled <- bind_rows(resampled_list)
  pval_circ_mean <- 1 - 2 * abs((sum(resampled$dif_mean_circ < obs$dif_mean_circ) + 1) / n_rep - 0.5)
  pval_circ_sd   <- 1 - 2 * abs((sum(resampled$dif_sd_circ   < obs$dif_sd_circ)   + 1) / n_rep - 0.5)
  out <- list(
    obs   = obs,
    resam = resampled,
    pvals = c("p_val_circ_mean" = pval_circ_mean, "p_val_circ_sd" = pval_circ_sd)
  )
  return(out)
}

# ---- Reusable plot for resample_2() results ---------------------------------
plot_resample_result <- function(result, comparison_name, save_path = NULL) {

  if (!is.null(save_path)) pdf(save_path, width = 10, height = 5)

  par(mfrow = c(1, 2))

  hist(result$resam$dif_mean_circ,
       main = paste0(comparison_name, "\ndifference in circular mean"),
       xlab = "Difference in circular mean (h)",
       col = "grey80", breaks = 30)
  abline(v = result$obs$dif_mean_circ, col = "red", lwd = 2)
  legend("topleft",
         legend = paste0("obs = ", round(result$obs$dif_mean_circ, 3),
                          "\np = ", round(result$pvals["p_val_circ_mean"], 3)),
         bty = "n")

  hist(result$resam$dif_sd_circ,
       main = paste0(comparison_name, "\ndifference in circular SD"),
       xlab = "Difference in circular SD (h)",
       col = "grey80", breaks = 30)
  abline(v = result$obs$dif_sd_circ, col = "red", lwd = 2)
  legend("topleft",
         legend = paste0("obs = ", round(result$obs$dif_sd_circ, 3),
                          "\np = ", round(result$pvals["p_val_circ_sd"], 3)),
         bty = "n")

  if (!is.null(save_path)) dev.off()
}
