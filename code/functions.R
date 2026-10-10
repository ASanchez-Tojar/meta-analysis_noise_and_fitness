#-----------------------heterogeneity interpretation-----------------------#
het_interpret <- function(observed_value, het_type, es_type, data) {
  
  if (!het_type %in% c("I2", "CVH", "M", "sigma2", "V_bar")) {
    stop("Invalid heterogeneity type. Choose from 'I2', 'CVH', 'M', 'sigma2', or 'V_bar'.")
  }
  
  if (!es_type %in% unique(data$es.type)) {
    stop("Invalid effect size type. Check the levels of 'es.type'.")
  }
  
  filtered_data <- data %>% filter(es.type == es_type)
  
  if (nrow(filtered_data) == 0) {
    stop("No data available for the selected effect size type.")
  }
  
  het_values <- filtered_data[[het_type]]
  
  percentiles <- quantile(het_values, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
  
  lower_bound <- max(percentiles[percentiles <= observed_value], na.rm = TRUE)
  upper_bound <- min(percentiles[percentiles >= observed_value], na.rm = TRUE)
  lower_percentile <- names(percentiles)[percentiles == lower_bound]
  upper_percentile <- names(percentiles)[percentiles == upper_bound]
  
  percentile_range <- paste0(lower_percentile, "-", upper_percentile, "th percentile")
  
  # return the results
  return(tibble::tibble(
    observed_value = observed_value,
    het_type = het_type,
    es_type = es_type,
    percentile_range = percentile_range
  ))
}

#-----------------------estimates extraction -----------------------#
extract_model <- function(m, moderator = NULL,
                          es = "ES_ID", study = "Study_ID") {
  
  # Coefficient table
  tab <- coef(summary(m))
  tab$term <- rownames(tab)
  tab <- tab[, c("term", "estimate", "pval", "ci.lb", "ci.ub")]
  
  # Data actually used by the model (drops rows with NAs)
  dat <- m$data[m$not.na, ]
  
  if (is.null(moderator)) {
    # Intercept-only: all data points belong to the single estimate
    tab$k_ES <- n_distinct(dat[[es]])
    tab$n_Study <- n_distinct(dat[[study]])
  } else {
    counts <- dat |>
      group_by(level = .data[[moderator]]) |>
      summarise(k_ES = n_distinct(.data[[es]]),
                n_Study = n_distinct(.data[[study]]),
                .groups = "drop") |>
      mutate(term = paste0(moderator, level)) |>
      select(term, k_ES, n_Study)
    tab <- left_join(tab, counts, by = "term")
  }
  rownames(tab) <- NULL
  tab
}