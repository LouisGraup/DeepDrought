# Join Pfyn crown condition observations with tree spatial and treatment information.

required_packages <- c("dplyr", "ggplot2", "patchwork", "readr", "sf")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
})

find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)

  repeat {
    if (dir.exists(file.path(path, ".git"))) {
      return(path)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the repository root from: ", start, call. = FALSE)
    }

    path <- parent
  }
}

repo_root <- find_repo_root()

csv_path <- file.path(repo_root, "Data/Pfyn/PFY_crown_condition.csv")
lai_path <- file.path(repo_root, "Data/Pfyn/Pfyn_LAI04_21.csv")
shp_path <- file.path(
  repo_root,
  "Data/Pfyn/pfy_v_tree_info_spaJoin_stop/pfy_v_tree_info_spaJoin_stop.shp"
)

output_dir <- file.path(repo_root, "Data/Pfyn")
plot_dir <- file.path(repo_root, "Figures/Pfyn")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

spatial_csv_path <- file.path(output_dir, "PFY_crown_condition_spatial.csv")
spatial_rds_path <- file.path(output_dir, "PFY_crown_condition_spatial_matched.rds")
unmatched_banr_path <- file.path(output_dir, "PFY_crown_condition_spatial_unmatched_banr.csv")
unmatched_xy_path <- file.path(output_dir, "PFY_crown_condition_spatial_unmatched_xy.csv")
summary_path <- file.path(output_dir, "PFY_crown_condition_summary_by_year_plot.csv")
crown_summary_path <- file.path(output_dir, "PFY_crown_transparency_by_treatment.csv")
crown_plot_path <- file.path(plot_dir, "PFY_crown_transparency_by_treatment.png")
crown_plot_se_path <- file.path(plot_dir, "PFY_crown_transparency_by_treatment_se.png")
defoliation_change_summary_path <- file.path(output_dir, "PFY_annual_defoliation_change_by_treatment.csv")
defoliation_change_plot_path <- file.path(plot_dir, "PFY_annual_defoliation_change_by_treatment.png")
defoliation_mortality_plot_path <- file.path(plot_dir, "PFY_defoliation_change_mortality_by_treatment.png")
defoliation_mortality_full_period_plot_path <- file.path(
  plot_dir,
  "PFY_defoliation_change_mortality_by_treatment_2004_2024.png"
)
mortality_summary_path <- file.path(output_dir, "PFY_mortality_fraction_by_treatment.csv")
mortality_plot_path <- file.path(plot_dir, "PFY_mortality_fraction_by_treatment.png")
lai_crown_plot_data_path <- file.path(output_dir, "PFY_LAI_crown_plot_level_join.csv")
lai_crown_summary_path <- file.path(output_dir, "PFY_LAI_crown_by_treatment_year.csv")
lai_crown_correlation_path <- file.path(output_dir, "PFY_LAI_crown_correlations_by_treatment.csv")
lai_crown_timeseries_path <- file.path(plot_dir, "PFY_LAI_crown_timeseries_by_treatment.png")
lai_crown_scatter_path <- file.path(plot_dir, "PFY_LAI_crown_correlation_by_treatment.png")

scenario_colors <- c(
  "Control" = "#E69F00",
  "Irrigation" = "#56B4E9",
  "Irrigation stop" = "#009E73"
)
pfyn_text_theme <- theme(
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  strip.text = element_text(size = 12, face = "bold")
)
lai_comparison_treatments <- c("Control", "Irrigation")
crown_timeseries_treatments <- c("Control", "Irrigation", "Irrigation stop")
crown_timeseries_start_year <- 2003
crown_timeseries_end_year <- 2024
plot_start_year <- 2014
plot_end_year <- 2024
defoliation_change_start_year <- plot_start_year
defoliation_change_end_year <- plot_end_year
mortality_start_year <- 2014
mortality_end_year <- plot_end_year
full_period_start_year <- 2004
full_period_end_year <- 2024
irrigation_stop_year <- 2014
source_crown_transparency_column <- paste0("ALL", "N", "BV")
source_crown_transparency_change_column <- paste0(source_crown_transparency_column, "_DIFF")
source_crown_transparency_mort_cor_column <- paste0(
  source_crown_transparency_column,
  "_MORT_COR"
)
source_crown_transparency_change_mort_cor_column <- paste0(
  source_crown_transparency_change_column,
  "_MORT_COR"
)
source_needle_loss_column <- paste0("N", "BV")
source_needle_loss_change_column <- paste0(source_needle_loss_column, "_DIFF")
source_needle_loss_mort_cor_column <- paste0(
  source_needle_loss_column,
  "_MORT_COR"
)
source_needle_loss_change_mort_cor_column <- paste0(
  source_needle_loss_change_column,
  "_MORT_COR"
)
source_condition_measure_columns <- c(
  source_crown_transparency_column,
  source_crown_transparency_change_column,
  source_crown_transparency_mort_cor_column,
  source_crown_transparency_change_mort_cor_column,
  source_needle_loss_column,
  source_needle_loss_change_column,
  source_needle_loss_mort_cor_column,
  source_needle_loss_change_mort_cor_column
)
source_shp_condition_column_renames <- c(
  "needle_loss",
  "crown_transparency",
  "crown_transparency_original"
)
names(source_shp_condition_column_renames) <- c(
  tolower(source_needle_loss_column),
  tolower(source_crown_transparency_column),
  paste0(tolower(source_crown_transparency_column), "o")
)

message("Reading crown condition table: ", csv_path)
crown_trans <- read_csv(
  csv_path,
  na = c("", "NA"),
  col_types = cols(
    BADATUM = col_character(),
    .default = col_guess()
  ),
  guess_max = 100000,
  show_col_types = FALSE
)

required_crown_columns <- c(
  "INV_KIND",
  "YEAR",
  "BANR",
  "X",
  "Y",
  source_condition_measure_columns,
  "MORTALITY_HEALTH",
  "CUT"
)
missing_crown_columns <- setdiff(required_crown_columns, names(crown_trans))
if (length(missing_crown_columns) > 0) {
  stop(
    "The crown condition table is missing required columns: ",
    paste(missing_crown_columns, collapse = ", "),
    call. = FALSE
  )
}

crown_trans <- crown_trans |>
  rename(
    crown_transparency = all_of(source_crown_transparency_column),
    crown_transparency_change = all_of(source_crown_transparency_change_column),
    crown_transparency_mort_cor = all_of(source_crown_transparency_mort_cor_column),
    crown_transparency_change_mort_cor = all_of(source_crown_transparency_change_mort_cor_column),
    needle_loss = all_of(source_needle_loss_column),
    needle_loss_change = all_of(source_needle_loss_change_column),
    needle_loss_mort_cor = all_of(source_needle_loss_mort_cor_column),
    needle_loss_change_mort_cor = all_of(source_needle_loss_change_mort_cor_column)
  )

raw_crown_condition_rows <- nrow(crown_trans)
crown_trans <- crown_trans |>
  filter(INV_KIND == "primary")
non_primary_crown_condition_rows <- raw_crown_condition_rows - nrow(crown_trans)

message("Reading LAI table: ", lai_path)
lai <- read_csv(
  lai_path,
  col_types = cols(
    year = col_double(),
    block = col_double(),
    plot = col_double(),
    treatment = col_character(),
    LAI = col_double()
  ),
  show_col_types = FALSE
)

message("Reading spatial tree shapefile: ", shp_path)
tree_spatial <- st_read(shp_path, quiet = TRUE)

stopifnot(all(c("year", "block", "plot", "treatment", "LAI") %in% names(lai)))
stopifnot(all(c("banr", "banreti", "x", "y", "PLOT", "TYPE") %in% names(tree_spatial)))

geom_xy <- st_coordinates(tree_spatial)

shp_treatment <- function(type) {
  case_when(
    type == paste0("nicht bew", "\u00e4", "ssert") ~ "Control",
    type == paste0("bew", "\u00e4", "ssert") ~ "Irrigation",
    type == paste0("Bew", "\u00e4", "sserungsstopp") ~ "Irrigation stop",
    type == "Pufferzone" ~ "Buffer",
    .default = NA_character_
  )
}

mortality_rate_from_flags <- function(mortality_flag, cut_flag) {
  mortality_flag <- as.numeric(mortality_flag)
  cut_flag <- !is.na(cut_flag) & cut_flag > 0
  valid <- !is.na(mortality_flag)
  disappeared_without_counted_death <- valid & cut_flag & mortality_flag == 0
  at_risk <- valid & !disappeared_without_counted_death

  if (sum(at_risk) == 0) {
    return(NA_real_)
  }

  sum(mortality_flag[at_risk] > 0, na.rm = TRUE) / sum(at_risk)
}

expand_range <- function(values, pad_fraction = 0.08) {
  range_values <- range(values, na.rm = TRUE)

  if (!all(is.finite(range_values))) {
    return(c(0, 1))
  }

  range_width <- diff(range_values)
  if (range_width == 0) {
    return(range_values + c(-0.5, 0.5))
  }

  range_values + c(-1, 1) * range_width * pad_fraction
}

mortality_column <- "MORTALITY_HEALTH"
mortality_definition_label <- "Health-related mortality"

# The shapefile geometry is WGS84 lon/lat, while its x/y attributes match
# the crown condition table's X/Y coordinates. Join by tree id because one known tree
# has slightly different coordinates but the same BANR.
tree_lookup <- tree_spatial |>
  rename_with(
    .fn = function(column_name) unname(source_shp_condition_column_renames[column_name]),
    .cols = any_of(names(source_shp_condition_column_renames))
  ) |>
  mutate(
    row_id = row_number(),
    treatment_broad = case_when(
      is.na(TYPE) & PLOT %in% c(0, 9) ~ "Buffer",
      .default = shp_treatment(TYPE)
    ),
    geom_lon = geom_xy[, "X"],
    geom_lat = geom_xy[, "Y"],
    geometry_wkt = st_as_text(geometry)
  ) |>
  st_drop_geometry() |>
  mutate(join_BANR = banr)

duplicate_banr <- tree_lookup |>
  count(join_BANR, name = "n") |>
  filter(n > 1)

if (nrow(duplicate_banr) > 0) {
  stop(
    "The shapefile has duplicate banr values. ",
    "Inspect `duplicate_banr` before using BANR as the join key.",
    call. = FALSE
  )
}

tree_lookup <- tree_lookup |>
  rename_with(.fn = function(column_name) paste0("shp_", column_name), .cols = -join_BANR)

crown_trans_spatial <- crown_trans |>
  left_join(tree_lookup, by = c("BANR" = "join_BANR")) |>
  mutate(
    spatial_match = !is.na(shp_row_id),
    coordinate_exact_match = spatial_match & X == shp_x & Y == shp_y,
    coordinate_distance_m = if_else(
      spatial_match,
      sqrt((X - shp_x)^2 + (Y - shp_y)^2),
      NA_real_
    ),
    analysis_Treatment = shp_treatment_broad,
    analysis_Treatment_source = case_when(
      !is.na(shp_treatment_broad) ~ "shapefile",
      .default = NA_character_
    ),
    analysis_Plot = case_when(
      analysis_Treatment == "Buffer" & shp_PLOT %in% c(0, 9) ~ 0,
      .default = shp_PLOT
    ),
    analysis_include = analysis_Treatment %in% c("Control", "Irrigation", "Irrigation stop")
  ) |>
  relocate(
    spatial_match,
    coordinate_exact_match,
    coordinate_distance_m,
    analysis_Plot,
    analysis_Treatment,
    analysis_Treatment_source,
    analysis_include,
    .after = Y
  )

unmatched_banr <- crown_trans_spatial |>
  filter(!spatial_match) |>
  group_by(BANR, X, Y) |>
  summarise(
    n_rows = n(),
    first_year = min(YEAR, na.rm = TRUE),
    last_year = max(YEAR, na.rm = TRUE),
    first_CLNR = first(CLNR),
    .groups = "drop"
  ) |>
  arrange(BANR)

unmatched_xy <- crown_trans_spatial |>
  filter(!coordinate_exact_match) |>
  group_by(
    BANR,
    X,
    Y,
    spatial_match,
    shp_banr,
    shp_x,
    shp_y,
    coordinate_distance_m
  ) |>
  summarise(
    n_rows = n(),
    first_year = min(YEAR, na.rm = TRUE),
    last_year = max(YEAR, na.rm = TRUE),
    first_CLNR = first(CLNR),
    .groups = "drop"
  ) |>
  arrange(BANR)

summary_by_year_plot <- crown_trans_spatial |>
  filter(spatial_match, analysis_include) |>
  group_by(YEAR, analysis_Treatment, analysis_Plot, shp_TEXT) |>
  summarise(
    n_records = n(),
    n_trees = n_distinct(BANR),
    mean_crown_transparency = mean(crown_transparency, na.rm = TRUE),
    sd_crown_transparency = sd(crown_transparency, na.rm = TRUE),
    mortality_rate_health_related = mortality_rate_from_flags(MORTALITY_HEALTH, CUT),
    n_cut = sum(CUT > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(YEAR, analysis_Treatment, analysis_Plot, shp_TEXT)

crown_summary <- crown_trans_spatial |>
  filter(
    spatial_match,
    analysis_include,
    YEAR >= plot_start_year,
    YEAR <= plot_end_year,
    !is.na(crown_transparency)
  ) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = c("Control", "Irrigation", "Irrigation stop")
    )
  ) |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_records = n(),
    n_trees = n_distinct(BANR),
    mean_crown_transparency = mean(crown_transparency, na.rm = TRUE),
    sd_crown_transparency = sd(crown_transparency, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(se_crown_transparency = sd_crown_transparency / sqrt(n_records)) |>
  arrange(YEAR, analysis_Treatment)

plot_crown_transparency <- function(summary_data, error_column) {
  ggplot(
    summary_data,
    aes(x = YEAR, y = mean_crown_transparency, color = analysis_Treatment, group = analysis_Treatment)
  ) +
    geom_errorbar(
      aes(
        ymin = mean_crown_transparency - .data[[error_column]],
        ymax = mean_crown_transparency + .data[[error_column]]
      ),
      width = 0.25,
      alpha = 0.75,
      linewidth = 0.5
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = scenario_colors, drop = FALSE) +
    scale_x_continuous(breaks = sort(unique(summary_data$YEAR))) +
    labs(
      x = "",
      y = "Crown transparency (%)",
      color = "Treatment"
    ) +
    theme_bw() +
    pfyn_text_theme +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

crown_transparency_plot <- plot_crown_transparency(
  crown_summary,
  error_column = "sd_crown_transparency"
)

crown_transparency_se_plot <- plot_crown_transparency(
  crown_summary,
  error_column = "se_crown_transparency"
)

defoliation_change_data <- crown_trans_spatial |>
  filter(
    spatial_match,
    analysis_include,
    YEAR >= defoliation_change_start_year,
    YEAR <= defoliation_change_end_year,
    !is.na(crown_transparency_change_mort_cor)
  ) |>
  transmute(
    YEAR,
    BANR,
    analysis_Treatment,
    annual_change = crown_transparency_change_mort_cor
  ) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = c("Control", "Irrigation", "Irrigation stop")
    )
  )

defoliation_change_summary <- defoliation_change_data |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_records = n(),
    n_trees = n_distinct(BANR),
    mean_change = mean(annual_change, na.rm = TRUE),
    sd_change = sd(annual_change, na.rm = TRUE),
    se_change = sd_change / sqrt(n_records),
    median_change = median(annual_change, na.rm = TRUE),
    q25_change = unname(quantile(annual_change, probs = 0.25, na.rm = TRUE)),
    q75_change = unname(quantile(annual_change, probs = 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  arrange(YEAR, analysis_Treatment)

defoliation_change_plot <- ggplot(
  defoliation_change_summary,
  aes(x = YEAR, y = mean_change, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_change - se_change, ymax = mean_change + se_change),
    width = 0.25,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(defoliation_change_start_year, defoliation_change_end_year, by = 1)
  ) +
  labs(
    x = "",
    y = "Annual defoliation change (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

mortality_summary <- crown_trans_spatial |>
  filter(
    spatial_match,
    analysis_include,
    YEAR >= mortality_start_year,
    YEAR <= mortality_end_year
  ) |>
  select(
    YEAR,
    BANR,
    analysis_Treatment,
    all_of(mortality_column),
    CUT
  ) |>
  rename(mortality_flag = all_of(mortality_column)) |>
  filter(!is.na(mortality_flag)) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = c("Control", "Irrigation", "Irrigation stop")
    ),
    mortality_column = mortality_column,
    mortality_definition = mortality_definition_label,
    mortality_flag = as.numeric(mortality_flag > 0),
    cut_flag = !is.na(CUT) & CUT > 0,
    disappeared_excluded = cut_flag & mortality_flag == 0,
    at_risk = !disappeared_excluded
  ) |>
  group_by(YEAR, analysis_Treatment, mortality_definition, mortality_column) |>
  summarise(
    n_records_with_flag = n(),
    n_trees_with_flag = n_distinct(BANR),
    n_disappeared_excluded = sum(disappeared_excluded),
    n_at_risk = sum(at_risk),
    n_dead = sum(mortality_flag > 0 & at_risk),
    mortality_fraction = if_else(n_at_risk > 0, n_dead / n_at_risk, NA_real_),
    mortality_percent = 100 * mortality_fraction,
    .groups = "drop"
  ) |>
  arrange(YEAR, analysis_Treatment)

mortality_plot <- ggplot(
  mortality_summary,
  aes(x = YEAR, y = mortality_percent, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(breaks = sort(unique(mortality_summary$YEAR))) +
  labs(
    x = "",
    y = "Annual health-related mortality rate (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

combined_year_breaks <- seq(
  defoliation_change_start_year,
  defoliation_change_end_year,
  by = 1
)

mortality_panel_plot <- ggplot(
  mortality_summary,
  aes(x = YEAR, y = mortality_percent, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(breaks = combined_year_breaks) +
  labs(
    x = "",
    y = "Annual\nmortality rate (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

defoliation_change_panel_plot <- ggplot(
  defoliation_change_summary,
  aes(x = YEAR, y = mean_change, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = mean_change - se_change, ymax = mean_change + se_change),
    width = 0.25,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(breaks = combined_year_breaks) +
  labs(
    x = "",
    y = "Annual\ndefoliation change (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

defoliation_mortality_plot <-
  mortality_panel_plot / defoliation_change_panel_plot +
  plot_layout(heights = c(1, 1))

defoliation_change_full_period_summary <- crown_trans_spatial |>
  filter(
    spatial_match,
    analysis_include,
    YEAR >= full_period_start_year,
    YEAR <= full_period_end_year,
    !is.na(crown_transparency_change_mort_cor)
  ) |>
  transmute(
    YEAR,
    BANR,
    analysis_Treatment,
    annual_change = crown_transparency_change_mort_cor
  ) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = c("Control", "Irrigation", "Irrigation stop")
    )
  ) |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_records = n(),
    n_trees = n_distinct(BANR),
    mean_change = mean(annual_change, na.rm = TRUE),
    sd_change = sd(annual_change, na.rm = TRUE),
    se_change = sd_change / sqrt(n_records),
    .groups = "drop"
  ) |>
  arrange(YEAR, analysis_Treatment)

mortality_full_period_summary <- crown_trans_spatial |>
  filter(
    spatial_match,
    analysis_include,
    YEAR >= full_period_start_year,
    YEAR <= full_period_end_year
  ) |>
  select(
    YEAR,
    BANR,
    analysis_Treatment,
    all_of(mortality_column),
    CUT
  ) |>
  rename(mortality_flag = all_of(mortality_column)) |>
  filter(!is.na(mortality_flag)) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = c("Control", "Irrigation", "Irrigation stop")
    ),
    mortality_flag = as.numeric(mortality_flag > 0),
    cut_flag = !is.na(CUT) & CUT > 0,
    disappeared_excluded = cut_flag & mortality_flag == 0,
    at_risk = !disappeared_excluded
  ) |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_at_risk = sum(at_risk),
    n_dead = sum(mortality_flag > 0 & at_risk),
    mortality_percent = if_else(n_at_risk > 0, 100 * n_dead / n_at_risk, NA_real_),
    .groups = "drop"
  ) |>
  arrange(YEAR, analysis_Treatment)

full_period_year_breaks <- seq(full_period_start_year, full_period_end_year, by = 1)

mortality_full_period_panel_plot <- ggplot(
  mortality_full_period_summary,
  aes(x = YEAR, y = mortality_percent, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_vline(
    xintercept = irrigation_stop_year,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = full_period_year_breaks,
    limits = c(full_period_start_year, full_period_end_year)
  ) +
  labs(
    x = "",
    y = "Annual\nmortality rate (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

defoliation_change_full_period_panel_plot <- ggplot(
  defoliation_change_full_period_summary,
  aes(x = YEAR, y = mean_change, color = analysis_Treatment, group = analysis_Treatment)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
  geom_vline(
    xintercept = irrigation_stop_year,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = mean_change - se_change, ymax = mean_change + se_change),
    width = 0.25,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = full_period_year_breaks,
    limits = c(full_period_start_year, full_period_end_year)
  ) +
  labs(
    x = "",
    y = "Annual\ndefoliation change (%)",
    color = "Treatment"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

defoliation_mortality_full_period_plot <-
  mortality_full_period_panel_plot / defoliation_change_full_period_panel_plot +
  plot_layout(heights = c(1, 1))

lai_lookup <- lai |>
  transmute(
    YEAR = year,
    analysis_Plot = plot,
    LAI_block = block,
    LAI_treatment_raw = treatment,
    LAI_treatment = case_when(
      tolower(treatment) == "control" ~ "Control",
      tolower(treatment) == "irrigated" ~ "Irrigation",
      .default = treatment
    ),
    LAI = LAI
  )

crown_plot_summary <- crown_trans_spatial |>
  filter(spatial_match, analysis_include, !is.na(crown_transparency)) |>
  group_by(YEAR, analysis_Treatment, analysis_Plot) |>
  summarise(
    n_trees = n_distinct(BANR),
    n_records = n(),
    plot_mean_crown_transparency = mean(crown_transparency, na.rm = TRUE),
    plot_sd_crown_transparency = sd(crown_transparency, na.rm = TRUE),
    .groups = "drop"
  )

lai_crown_plot_data <- crown_plot_summary |>
  inner_join(lai_lookup, by = c("YEAR", "analysis_Plot")) |>
  filter(
    analysis_Treatment %in% lai_comparison_treatments,
    analysis_Treatment == LAI_treatment
  ) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = lai_comparison_treatments
    )
  ) |>
  arrange(YEAR, analysis_Treatment, analysis_Plot)

lai_timeseries_summary <- lai_crown_plot_data |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_LAI_plots = n(),
    mean_LAI = mean(LAI, na.rm = TRUE),
    sd_LAI = sd(LAI, na.rm = TRUE),
    se_LAI = sd_LAI / sqrt(n_LAI_plots),
    .groups = "drop"
  )

crown_timeseries_summary <- crown_plot_summary |>
  filter(
    YEAR >= crown_timeseries_start_year,
    YEAR <= crown_timeseries_end_year,
    analysis_Treatment %in% crown_timeseries_treatments
  ) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = crown_timeseries_treatments
    )
  ) |>
  group_by(YEAR, analysis_Treatment) |>
  summarise(
    n_crown_plots = n(),
    n_trees = sum(n_trees),
    mean_crown_transparency = mean(plot_mean_crown_transparency, na.rm = TRUE),
    sd_crown_transparency = sd(plot_mean_crown_transparency, na.rm = TRUE),
    se_crown_transparency = sd_crown_transparency / sqrt(n_crown_plots),
    .groups = "drop"
  )

lai_crown_summary <- full_join(
  crown_timeseries_summary,
  lai_timeseries_summary,
  by = c("YEAR", "analysis_Treatment")
) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = crown_timeseries_treatments
    )
  ) |>
  arrange(YEAR, analysis_Treatment)

lai_axis_range <- expand_range(c(
  lai_crown_summary$mean_LAI - lai_crown_summary$se_LAI,
  lai_crown_summary$mean_LAI + lai_crown_summary$se_LAI
))
crown_transparency_axis_range <- expand_range(c(
  crown_timeseries_summary$mean_crown_transparency - crown_timeseries_summary$se_crown_transparency,
  crown_timeseries_summary$mean_crown_transparency + crown_timeseries_summary$se_crown_transparency
))

scale_crown_transparency_to_lai <- function(crown_transparency) {
  (crown_transparency_axis_range[2] - crown_transparency) /
    diff(crown_transparency_axis_range) *
    diff(lai_axis_range) +
    lai_axis_range[1]
}

scale_lai_to_crown_transparency <- function(lai_value) {
  crown_transparency_axis_range[2] -
    (lai_value - lai_axis_range[1]) /
    diff(lai_axis_range) *
    diff(crown_transparency_axis_range)
}

lai_crown_timeseries <- bind_rows(
  lai_timeseries_summary |>
    transmute(
      YEAR,
      analysis_Treatment,
      metric = "LAI",
      mean_value = mean_LAI,
      se_value = se_LAI,
      plot_value = mean_LAI,
      plot_lower = mean_LAI - se_LAI,
      plot_upper = mean_LAI + se_LAI
    ),
  crown_timeseries_summary |>
    transmute(
      YEAR,
      analysis_Treatment,
      metric = "Crown transparency",
      mean_value = mean_crown_transparency,
      se_value = se_crown_transparency,
      plot_value = scale_crown_transparency_to_lai(mean_value),
      plot_lower = pmin(
        scale_crown_transparency_to_lai(mean_value - se_value),
        scale_crown_transparency_to_lai(mean_value + se_value)
      ),
      plot_upper = pmax(
        scale_crown_transparency_to_lai(mean_value - se_value),
        scale_crown_transparency_to_lai(mean_value + se_value)
      )
    )
) |>
  mutate(
    analysis_Treatment = factor(
      analysis_Treatment,
      levels = crown_timeseries_treatments
    ),
    metric = factor(metric, levels = c("LAI", "Crown transparency"))
  )

correlation_for_group <- function(data) {
  data <- data |>
    filter(!is.na(LAI), !is.na(plot_mean_crown_transparency))

  if (nrow(data) < 3) {
    return(tibble(
      n = nrow(data),
      pearson_r = NA_real_,
      p_value = NA_real_,
      intercept = NA_real_,
      slope = NA_real_
    ))
  }

  test <- cor.test(data$LAI, data$plot_mean_crown_transparency)
  fit <- lm(plot_mean_crown_transparency ~ LAI, data = data)

  tibble(
    n = nrow(data),
    pearson_r = unname(test$estimate),
    p_value = test$p.value,
    intercept = unname(coef(fit)[1]),
    slope = unname(coef(fit)[2])
  )
}

lai_crown_correlation <- lai_crown_plot_data |>
  group_by(analysis_Treatment) |>
  group_modify(~correlation_for_group(.x)) |>
  ungroup() |>
  mutate(
    label = paste0(
      "r = ", round(pearson_r, 2),
      "\np = ", if_else(p_value < 0.001, "<0.001", format(signif(p_value, 2), scientific = FALSE))
    )
  )

lai_crown_timeseries_plot <- ggplot(
  lai_crown_timeseries,
  aes(
    x = YEAR,
    y = plot_value,
    color = analysis_Treatment,
    linetype = metric,
    group = interaction(analysis_Treatment, metric)
  )
) +
  geom_vline(
    data = tibble(
      analysis_Treatment = factor(
        "Irrigation stop",
        levels = crown_timeseries_treatments
      ),
      YEAR = irrigation_stop_year
    ),
    aes(xintercept = YEAR),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = plot_lower, ymax = plot_upper),
    width = 0.25,
    alpha = 0.75,
    linewidth = 0.45
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~analysis_Treatment, ncol = 1, scales = "free_y") +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_linetype_manual(values = c("LAI" = "solid", "Crown transparency" = "dashed")) +
  scale_x_continuous(
    breaks = seq(crown_timeseries_start_year, crown_timeseries_end_year, by = 1),
    limits = c(crown_timeseries_start_year, crown_timeseries_end_year)
  ) +
  scale_y_continuous(
    name = "LAI",
    sec.axis = sec_axis(
      ~scale_lai_to_crown_transparency(.),
      name = "Crown transparency (%)"
    )
  ) +
  labs(
    x = "",
    color = "Treatment",
    linetype = ""
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

lai_crown_scatter_plot <- ggplot(
  lai_crown_plot_data,
  aes(x = LAI, y = plot_mean_crown_transparency, color = analysis_Treatment)
) +
  geom_point(aes(shape = factor(analysis_Plot)), size = 2.4, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  geom_text(
    data = lai_crown_correlation,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    size = 3.5
  ) +
  facet_wrap(~analysis_Treatment, nrow = 1) +
  scale_color_manual(values = scenario_colors, drop = FALSE) +
  scale_shape_manual(values = c(
    "1" = 16,
    "2" = 17,
    "3" = 15,
    "4" = 3,
    "5" = 7,
    "6" = 8,
    "7" = 0,
    "8" = 2
  )) +
  labs(
    x = "LAI",
    y = "Mean crown transparency (%)",
    color = "Treatment",
    shape = "Plot"
  ) +
  theme_bw() +
  pfyn_text_theme +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

crown_trans_spatial_sf <- crown_trans_spatial |>
  filter(spatial_match) |>
  st_as_sf(wkt = "shp_geometry_wkt", crs = st_crs(tree_spatial), remove = FALSE)

write_csv(crown_trans_spatial, spatial_csv_path)
write_csv(unmatched_banr, unmatched_banr_path)
write_csv(unmatched_xy, unmatched_xy_path)
write_csv(summary_by_year_plot, summary_path)
write_csv(crown_summary, crown_summary_path)
write_csv(defoliation_change_summary, defoliation_change_summary_path)
write_csv(mortality_summary, mortality_summary_path)
write_csv(lai_crown_plot_data, lai_crown_plot_data_path)
write_csv(lai_crown_summary, lai_crown_summary_path)
write_csv(lai_crown_correlation, lai_crown_correlation_path)
saveRDS(crown_trans_spatial_sf, spatial_rds_path)
ggsave(crown_plot_path, crown_transparency_plot, width = 9, height = 5, dpi = 300)
ggsave(crown_plot_se_path, crown_transparency_se_plot, width = 9, height = 5, dpi = 300)
ggsave(defoliation_change_plot_path, defoliation_change_plot, width = 10, height = 6, dpi = 300)
ggsave(mortality_plot_path, mortality_plot, width = 9, height = 4, dpi = 300)
ggsave(defoliation_mortality_plot_path, defoliation_mortality_plot, width = 10, height = 6, dpi = 300)
ggsave(
  defoliation_mortality_full_period_plot_path,
  defoliation_mortality_full_period_plot,
  width = 12,
  height = 6,
  dpi = 300
)
ggsave(lai_crown_timeseries_path, lai_crown_timeseries_plot, width = 10, height = 10, dpi = 300)
ggsave(lai_crown_scatter_path, lai_crown_scatter_plot, width = 10, height = 4, dpi = 300)

message("Rows in crown condition table: ", raw_crown_condition_rows)
message("Rows after primary inventory filter: ", nrow(crown_trans))
message("Non-primary rows excluded: ", non_primary_crown_condition_rows)
message("Rows matched to shapefile: ", sum(crown_trans_spatial$spatial_match))
message("Rows without a shapefile match: ", sum(!crown_trans_spatial$spatial_match))
message("Rows included in three-treatment analysis: ", sum(crown_trans_spatial$analysis_include))
message("Rows excluded from three-treatment analysis: ", sum(!crown_trans_spatial$analysis_include))
message("Distinct unmatched BANR values: ", nrow(unmatched_banr))
message("Distinct non-exact coordinate pairs: ", nrow(unmatched_xy))
message("Wrote enriched CSV: ", spatial_csv_path)
message("Wrote matched sf RDS: ", spatial_rds_path)
message("Wrote unmatched BANR report: ", unmatched_banr_path)
message("Wrote coordinate diagnostic report: ", unmatched_xy_path)
message("Wrote summary table: ", summary_path)
message("Wrote crown transparency summary: ", crown_summary_path)
message("Wrote crown transparency plot: ", crown_plot_path)
message("Wrote crown transparency SE plot: ", crown_plot_se_path)
message("Wrote defoliation change summary: ", defoliation_change_summary_path)
message("Wrote defoliation change plot: ", defoliation_change_plot_path)
message("Wrote mortality summary: ", mortality_summary_path)
message("Wrote mortality plot: ", mortality_plot_path)
message("Wrote defoliation change and mortality plot: ", defoliation_mortality_plot_path)
message(
  "Wrote full-period defoliation change and mortality plot: ",
  defoliation_mortality_full_period_plot_path
)
message("Wrote LAI/crown plot-level data: ", lai_crown_plot_data_path)
message("Wrote LAI/crown treatment-year summary: ", lai_crown_summary_path)
message("Wrote LAI/crown correlation table: ", lai_crown_correlation_path)
message("Wrote LAI/crown time-series plot: ", lai_crown_timeseries_path)
message("Wrote LAI/crown correlation plot: ", lai_crown_scatter_path)
