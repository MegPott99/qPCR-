# qPCR Time Series Analysis Script
# For analyzing bacterial populations over time in gut model experiments
#
# Sample format: Donor_Vessel_Day (e.g., LC01_V1_D3) or Donor_Slurry (baseline)
# Donors: LC01-04 (Long COVID), Rec01-03 (Recovered), POOL
# Vessels: V1 (proximal colon), V3 (distal colon) - analyzed separately
# Days: Slurry=D0, D1, D2, D3, D6, D8, D10, D13, D17, D20

# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(broom)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Standard curve settings (fixed for this experiment)
STANDARD_START_CONC <- 5e9
DILUTION_FACTOR <- 10

# Define donor groups for analysis
DONOR_GROUPS <- list(

Long_COVID = c("LC01", "LC02", "LC03", "LC04"),
  Recovered = c("Rec01", "Rec02", "Rec03"),
  Pool = c("POOL")
)

# Day order for proper time series plotting
DAY_ORDER <- c("Slurry", "D1", "D2", "D3", "D6", "D8", "D10", "D13", "D17", "D20")

# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

#' Parse sample IDs in Donor_Vessel_Day or Donor_Slurry format
#'
#' @param sample_ids Vector of sample ID strings
#' @return Data frame with parsed components
parse_sample_ids <- function(sample_ids) {

  parsed <- data.frame(
    Original_Sample_ID = sample_ids,
    Donor = NA_character_,
    Vessel = NA_character_,
    Day = NA_character_,
    Day_Numeric = NA_real_,
    Donor_Group = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(sample_ids)) {
    sample <- sample_ids[i]

    if (is.na(sample) || sample == "") next

    parts <- str_split(sample, "_")[[1]]

    # Check for Slurry samples (Donor_Slurry format)
    if (length(parts) == 2 && tolower(parts[2]) == "slurry") {
      parsed$Donor[i] <- parts[1]
      parsed$Vessel[i] <- "Slurry"
      parsed$Day[i] <- "Slurry"
      parsed$Day_Numeric[i] <- 0

    # Check for standard format (Donor_Vessel_Day)
    } else if (length(parts) == 3) {
      parsed$Donor[i] <- parts[1]
      parsed$Vessel[i] <- parts[2]
      parsed$Day[i] <- parts[3]

      # Extract numeric day value
      day_num <- as.numeric(str_extract(parts[3], "\\d+"))
      parsed$Day_Numeric[i] <- ifelse(is.na(day_num), NA, day_num)
    }

    # Assign donor group
    donor <- parsed$Donor[i]
    if (!is.na(donor)) {
      if (donor %in% DONOR_GROUPS$Long_COVID) {
        parsed$Donor_Group[i] <- "Long COVID"
      } else if (donor %in% DONOR_GROUPS$Recovered) {
        parsed$Donor_Group[i] <- "Recovered"
      } else if (donor %in% DONOR_GROUPS$Pool) {
        parsed$Donor_Group[i] <- "Pool"
      } else {
        parsed$Donor_Group[i] <- "Unknown"
      }
    }
  }

  # Convert Day to ordered factor for proper plotting
  parsed$Day <- factor(parsed$Day, levels = DAY_ORDER, ordered = TRUE)

  return(parsed)
}

#' Extract bacterial target from filename
#'
#' @param filename The CSV filename (format: date_BacterialTarget_plate_output.csv)
#' @return The bacterial target name
extract_bacterial_target <- function(filename) {
  # Remove path and extension
  base_name <- tools::file_path_sans_ext(basename(filename))

  # Split by underscore
  parts <- str_split(base_name, "_")[[1]]

  # Target should be the second element (after date)
  if (length(parts) >= 2) {
    return(parts[2])
  } else {
    return(base_name)
  }
}

# =============================================================================
# STANDARD CURVE FUNCTIONS
# =============================================================================

#' Process standard curve data and fit linear model
#'
#' @param data Raw data frame containing standards
#' @param start_conc Starting concentration (default 5e9)
#' @param dilution_factor Dilution factor between standards (default 10)
#' @return List containing model, slope, intercept, R², and efficiency
process_standard_curve <- function(data, start_conc = STANDARD_START_CONC,
                                    dilution_factor = DILUTION_FACTOR) {

  # Filter for standards
  standard_data <- data %>%
    filter(grepl("Standard", `Sample name`, ignore.case = TRUE)) %>%
    filter(!is.na(Ct) & Ct != "No Ct" & Ct != "" & Ct != "Undetermined") %>%
    mutate(Ct = as.numeric(as.character(Ct))) %>%
    filter(!is.na(Ct) & !is.infinite(Ct))

  if (nrow(standard_data) == 0) {
    stop("No valid standard data found")
  }

  # Extract standard numbers and assign concentrations
  standard_data <- standard_data %>%
    mutate(
      Standard_Number = as.numeric(str_extract(`Sample name`, "\\d+")),
      Concentration = start_conc / (dilution_factor ^ (Standard_Number - 1)),
      Log_Concentration = log10(Concentration)
    ) %>%
    filter(!is.na(Standard_Number))

  # Calculate mean Ct per standard
  standards_summary <- standard_data %>%
    group_by(Standard_Number, Concentration, Log_Concentration) %>%
    summarise(
      Mean_Ct = mean(Ct, na.rm = TRUE),
      SD_Ct = sd(Ct, na.rm = TRUE),
      N = n(),
      .groups = 'drop'
    )

  # Fit linear model
  model <- lm(Mean_Ct ~ Log_Concentration, data = standards_summary)
  model_summary <- glance(model)
  coefficients <- tidy(model)

  slope <- coefficients$estimate[coefficients$term == "Log_Concentration"]
  intercept <- coefficients$estimate[coefficients$term == "(Intercept)"]
  r_squared <- model_summary$r.squared
  efficiency <- (10^(-1/slope) - 1) * 100

  cat(sprintf("  Standard curve: R² = %.4f, Efficiency = %.1f%%\n", r_squared, efficiency))

  return(list(
    model = model,
    slope = slope,
    intercept = intercept,
    r_squared = r_squared,
    efficiency = efficiency,
    data = standards_summary,
    start_concentration = start_conc
  ))
}

#' Calculate copy numbers from Ct values using standard curve
#'
#' @param ct_values Vector of Ct values
#' @param slope Slope from standard curve
#' @param intercept Intercept from standard curve
#' @return Data frame with Ct, Log_Copy_Number, and Copy_Number
calculate_copy_numbers <- function(ct_values, slope, intercept) {
  log_copy <- (ct_values - intercept) / slope
  copy_number <- 10^log_copy

  return(data.frame(
    Ct = ct_values,
    Log_Copy_Number = log_copy,
    Copy_Number = copy_number
  ))
}

# =============================================================================
# FILE PROCESSING FUNCTIONS
# =============================================================================

#' Read qPCR data file with automatic header detection
#'
#' @param file_path Path to the CSV file
#' @return Data frame with cleaned data
read_qpcr_file <- function(file_path) {

  # Try different encodings
  raw_data <- NULL
  encodings <- c("UTF-8", "latin1", "cp1252", "ISO-8859-1")

  for (enc in encodings) {
    tryCatch({
      raw_data <- read_lines(file_path, locale = locale(encoding = enc))
      break
    }, error = function(e) NULL)
  }

  if (is.null(raw_data)) {
    stop("Could not read file with any supported encoding")
  }

  # Find header line (contains "Well" and "Sample")
  header_line <- NULL
  for (i in seq_along(raw_data)) {
    if (grepl("Well", raw_data[i], ignore.case = TRUE) &&
        grepl("Sample", raw_data[i], ignore.case = TRUE)) {
      header_line <- i
      break
    }
  }

  if (is.null(header_line)) {
    stop("Could not find data table header")
  }

  # Read data from header line
  data <- read_csv(file_path, skip = header_line - 1, show_col_types = FALSE)

  # Clean column names
  data <- data %>%
    select_if(~!all(is.na(.))) %>%
    rename_with(~trimws(.))

  # Standardize column names
  if (!"Sample name" %in% names(data)) {
    sample_col <- names(data)[grepl("sample", names(data), ignore.case = TRUE)][1]
    if (!is.na(sample_col)) names(data)[names(data) == sample_col] <- "Sample name"
  }

  if (!"Ct" %in% names(data)) {
    ct_col <- names(data)[grepl("^Ct$", names(data), ignore.case = TRUE)][1]
    if (!is.na(ct_col)) names(data)[names(data) == ct_col] <- "Ct"
  }

  return(data)
}

#' Process a single qPCR file
#'
#' @param file_path Path to the CSV file
#' @return List containing all results
process_qpcr_file <- function(file_path) {

  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Processing:", basename(file_path), "\n")

  # Extract bacterial target from filename
  bacterial_target <- extract_bacterial_target(file_path)
  cat("  Bacterial target:", bacterial_target, "\n")

  # Read data
  data <- read_qpcr_file(file_path)

  # Identify sample types
  all_samples <- unique(data$`Sample name`[!is.na(data$`Sample name`)])
  standards <- all_samples[grepl("Standard", all_samples, ignore.case = TRUE)]
  ntcs <- all_samples[grepl("NTC|No Template|Negative", all_samples, ignore.case = TRUE)]
  experimental <- all_samples[!all_samples %in% c(standards, ntcs)]

  cat(sprintf("  Found: %d standards, %d NTCs, %d experimental samples\n",
              length(standards), length(ntcs), length(experimental)))

  # Process standard curve
  curve <- process_standard_curve(data)

  # Process experimental samples
  if (length(experimental) == 0) {
    cat("  No experimental samples found\n")
    return(list(
      file_name = basename(file_path),
      bacterial_target = bacterial_target,
      standard_curve = curve,
      samples = NULL,
      summary = NULL
    ))
  }

  # Filter and clean experimental data
  sample_data <- data %>%
    filter(`Sample name` %in% experimental) %>%
    filter(!is.na(Ct) & Ct != "No Ct" & Ct != "" & Ct != "Undetermined") %>%
    mutate(Ct = as.numeric(as.character(Ct))) %>%
    filter(!is.na(Ct) & !is.infinite(Ct))

  # Calculate copy numbers for each Ct value
  copy_results <- calculate_copy_numbers(sample_data$Ct, curve$slope, curve$intercept)

  # Parse sample IDs
  parsed <- parse_sample_ids(sample_data$`Sample name`)

  # Combine results
  sample_results <- sample_data %>%
    select(Well, `Sample name`, Ct) %>%
    bind_cols(copy_results %>% select(-Ct)) %>%
    bind_cols(parsed %>% select(-Original_Sample_ID)) %>%
    mutate(
      Bacterial_Target = bacterial_target,
      Source_File = basename(file_path)
    )

  # Calculate summary statistics (average replicates)
  summary_stats <- sample_results %>%
    group_by(`Sample name`, Donor, Vessel, Day, Day_Numeric, Donor_Group, Bacterial_Target) %>%
    summarise(
      N_Replicates = n(),
      Mean_Ct = mean(Ct, na.rm = TRUE),
      SD_Ct = sd(Ct, na.rm = TRUE),
      Mean_Log_Copy = mean(Log_Copy_Number, na.rm = TRUE),
      SD_Log_Copy = sd(Log_Copy_Number, na.rm = TRUE),
      Mean_Copy_Number = mean(Copy_Number, na.rm = TRUE),
      SD_Copy_Number = sd(Copy_Number, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(Source_File = basename(file_path))

  cat(sprintf("  Quantified: %d unique samples from %d Ct values\n",
              nrow(summary_stats), nrow(sample_results)))

  return(list(
    file_name = basename(file_path),
    bacterial_target = bacterial_target,
    standard_curve = curve,
    individual_results = sample_results,
    summary = summary_stats
  ))
}

# =============================================================================
# BATCH PROCESSING
# =============================================================================

#' Process multiple qPCR files
#'
#' @param file_paths Vector of file paths (or NULL to use file dialog)
#' @return List of all results
process_batch <- function(file_paths = NULL) {

  # If no files provided, use file dialog

if (is.null(file_paths)) {
    cat("Select your qPCR CSV files. Click 'Cancel' when done.\n\n")
    file_paths <- c()

    repeat {
      selected <- tryCatch(file.choose(), error = function(e) NULL)
      if (is.null(selected)) break
      if (grepl("\\.csv$", selected, ignore.case = TRUE)) {
        file_paths <- c(file_paths, selected)
        cat("Added:", basename(selected), "\n")
      }
    }
  }

  if (length(file_paths) == 0) {
    stop("No files selected")
  }

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PROCESSING", length(file_paths), "FILES\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  # Process each file
  all_results <- list()
  failed <- c()

  for (i in seq_along(file_paths)) {
    tryCatch({
      result <- process_qpcr_file(file_paths[i])
      all_results[[i]] <- result
    }, error = function(e) {
      cat("  ERROR:", e$message, "\n")
      failed <- c(failed, basename(file_paths[i]))
      all_results[[i]] <- NULL
    })
  }

  # Remove failed files
  all_results <- all_results[!sapply(all_results, is.null)]

  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("BATCH COMPLETE\n")
  cat(sprintf("  Successful: %d files\n", length(all_results)))
  if (length(failed) > 0) {
    cat(sprintf("  Failed: %d files (%s)\n", length(failed), paste(failed, collapse = ", ")))
  }

  return(all_results)
}

#' Combine results from multiple files into master data frames
#'
#' @param results_list List of results from process_batch()
#' @return List with combined individual and summary data
combine_results <- function(results_list) {

  # Combine all individual results
  all_individual <- bind_rows(lapply(results_list, function(x) x$individual_results))

  # Combine all summaries
  all_summary <- bind_rows(lapply(results_list, function(x) x$summary))

  # Create curve info summary
  curve_info <- data.frame(
    File = sapply(results_list, function(x) x$file_name),
    Bacterial_Target = sapply(results_list, function(x) x$bacterial_target),
    R_Squared = sapply(results_list, function(x) x$standard_curve$r_squared),
    Efficiency = sapply(results_list, function(x) x$standard_curve$efficiency),
    Slope = sapply(results_list, function(x) x$standard_curve$slope),
    Intercept = sapply(results_list, function(x) x$standard_curve$intercept)
  )

  return(list(
    individual = all_individual,
    summary = all_summary,
    curves = curve_info
  ))
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

#' Create time series line plot for a single bacterial target
#'
#' @param data Summary data frame (from combine_results()$summary)
#' @param target Bacterial target to plot
#' @param vessel Which vessel to plot ("V1", "V3", or "both")
#' @param log_scale Use log10 scale for y-axis (default TRUE)
#' @param include_slurry Include slurry samples as Day 0 (default TRUE)
#' @return ggplot object
plot_timeseries <- function(data, target, vessel = "both",
                            log_scale = TRUE, include_slurry = TRUE) {

  # Filter for target
  plot_data <- data %>%
    filter(Bacterial_Target == target)

  if (!include_slurry) {
    plot_data <- plot_data %>% filter(Vessel != "Slurry")
  }

  # Handle slurry samples - they apply to both vessels
  if (include_slurry && "Slurry" %in% plot_data$Vessel) {
    slurry_data <- plot_data %>% filter(Vessel == "Slurry")

    # Duplicate slurry for both V1 and V3
    slurry_v1 <- slurry_data %>% mutate(Vessel = "V1")
    slurry_v3 <- slurry_data %>% mutate(Vessel = "V3")

    plot_data <- plot_data %>%
      filter(Vessel != "Slurry") %>%
      bind_rows(slurry_v1, slurry_v3)
  }

  # Filter for vessel
  if (vessel != "both") {
    plot_data <- plot_data %>% filter(Vessel == vessel)
  }

  # Create plot
  p <- ggplot(plot_data, aes(x = Day_Numeric, y = Mean_Copy_Number,
                              color = Donor, group = Donor)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Mean_Copy_Number - SD_Copy_Number,
                      ymax = Mean_Copy_Number + SD_Copy_Number),
                  width = 0.3, alpha = 0.5) +
    labs(
      title = paste(target, "- Population Over Time"),
      x = "Day",
      y = "Copy Number",
      color = "Donor"
    ) +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 8, 10, 13, 17, 20)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  # Add log scale if requested
  if (log_scale) {
    p <- p +
      scale_y_log10(labels = scales::scientific) +
      labs(y = "Copy Number (log scale)")
  }

  # Facet by vessel if showing both
  if (vessel == "both") {
    p <- p + facet_wrap(~Vessel, ncol = 2)
  }

  return(p)
}

#' Create time series plot colored by donor group
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot
#' @param vessel Which vessel to plot ("V1", "V3", or "both")
#' @param show_individual_donors Show individual donor lines (default TRUE)
#' @return ggplot object
plot_timeseries_by_group <- function(data, target, vessel = "both",
                                      show_individual_donors = TRUE) {

  plot_data <- data %>%
    filter(Bacterial_Target == target)

  # Handle slurry samples
  if ("Slurry" %in% plot_data$Vessel) {
    slurry_data <- plot_data %>% filter(Vessel == "Slurry")
    slurry_v1 <- slurry_data %>% mutate(Vessel = "V1")
    slurry_v3 <- slurry_data %>% mutate(Vessel = "V3")

    plot_data <- plot_data %>%
      filter(Vessel != "Slurry") %>%
      bind_rows(slurry_v1, slurry_v3)
  }

  if (vessel != "both") {
    plot_data <- plot_data %>% filter(Vessel == vessel)
  }

  # Define colors for donor groups
  group_colors <- c(
    "Long COVID" = "#E41A1C",
    "Recovered" = "#377EB8",
    "Pool" = "#4DAF4A"
  )

  p <- ggplot(plot_data, aes(x = Day_Numeric, y = Mean_Copy_Number))

  if (show_individual_donors) {
    p <- p +
      geom_line(aes(color = Donor_Group, group = Donor),
                linewidth = 0.6, alpha = 0.7) +
      geom_point(aes(color = Donor_Group), size = 2, alpha = 0.7)
  }

  # Add group means
  group_means <- plot_data %>%
    group_by(Donor_Group, Vessel, Day_Numeric) %>%
    summarise(
      Group_Mean = mean(Mean_Copy_Number, na.rm = TRUE),
      Group_SD = sd(Mean_Copy_Number, na.rm = TRUE),
      .groups = 'drop'
    )

  p <- p +
    geom_line(data = group_means,
              aes(y = Group_Mean, color = Donor_Group, group = Donor_Group),
              linewidth = 1.5) +
    geom_point(data = group_means,
               aes(y = Group_Mean, color = Donor_Group),
               size = 3) +
    scale_color_manual(values = group_colors) +
    scale_y_log10(labels = scales::scientific) +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 8, 10, 13, 17, 20)) +
    labs(
      title = paste(target, "- By Donor Group"),
      x = "Day",
      y = "Copy Number (log scale)",
      color = "Group"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  if (vessel == "both") {
    p <- p + facet_wrap(~Vessel, ncol = 2)
  }

  return(p)
}

#' Create multi-panel plot for all bacterial targets
#'
#' @param data Summary data frame
#' @param vessel Which vessel to plot
#' @return ggplot object
plot_all_targets <- function(data, vessel = "V1") {

  plot_data <- data

  # Handle slurry
  if ("Slurry" %in% plot_data$Vessel) {
    slurry_data <- plot_data %>% filter(Vessel == "Slurry")
    slurry_vessel <- slurry_data %>% mutate(Vessel = vessel)

    plot_data <- plot_data %>%
      filter(Vessel != "Slurry") %>%
      bind_rows(slurry_vessel)
  }

  plot_data <- plot_data %>% filter(Vessel == vessel)

  group_colors <- c(
    "Long COVID" = "#E41A1C",
    "Recovered" = "#377EB8",
    "Pool" = "#4DAF4A"
  )

  p <- ggplot(plot_data, aes(x = Day_Numeric, y = Mean_Copy_Number,
                              color = Donor_Group, group = Donor)) +
    geom_line(linewidth = 0.5, alpha = 0.6) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = group_colors) +
    scale_y_log10(labels = scales::scientific) +
    scale_x_continuous(breaks = c(0, 3, 6, 10, 13, 17, 20)) +
    facet_wrap(~Bacterial_Target, scales = "free_y", ncol = 3) +
    labs(
      title = paste("All Bacterial Targets -", vessel),
      x = "Day",
      y = "Copy Number (log scale)",
      color = "Group"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  return(p)
}

# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

#' Save results to CSV files
#'
#' @param combined_results Output from combine_results()
#' @param output_dir Directory to save files (created if doesn't exist)
save_results <- function(combined_results, output_dir = "analysis_results") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save individual results
  individual_file <- file.path(output_dir,
                                paste0("individual_ct_values_", timestamp, ".csv"))
  write_csv(combined_results$individual, individual_file)
  cat("Saved:", individual_file, "\n")

  # Save summary results
  summary_file <- file.path(output_dir,
                            paste0("summary_by_sample_", timestamp, ".csv"))
  write_csv(combined_results$summary, summary_file)
  cat("Saved:", summary_file, "\n")

  # Save curve info
  curves_file <- file.path(output_dir,
                           paste0("standard_curves_", timestamp, ".csv"))
  write_csv(combined_results$curves, curves_file)
  cat("Saved:", curves_file, "\n")

  # Create wide format for easy viewing
  wide_data <- combined_results$summary %>%
    select(Bacterial_Target, Donor, Vessel, Day, Mean_Copy_Number) %>%
    pivot_wider(
      names_from = c(Vessel, Day),
      values_from = Mean_Copy_Number,
      names_sep = "_"
    )

  wide_file <- file.path(output_dir,
                         paste0("wide_format_", timestamp, ".csv"))
  write_csv(wide_data, wide_file)
  cat("Saved:", wide_file, "\n")

  return(list(
    individual = individual_file,
    summary = summary_file,
    curves = curves_file,
    wide = wide_file
  ))
}

#' Save all plots to PDF
#'
#' @param combined_results Output from combine_results()
#' @param output_dir Directory to save plots
save_plots <- function(combined_results, output_dir = "analysis_results") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  targets <- unique(combined_results$summary$Bacterial_Target)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save individual target plots
  for (target in targets) {
    # V1 plot
    p_v1 <- plot_timeseries(combined_results$summary, target, vessel = "V1")
    ggsave(file.path(output_dir, paste0(target, "_V1_", timestamp, ".pdf")),
           p_v1, width = 10, height = 6)

    # V3 plot
    p_v3 <- plot_timeseries(combined_results$summary, target, vessel = "V3")
    ggsave(file.path(output_dir, paste0(target, "_V3_", timestamp, ".pdf")),
           p_v3, width = 10, height = 6)

    # Combined plot
    p_both <- plot_timeseries(combined_results$summary, target, vessel = "both")
    ggsave(file.path(output_dir, paste0(target, "_both_vessels_", timestamp, ".pdf")),
           p_both, width = 12, height = 6)

    # By group plot
    p_group <- plot_timeseries_by_group(combined_results$summary, target, vessel = "both")
    ggsave(file.path(output_dir, paste0(target, "_by_group_", timestamp, ".pdf")),
           p_group, width = 12, height = 6)
  }

  # Save multi-panel overview
  p_all_v1 <- plot_all_targets(combined_results$summary, "V1")
  ggsave(file.path(output_dir, paste0("ALL_targets_V1_", timestamp, ".pdf")),
         p_all_v1, width = 14, height = 10)

  p_all_v3 <- plot_all_targets(combined_results$summary, "V3")
  ggsave(file.path(output_dir, paste0("ALL_targets_V3_", timestamp, ".pdf")),
         p_all_v3, width = 14, height = 10)

  cat("Saved all plots to:", output_dir, "\n")
}

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

#' Run complete analysis pipeline
#'
#' @param file_paths Vector of file paths (NULL to use dialog)
#' @param output_dir Directory for output files
#' @param save_outputs Whether to save CSV and plots
#' @return List with all results
run_analysis <- function(file_paths = NULL, output_dir = "analysis_results",
                         save_outputs = TRUE) {

  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("qPCR TIME SERIES ANALYSIS\n")
  cat("Gut Model Bacterial Population Dynamics\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  # Process files
  results <- process_batch(file_paths)

  # Combine results
  combined <- combine_results(results)

  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  cat("\nBacterial targets processed:\n")
  for (target in unique(combined$summary$Bacterial_Target)) {
    n_samples <- sum(combined$summary$Bacterial_Target == target)
    cat(sprintf("  - %s: %d samples\n", target, n_samples))
  }

  cat("\nDonors found:\n")
  for (group in names(DONOR_GROUPS)) {
    donors_found <- intersect(DONOR_GROUPS[[group]], unique(combined$summary$Donor))
    if (length(donors_found) > 0) {
      cat(sprintf("  %s: %s\n", group, paste(donors_found, collapse = ", ")))
    }
  }

  cat("\nVessels:", paste(unique(combined$summary$Vessel), collapse = ", "), "\n")
  cat("Days:", paste(sort(unique(combined$summary$Day_Numeric)), collapse = ", "), "\n")

  # Save outputs
  if (save_outputs) {
    cat("\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    cat("SAVING RESULTS\n")
    cat(paste(rep("=", 60), collapse = ""), "\n\n")

    saved_files <- save_results(combined, output_dir)
    save_plots(combined, output_dir)
  }

  cat("\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  return(list(
    raw_results = results,
    combined = combined
  ))
}

# =============================================================================
# QUICK START
# =============================================================================

# To run the analysis:
#
# Option 1: Interactive file selection
#   results <- run_analysis()
#
# Option 2: Specify files directly
#   files <- c("path/to/file1.csv", "path/to/file2.csv")
#   results <- run_analysis(files)
#
# After running, you can create additional plots:
#   plot_timeseries(results$combined$summary, "Bacteroides", vessel = "V1")
#   plot_timeseries_by_group(results$combined$summary, "Lactobacillus", vessel = "both")
#   plot_all_targets(results$combined$summary, "V1")
