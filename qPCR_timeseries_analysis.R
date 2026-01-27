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

    # Split by underscore and remove empty parts (handles double underscores)
    parts <- str_split(sample, "_")[[1]]
    parts <- parts[parts != ""]  # Remove empty strings from double underscores

    # Check for Slurry samples - look for "slurry" anywhere in the parts
    is_slurry <- any(tolower(parts) == "slurry")

    if (is_slurry) {
      # Slurry sample: first part is donor
      parsed$Donor[i] <- parts[1]
      parsed$Vessel[i] <- "Slurry"
      parsed$Day[i] <- "Slurry"
      parsed$Day_Numeric[i] <- 0

    # Check for standard format (Donor_Vessel_Day)
    } else if (length(parts) >= 3) {
      parsed$Donor[i] <- parts[1]
      parsed$Vessel[i] <- parts[2]
      parsed$Day[i] <- parts[3]

      # Extract numeric day value
      day_num <- as.numeric(str_extract(parts[3], "\\d+"))
      parsed$Day_Numeric[i] <- ifelse(is.na(day_num), NA, day_num)

    # Handle 2-part format (might be Donor_Day without vessel?)
    } else if (length(parts) == 2) {
      parsed$Donor[i] <- parts[1]
      # Check if second part looks like a day
      if (grepl("^D\\d+$", parts[2], ignore.case = TRUE)) {
        parsed$Day[i] <- parts[2]
        day_num <- as.numeric(str_extract(parts[2], "\\d+"))
        parsed$Day_Numeric[i] <- ifelse(is.na(day_num), NA, day_num)
      } else {
        parsed$Vessel[i] <- parts[2]
      }
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

#' Process standard curve data with Cook's distance outlier detection
#'
#' Fits both original and outlier-removed curves, lets user choose.
#' Cook's distance threshold: 2/n
#'
#' @param data Raw data frame containing standards
#' @param file_name Name of file (for display)
#' @param start_conc Starting concentration (default 5e9)
#' @param dilution_factor Dilution factor between standards (default 10)
#' @return List containing selected model, slope, intercept, R², efficiency
process_standard_curve <- function(data, file_name = "",
                                    start_conc = STANDARD_START_CONC,
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

  # --- FIT ORIGINAL MODEL ---
  model_original <- lm(Mean_Ct ~ Log_Concentration, data = standards_summary)
  summary_original <- glance(model_original)
  coeff_original <- tidy(model_original)

  slope_original <- coeff_original$estimate[coeff_original$term == "Log_Concentration"]
  intercept_original <- coeff_original$estimate[coeff_original$term == "(Intercept)"]
  r2_original <- summary_original$r.squared
  efficiency_original <- (10^(-1/slope_original) - 1) * 100

  # --- COOK'S DISTANCE OUTLIER DETECTION (threshold = 2/n) ---
  n <- nrow(standards_summary)
  cooks_d <- cooks.distance(model_original)
  cooks_threshold <- 2 / n
  outlier_indices <- which(cooks_d > cooks_threshold)

  has_outliers <- length(outlier_indices) > 0

  if (has_outliers) {
    outlier_standards <- standards_summary$Standard_Number[outlier_indices]
    standards_clean <- standards_summary[-outlier_indices, ]

    if (nrow(standards_clean) >= 3) {
      # Fit cleaned model
      model_clean <- lm(Mean_Ct ~ Log_Concentration, data = standards_clean)
      summary_clean <- glance(model_clean)
      coeff_clean <- tidy(model_clean)

      slope_clean <- coeff_clean$estimate[coeff_clean$term == "Log_Concentration"]
      intercept_clean <- coeff_clean$estimate[coeff_clean$term == "(Intercept)"]
      r2_clean <- summary_clean$r.squared
      efficiency_clean <- (10^(-1/slope_clean) - 1) * 100
    } else {
      has_outliers <- FALSE  # Not enough points after removal
    }
  }

  # --- DISPLAY OPTIONS AND LET USER CHOOSE ---
  cat(sprintf("\n  === STANDARD CURVE: %s ===\n", file_name))
  cat(sprintf("  Original curve (%d points):\n", n))
  cat(sprintf("    R² = %.4f | Efficiency = %.1f%% | Slope = %.3f\n",
              r2_original, efficiency_original, slope_original))

  if (has_outliers) {
    cat(sprintf("  Outlier-removed curve (%d points, removed Standard(s) %s):\n",
                nrow(standards_clean), paste(outlier_standards, collapse = ", ")))
    cat(sprintf("    R² = %.4f | Efficiency = %.1f%% | Slope = %.3f\n",
                r2_clean, efficiency_clean, slope_clean))
    cat(sprintf("    Cook's distance threshold: %.4f (2/n)\n", cooks_threshold))

    cat("  Choose curve: 1 = Original, 2 = Outlier-removed: ")
    choice <- readline(prompt = "")

    if (choice == "2") {
      cat("  -> Using outlier-removed curve\n")
      return(list(
        model = model_clean,
        slope = slope_clean,
        intercept = intercept_clean,
        r_squared = r2_clean,
        efficiency = efficiency_clean,
        data = standards_clean,
        start_concentration = start_conc,
        curve_version = "Outlier-removed",
        n_outliers = length(outlier_indices),
        outlier_standards = outlier_standards
      ))
    }
  } else {
    cat("  No outliers detected (Cook's distance, threshold = 2/n)\n")
  }

  cat("  -> Using original curve\n")
  return(list(
    model = model_original,
    slope = slope_original,
    intercept = intercept_original,
    r_squared = r2_original,
    efficiency = efficiency_original,
    data = standards_summary,
    start_concentration = start_conc,
    curve_version = "Original",
    n_outliers = 0,
    outlier_standards = integer(0)
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

  # Process standard curve (with Cook's distance outlier option)
  curve <- process_standard_curve(data, file_name = basename(file_path))

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

#' Report missing data in the dataset
#'
#' @param summary_data Summary data frame
#' @return Invisibly returns missing data summary
report_missing_data <- function(summary_data) {

  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("MISSING DATA CHECK\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # Expected donors, vessels, days
expected_donors <- c("LC01", "LC02", "LC03", "LC04", "Rec01", "Rec02", "Rec03", "POOL")
  expected_vessels <- c("Slurry", "V1", "V3")
  expected_days <- c(0, 1, 2, 3, 6, 8, 10, 13, 17, 20)

  # Get unique values in data
  targets <- unique(summary_data$Bacterial_Target)

  missing_summary <- list()

  for (target in targets) {
    target_data <- summary_data %>% filter(Bacterial_Target == target)

    # Check which donors are present
    donors_present <- unique(target_data$Donor)
    donors_missing <- setdiff(expected_donors, donors_present)

    # Check coverage for each vessel
    for (vessel in c("V1", "V3")) {
      vessel_data <- target_data %>% filter(Vessel == vessel | Vessel == "Slurry")

      if (nrow(vessel_data) == 0) next

      # For each donor, which days are missing?
      for (donor in donors_present) {
        donor_data <- vessel_data %>% filter(Donor == donor)
        days_present <- unique(donor_data$Day_Numeric)
        days_missing <- setdiff(expected_days, days_present)

        if (length(days_missing) > 0) {
          key <- paste(target, vessel, donor, sep = "_")
          missing_summary[[key]] <- list(
            target = target,
            vessel = vessel,
            donor = donor,
            missing_days = days_missing
          )
        }
      }
    }

    # Report donors completely missing for this target
    if (length(donors_missing) > 0) {
      cat("  ", target, ": Missing donors -", paste(donors_missing, collapse = ", "), "\n")
    }
  }

  # Summarize missing timepoints
  if (length(missing_summary) > 0) {
    cat("\n  Missing timepoints detected:\n")

    # Group by target
    for (target in targets) {
      target_missing <- missing_summary[grepl(paste0("^", target, "_"), names(missing_summary))]

      if (length(target_missing) > 0) {
        cat("  ", target, ":\n")

        for (item in target_missing) {
          cat("    ", item$donor, "(", item$vessel, "): Days",
              paste(item$missing_days, collapse = ", "), "\n")
        }
      }
    }

    cat("\n  NOTE: Missing data points are excluded from group averages.\n")
    cat("  Group means at each timepoint only include donors with data.\n")
  } else {
    cat("  No missing timepoints detected for donors present in data.\n")
  }

  cat(paste(rep("-", 60), collapse = ""), "\n\n")

  invisible(missing_summary)
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

  # Report missing data
  report_missing_data(all_summary)

  # Create curve info summary (now includes curve version and outlier info)
  curve_info <- data.frame(
    File = sapply(results_list, function(x) x$file_name),
    Bacterial_Target = sapply(results_list, function(x) x$bacterial_target),
    Curve_Version = sapply(results_list, function(x) {
      v <- x$standard_curve$curve_version
      if (is.null(v)) "Original" else v
    }),
    R_Squared = sapply(results_list, function(x) x$standard_curve$r_squared),
    Efficiency = sapply(results_list, function(x) x$standard_curve$efficiency),
    Slope = sapply(results_list, function(x) x$standard_curve$slope),
    Intercept = sapply(results_list, function(x) x$standard_curve$intercept),
    Outliers_Removed = sapply(results_list, function(x) {
      n <- x$standard_curve$n_outliers
      if (is.null(n)) 0 else n
    })
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

# Custom color palette: Selected from user's palette
# Using distinct colors across the full palette, avoiding beige/grey
DONOR_COLORS <- c(
  # Long COVID donors - varied colors from palette
  "LC01" = "#1D3557",
  "LC02" = "#C9787A",
  "LC03" = "#7B68A6",
  "LC04" = "#E07A5F",
  # Recovered donors - varied colors from palette
  "Rec01" = "#4A90A4",
  "Rec02" = "#81B29A",
  "Rec03" = "#264653",
  # Pool - distinct teal
  "POOL" = "#2A9D8F"
)

# Group colors (for group-level summaries)
GROUP_COLORS <- c(
  "Long COVID" = "#C9787A",
  "Recovered" = "#4A90A4",
  "Pool" = "#2A9D8F"
)

#' Create a consistent log10 y-axis scale
#' Computes shared limits across all data for comparable plots
#'
#' @param data Summary data frame
#' @param padding Padding factor for limits (default 2 = half decade below/above)
#' @return scale_y_log10 object
consistent_y_scale <- function(data, padding = 2) {
  all_values <- data$Mean_Copy_Number[!is.na(data$Mean_Copy_Number) & data$Mean_Copy_Number > 0]

  if (length(all_values) == 0) {
    return(scale_y_log10(labels = scales::scientific))
  }

  y_min <- min(all_values) / padding
  y_max <- max(all_values) * padding

  # Round to nearest power of 10 for clean limits
  y_min <- 10^floor(log10(y_min))
  y_max <- 10^ceiling(log10(y_max))

  # Create clean breaks at each power of 10
  breaks <- 10^seq(log10(y_min), log10(y_max))

  scale_y_log10(
    limits = c(y_min, y_max),
    breaks = breaks,
    labels = scales::scientific
  )
}

#' Prepare plot data - handles slurry and filters vessels
#' @param data Summary data frame
#' @param target Bacterial target (or NULL for all)
#' @param vessel Which vessel ("V1", "V3", or "both")
#' @param include_slurry Include slurry as Day 0
#' @return Filtered data frame
prepare_plot_data <- function(data, target = NULL, vessel = "both", include_slurry = TRUE) {

  plot_data <- data

  # Filter for target if specified

if (!is.null(target)) {
    plot_data <- plot_data %>% filter(Bacterial_Target == target)
  }

  # Remove rows with missing data
  plot_data <- plot_data %>%
    filter(!is.na(Mean_Copy_Number) & !is.na(Day_Numeric) & !is.na(Donor))

  if (nrow(plot_data) == 0) return(plot_data)

  # Handle slurry samples - duplicate for both vessels
  if (include_slurry && "Slurry" %in% plot_data$Vessel) {
    slurry_data <- plot_data %>% filter(Vessel == "Slurry")
    slurry_v1 <- slurry_data %>% mutate(Vessel = "V1")
    slurry_v3 <- slurry_data %>% mutate(Vessel = "V3")

    plot_data <- plot_data %>%
      filter(Vessel != "Slurry") %>%
      bind_rows(slurry_v1, slurry_v3)
  } else if (!include_slurry) {
    plot_data <- plot_data %>% filter(Vessel != "Slurry")
  }

  # Filter for vessel
  if (vessel != "both") {
    plot_data <- plot_data %>% filter(Vessel == vessel)
  }

  # Only keep vessels that have data
  plot_data <- plot_data %>%
    group_by(Vessel) %>%
    filter(n() > 0) %>%
    ungroup()

  return(plot_data)
}

#' Create time series plot showing INDIVIDUAL DONORS
#' Each donor gets their own color (pink shades for LC, blue shades for Rec)
#'
#' @param data Summary data frame (from combine_results()$summary)
#' @param target Bacterial target to plot
#' @param vessel Which vessel to plot ("V1", "V3", or "both")
#' @param log_scale Use log10 scale for y-axis (default TRUE)
#' @param include_slurry Include slurry samples as Day 0 (default TRUE)
#' @return ggplot object
plot_individual_donors <- function(data, target, vessel = "both",
                                    log_scale = TRUE, include_slurry = TRUE) {

  plot_data <- prepare_plot_data(data, target, vessel, include_slurry)

  if (nrow(plot_data) == 0) {
    message("No data available for ", target, " in vessel ", vessel)
    return(NULL)
  }

  # Get donors present in data and their colors
  donors_present <- unique(plot_data$Donor)
  colors_to_use <- DONOR_COLORS[donors_present]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  p <- ggplot(plot_data, aes(x = Day_Numeric, y = Mean_Copy_Number,
                              color = Donor, group = Donor)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = colors_to_use) +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 8, 10, 13, 17, 20)) +
    labs(
      title = paste(target, "- Individual Donors"),
      x = "Day",
      y = "Copy Number",
      color = "Donor"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  if (log_scale) {
    p <- p +
      consistent_y_scale(plot_data) +
      labs(y = "Copy Number (log scale)")
  }

  # Facet by vessel if showing both - only show panels with data
  # Using shared y-axis for comparability across vessels
  if (vessel == "both" && length(unique(plot_data$Vessel)) > 1) {
    p <- p + facet_wrap(~Vessel, ncol = 2)
  } else if (vessel == "both") {
    v <- unique(plot_data$Vessel)[1]
    p <- p + labs(title = paste(target, "-", v, "- Individual Donors"))
  }

  return(p)
}

#' Create time series plot showing GROUP MEANS with confidence ribbons
#' Shows average across donors within each group (LC, Recovered, Pool)
#' Error bars show BIOLOGICAL variation (SD across donors in group)
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot
#' @param vessel Which vessel to plot ("V1", "V3", or "both")
#' @param show_ribbon Show SD ribbon around mean (default TRUE)
#' @param show_n Show sample size (N) labels on plot (default TRUE)
#' @return ggplot object
plot_group_means <- function(data, target, vessel = "both", show_ribbon = TRUE, show_n = TRUE) {

  plot_data <- prepare_plot_data(data, target, vessel, include_slurry = TRUE)

  if (nrow(plot_data) == 0) {
    message("No data available for ", target, " in vessel ", vessel)
    return(NULL)
  }

  # Calculate group means and SD (BIOLOGICAL variation - across donors)
  group_summary <- plot_data %>%
    group_by(Donor_Group, Vessel, Day_Numeric) %>%
    summarise(
      Group_Mean = mean(Mean_Copy_Number, na.rm = TRUE),
      Group_SD = sd(Mean_Copy_Number, na.rm = TRUE),
      Group_SE = sd(Mean_Copy_Number, na.rm = TRUE) / sqrt(n()),
      N = n(),
      .groups = 'drop'
    ) %>%
    filter(!is.na(Group_Mean))

  # Replace NA SD with 0 (for single observations)
  group_summary <- group_summary %>%
    mutate(
      Group_SD = ifelse(is.na(Group_SD), 0, Group_SD),
      Group_SE = ifelse(is.na(Group_SE), 0, Group_SE)
    )

  # Check for varying N (indicates missing data)
  n_varies <- group_summary %>%
    group_by(Donor_Group, Vessel) %>%
    summarise(n_unique = n_distinct(N), .groups = 'drop') %>%
    pull(n_unique) %>%
    any(. > 1)

  p <- ggplot(group_summary, aes(x = Day_Numeric, y = Group_Mean,
                                  color = Donor_Group, fill = Donor_Group,
                                  group = Donor_Group))

  if (show_ribbon && any(group_summary$Group_SD > 0)) {
    p <- p +
      geom_ribbon(aes(ymin = pmax(Group_Mean - Group_SD, 1),
                      ymax = Group_Mean + Group_SD),
                  alpha = 0.2, color = NA)
  }

  p <- p +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4)

  # Add N labels if requested and N varies (indicates missing data)
  if (show_n && n_varies) {
    p <- p +
      geom_text(aes(label = paste0("n=", N)),
                vjust = -1.5, size = 2.5, show.legend = FALSE)
  }

  p <- p +
    scale_color_manual(values = GROUP_COLORS) +
    scale_fill_manual(values = GROUP_COLORS) +
    consistent_y_scale(plot_data) +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 6, 8, 10, 13, 17, 20)) +
    labs(
      title = paste(target, "- Group Averages"),
      subtitle = if(n_varies) "Ribbon = SD across donors; n varies due to missing data" else "Ribbon = SD across donors",
      x = "Day",
      y = "Copy Number (log scale)",
      color = "Group",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  # Facet by vessel - shared y-axis for comparability
  if (vessel == "both" && length(unique(group_summary$Vessel)) > 1) {
    p <- p + facet_wrap(~Vessel, ncol = 2)
  } else if (vessel == "both") {
    v <- unique(group_summary$Vessel)[1]
    p <- p + labs(title = paste(target, "-", v, "- Group Averages"))
  }

  return(p)
}

#' Create multi-panel plot for all bacterial targets
#'
#' @param data Summary data frame
#' @param vessel Which vessel to plot
#' @return ggplot object
plot_all_targets <- function(data, vessel = "V1") {

  plot_data <- prepare_plot_data(data, target = NULL, vessel = vessel, include_slurry = TRUE)

  if (nrow(plot_data) == 0) {
    message("No data available for vessel ", vessel)
    return(NULL)
  }

  # Filter to only the specified vessel
  if (vessel != "both") {
    plot_data <- plot_data %>% filter(Vessel == vessel)
  }

  # Get donors present and their colors
  donors_present <- unique(plot_data$Donor)
  colors_to_use <- DONOR_COLORS[donors_present]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  p <- ggplot(plot_data, aes(x = Day_Numeric, y = Mean_Copy_Number,
                              color = Donor, group = Donor)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2) +
    scale_color_manual(values = colors_to_use) +
    consistent_y_scale(plot_data) +
    scale_x_continuous(breaks = c(0, 3, 6, 10, 13, 17, 20)) +
    facet_wrap(~Bacterial_Target, ncol = 3) +
    labs(
      title = paste("All Bacterial Targets -", vessel),
      x = "Day",
      y = "Copy Number (log scale)",
      color = "Donor"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2))

  return(p)
}

#' Create bar chart comparing baseline (slurry) copy numbers across donors
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot (or NULL for all targets faceted)
#' @return ggplot object
plot_slurry_comparison <- function(data, target = NULL) {

  # Filter for slurry samples only
  plot_data <- data %>%
    filter(Vessel == "Slurry" | Day == "Slurry" | Day_Numeric == 0)

  if (!is.null(target)) {
    plot_data <- plot_data %>% filter(Bacterial_Target == target)
  }

  if (nrow(plot_data) == 0) {
    message("No slurry/baseline data available")
    return(NULL)
  }

  # Get donors and colors
  donors_present <- unique(plot_data$Donor)
  colors_to_use <- DONOR_COLORS[donors_present]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  # Order donors by group
  plot_data <- plot_data %>%
    mutate(Donor = factor(Donor, levels = c("LC01", "LC02", "LC03", "LC04",
                                             "Rec01", "Rec02", "Rec03", "POOL")))

  p <- ggplot(plot_data, aes(x = Donor, y = Mean_Copy_Number, fill = Donor)) +
    geom_col(width = 0.7, color = "black", linewidth = 1.0) +
    geom_errorbar(aes(ymin = pmax(Mean_Copy_Number - SD_Copy_Number, 1),
                      ymax = Mean_Copy_Number + SD_Copy_Number),
                  width = 0.2, linewidth = 0.7) +
    scale_fill_manual(values = colors_to_use) +
    consistent_y_scale(plot_data) +
    labs(
      title = if(is.null(target)) "Baseline (Slurry) Comparison - All Targets"
              else paste(target, "- Baseline (Slurry) Comparison"),
      x = "Donor",
      y = "Copy Number (log scale)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  # Facet by target if showing all - shared y-axis
  if (is.null(target) && length(unique(plot_data$Bacterial_Target)) > 1) {
    p <- p + facet_wrap(~Bacterial_Target, ncol = 3)
  }

  return(p)
}

#' Create start vs end comparison plot (dumbbell/slope chart)
#' Shows change from baseline (Day 0) to final timepoint
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot
#' @param vessel Which vessel ("V1" or "V3")
#' @param end_day The final day to compare (default finds max day in data)
#' @return ggplot object
plot_start_vs_end <- function(data, target, vessel = "V1", end_day = NULL) {

  # Prepare data
  plot_data <- prepare_plot_data(data, target, vessel, include_slurry = TRUE)

  if (nrow(plot_data) == 0) {
    message("No data available for ", target, " in vessel ", vessel)
    return(NULL)
  }

  # Find end day if not specified
  if (is.null(end_day)) {
    end_day <- max(plot_data$Day_Numeric, na.rm = TRUE)
  }

  # Get start (Day 0) and end data
  start_data <- plot_data %>%
    filter(Day_Numeric == 0) %>%
    select(Donor, Donor_Group, Start_Copy = Mean_Copy_Number)

  end_data <- plot_data %>%
    filter(Day_Numeric == end_day) %>%
    select(Donor, End_Copy = Mean_Copy_Number)

  # Combine
  comparison_data <- start_data %>%
    inner_join(end_data, by = "Donor") %>%
    mutate(
      Fold_Change = End_Copy / Start_Copy,
      Log2_FC = log2(Fold_Change),
      Direction = ifelse(End_Copy > Start_Copy, "Increased", "Decreased")
    )

  if (nrow(comparison_data) == 0) {
    message("No paired start/end data available")
    return(NULL)
  }

  # Get colors
  donors_present <- unique(comparison_data$Donor)
  colors_to_use <- DONOR_COLORS[donors_present]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  # Order donors
  comparison_data <- comparison_data %>%
    mutate(Donor = factor(Donor, levels = c("LC01", "LC02", "LC03", "LC04",
                                             "Rec01", "Rec02", "Rec03", "POOL")))

  # Create dumbbell plot
  p <- ggplot(comparison_data, aes(y = Donor)) +
    # Line connecting start to end
    geom_segment(aes(x = Start_Copy, xend = End_Copy, yend = Donor, color = Donor),
                 linewidth = 1.5, alpha = 0.7) +
    # Start point (circle)
    geom_point(aes(x = Start_Copy, color = Donor), size = 4, shape = 16) +
    # End point (triangle)
    geom_point(aes(x = End_Copy, color = Donor), size = 4, shape = 17) +
    scale_color_manual(values = colors_to_use) +
    scale_x_log10(
      labels = scales::scientific,
      breaks = scales::trans_breaks("log10", function(x) 10^x)
    ) +
    labs(
      title = paste(target, "-", vessel, ": Day 0 vs Day", end_day),
      subtitle = "Circle = Start (Day 0), Triangle = End",
      x = "Copy Number (log scale)",
      y = "Donor"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Create fold change bar chart (start to end)
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot
#' @param vessel Which vessel ("V1" or "V3")
#' @param end_day The final day to compare
#' @return ggplot object
plot_fold_change <- function(data, target, vessel = "V1", end_day = NULL) {

  # Prepare data
  plot_data <- prepare_plot_data(data, target, vessel, include_slurry = TRUE)

  if (nrow(plot_data) == 0) {
    message("No data available for ", target, " in vessel ", vessel)
    return(NULL)
  }

  # Find end day if not specified
  if (is.null(end_day)) {
    end_day <- max(plot_data$Day_Numeric, na.rm = TRUE)
  }

  # Get start and end data
  start_data <- plot_data %>%
    filter(Day_Numeric == 0) %>%
    select(Donor, Donor_Group, Start_Copy = Mean_Copy_Number)

  end_data <- plot_data %>%
    filter(Day_Numeric == end_day) %>%
    select(Donor, End_Copy = Mean_Copy_Number)

  # Calculate fold change
  fc_data <- start_data %>%
    inner_join(end_data, by = "Donor") %>%
    mutate(
      Log2_FC = log2(End_Copy / Start_Copy),
      Direction = ifelse(Log2_FC > 0, "Increased", "Decreased"),
      Donor = factor(Donor, levels = c("LC01", "LC02", "LC03", "LC04",
                                        "Rec01", "Rec02", "Rec03", "POOL"))
    )

  if (nrow(fc_data) == 0) {
    message("No paired start/end data available")
    return(NULL)
  }

  # Get colors
  donors_present <- unique(fc_data$Donor)
  colors_to_use <- DONOR_COLORS[as.character(donors_present)]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  p <- ggplot(fc_data, aes(x = Donor, y = Log2_FC, fill = Donor)) +
    geom_col(width = 0.7, color = "black", linewidth = 1.0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    scale_fill_manual(values = colors_to_use) +
    labs(
      title = paste(target, "-", vessel, ": Fold Change (Day 0 to Day", end_day, ")"),
      x = "Donor",
      y = expression(Log[2]~Fold~Change),
      caption = "Above 0 = increased, Below 0 = decreased"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  return(p)
}

#' Create comparison plot showing both vessels side by side for start vs end
#'
#' @param data Summary data frame
#' @param target Bacterial target to plot
#' @param end_day The final day to compare
#' @return ggplot object
plot_start_vs_end_both_vessels <- function(data, target, end_day = NULL) {

  # Get data for both vessels
  plot_data <- data %>%
    filter(Bacterial_Target == target) %>%
    filter(!is.na(Mean_Copy_Number) & !is.na(Day_Numeric) & !is.na(Donor))

  # Handle slurry - duplicate for both vessels
  if ("Slurry" %in% plot_data$Vessel) {
    slurry_data <- plot_data %>% filter(Vessel == "Slurry")
    slurry_v1 <- slurry_data %>% mutate(Vessel = "V1")
    slurry_v3 <- slurry_data %>% mutate(Vessel = "V3")

    plot_data <- plot_data %>%
      filter(Vessel != "Slurry") %>%
      bind_rows(slurry_v1, slurry_v3)
  }

  if (nrow(plot_data) == 0) {
    message("No data available for ", target)
    return(NULL)
  }

  # Find end day
  if (is.null(end_day)) {
    end_day <- max(plot_data$Day_Numeric, na.rm = TRUE)
  }

  # Calculate fold change for each vessel
  fc_data <- plot_data %>%
    filter(Day_Numeric %in% c(0, end_day)) %>%
    select(Donor, Donor_Group, Vessel, Day_Numeric, Mean_Copy_Number) %>%
    pivot_wider(names_from = Day_Numeric, values_from = Mean_Copy_Number,
                names_prefix = "Day_") %>%
    mutate(
      Log2_FC = log2(get(paste0("Day_", end_day)) / Day_0),
      Donor = factor(Donor, levels = c("LC01", "LC02", "LC03", "LC04",
                                        "Rec01", "Rec02", "Rec03", "POOL"))
    ) %>%
    filter(!is.na(Log2_FC) & !is.infinite(Log2_FC))

  if (nrow(fc_data) == 0) {
    message("No paired start/end data available")
    return(NULL)
  }

  # Get colors
  donors_present <- unique(fc_data$Donor)
  colors_to_use <- DONOR_COLORS[as.character(donors_present)]
  colors_to_use <- colors_to_use[!is.na(colors_to_use)]

  p <- ggplot(fc_data, aes(x = Donor, y = Log2_FC, fill = Donor)) +
    geom_col(width = 0.7, color = "black", linewidth = 1.0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    scale_fill_manual(values = colors_to_use) +
    facet_wrap(~Vessel, ncol = 2) +
    labs(
      title = paste(target, ": Fold Change (Day 0 to Day", end_day, ")"),
      x = "Donor",
      y = expression(Log[2]~Fold~Change),
      caption = "Above 0 = increased, Below 0 = decreased"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  return(p)
}

# Keep old function names as aliases for compatibility
plot_timeseries <- plot_individual_donors
plot_timeseries_by_group <- plot_group_means

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
  plots_saved <- 0

  # Save individual target plots
  for (target in targets) {
    cat("  Creating plots for:", target, "\n")

    # Individual donors - V1
    p <- plot_individual_donors(combined_results$summary, target, vessel = "V1")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_V1_individuals_", timestamp, ".pdf")),
             p, width = 10, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Individual donors - V3
    p <- plot_individual_donors(combined_results$summary, target, vessel = "V3")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_V3_individuals_", timestamp, ".pdf")),
             p, width = 10, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Individual donors - both vessels
    p <- plot_individual_donors(combined_results$summary, target, vessel = "both")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_both_individuals_", timestamp, ".pdf")),
             p, width = 12, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Group means - both vessels
    p <- plot_group_means(combined_results$summary, target, vessel = "both")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_group_means_", timestamp, ".pdf")),
             p, width = 12, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Slurry baseline comparison
    p <- plot_slurry_comparison(combined_results$summary, target)
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_slurry_baseline_", timestamp, ".pdf")),
             p, width = 10, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Start vs end fold change - both vessels
    p <- plot_start_vs_end_both_vessels(combined_results$summary, target)
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_fold_change_", timestamp, ".pdf")),
             p, width = 12, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Start vs end dumbbell - V1
    p <- plot_start_vs_end(combined_results$summary, target, vessel = "V1")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_V1_start_vs_end_", timestamp, ".pdf")),
             p, width = 10, height = 6)
      plots_saved <- plots_saved + 1
    }

    # Start vs end dumbbell - V3
    p <- plot_start_vs_end(combined_results$summary, target, vessel = "V3")
    if (!is.null(p)) {
      ggsave(file.path(output_dir, paste0(target, "_V3_start_vs_end_", timestamp, ".pdf")),
             p, width = 10, height = 6)
      plots_saved <- plots_saved + 1
    }
  }

  # Save multi-panel overview
  cat("  Creating overview plots...\n")

  p <- plot_all_targets(combined_results$summary, "V1")
  if (!is.null(p)) {
    ggsave(file.path(output_dir, paste0("ALL_targets_V1_", timestamp, ".pdf")),
           p, width = 14, height = 10)
    plots_saved <- plots_saved + 1
  }

  p <- plot_all_targets(combined_results$summary, "V3")
  if (!is.null(p)) {
    ggsave(file.path(output_dir, paste0("ALL_targets_V3_", timestamp, ".pdf")),
           p, width = 14, height = 10)
    plots_saved <- plots_saved + 1
  }

  # Save slurry comparison for all targets
  p <- plot_slurry_comparison(combined_results$summary, target = NULL)
  if (!is.null(p)) {
    ggsave(file.path(output_dir, paste0("ALL_slurry_baseline_", timestamp, ".pdf")),
           p, width = 14, height = 10)
    plots_saved <- plots_saved + 1
  }

  cat("Saved", plots_saved, "plots to:", output_dir, "\n")
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
