# cleanRfish V2: Uncertainty-Weighted Trajectory Reconstruction
# ============================================================
# Bidirectional Spline-Based Segment Linking with 
# Ornstein-Uhlenbeck Uncertainty Modelling

# Suppress R CMD check notes for dplyr/ggplot2 NSE
utils::globalVariables(c(
  "time", "x", "y", "segment", "jump_flag", "delta_x", "delta_y", "delta_t",
  "vx", "vy", "speed", "distance", "direction", "delta_theta", "delta_speed",
  "z_speed", "z_dir", "z_combo", "start_time", "end_time", "tracking_group",
  "gap_duration", "segment_before", "segment_after", "n_points", "total_movement",
  "x_original", "y_original", "x_reconstructed", "y_reconstructed", "x_smooth", "y_smooth",
  "individual_number", "start_x", "start_y", "end_x", "end_y", "pred_x", "pred_y",
  "residual_gap", "spatial_gap", "theoretical_velocity", "sigma_t", "log_liklihood",
  "feasible", "rank", "positional_deviation", "speed_deviation", "directional_deviation",
  "mean_positional", "sd_positional", "mean_speed", "mean_direction", "n", "is_projection",
  "fit_xy", "fit_x", "fit_y", "uncertainty", "x_filled", "y_filled", "recon_x", "recon_y",
  "next_start", "gap_start", "gap_end", "segment_duration", "segment_length",
  ".", "chosen", "x.x", "x.y", "y.y", "xy_pred_candidate", "thr", ".data",
  "acceleration", "residuals_x", "residuals_y", "theoretical_acceleration"
))

# 1. TRAJECTORY SEGMENTATION ====================================

#' Detect Jumps in Trajectory Data
#'
#' Identifies abrupt changes in trajectory that indicate tracking errors or identity swaps
#' using combined standardised velocity and directional change metrics.
#'
#' @param df Data frame with columns: time, x, y
#' @param z_thresh Z-score threshold for jump detection (default: 3)
#' @return Data frame with additional columns for jump detection metrics
#' @export
detect_jumps <- function(df, z_thresh = 3) {
  df <- df |> 
    dplyr::arrange(time) |> 
    dplyr::mutate(jump_flag = ifelse(is.na(dplyr::lag(x)) & is.na(dplyr::lag(y)), TRUE, FALSE)) |>
    dplyr::filter(!is.na(x), !is.na(y)) |>
    dplyr::mutate(
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),
      vx = delta_x / delta_t,
      vy = delta_y / delta_t,
      speed = sqrt(vx^2 + vy^2),
      distance = sqrt(delta_x^2 + delta_y^2),
      direction = atan2(vy, vx)
    )
  
  df <- df |> mutate(
    delta_theta = abs(direction - dplyr::lag(direction)),
    delta_theta = ifelse(delta_theta > pi, 2*pi - delta_theta, delta_theta),
    delta_speed = abs(speed - dplyr::lag(speed))
  )
  
  df <- df |> mutate(
    z_speed = scale(speed)[,1],
    z_dir   = scale(delta_theta)[,1],
    z_combo = sqrt(z_speed^2 + z_dir^2)
  )
  
  df <- df |> mutate(jump_flag = ifelse(jump_flag == FALSE & z_combo > z_thresh, TRUE, jump_flag))
  df
}

#' Assign Segment IDs Based on Jump Flags
#'
#' Converts jump flags into contiguous segment identifiers.
#'
#' @param df Data frame with jump_flag column
#' @return Data frame with segment column added
#' @export
assign_segments <- function(df){
  df <- df |> arrange(time)
  js <- df$jump_flag
  js[is.na(js)] <- FALSE
  if (length(js) > 0) js[1] <- FALSE
  df$segment <- 1L + cumsum(js)
  df
}

#' Get Speed Threshold for Feasibility Filtering
#'
#' Calculate biologically plausible maximum velocity using log-normal distribution.
#'
#' @param df Data frame with speed column
#' @param prob Quantile probability (default: 0.999)
#' @param eps Small constant to prevent log(0) (default: 1e-10)
#' @return Numeric speed threshold
#' @export
get_speed_threshold <- function(df, prob = 0.999, eps = 1e-10) {
  z <- log(pmax(df$speed, 0) + eps)
  exp(stats::quantile(z, probs = prob, na.rm = TRUE))
}

# 2. TRACKING GAP MANAGEMENT ====================================

#' Identify Large Tracking Gaps
#'
#' Detect large temporal gaps between segments that indicate complete tracking loss.
#'
#' @param df_all Data frame with segment and time columns
#' @param window Time window threshold for identifying gaps
#' @return Data frame of identified gaps
#' @export
identify_tracking_gaps <- function(df_all, window) {
  # Get time range for each segment
  segment_ranges <- df_all |>
    dplyr::group_by(segment) |>
    dplyr::summarise(
      start_time = min(time),
      end_time = max(time),
      .groups = "drop"
    ) |>
    dplyr::arrange(start_time)
  
  # Calculate gaps between consecutive segments
  gaps <- segment_ranges |>
    dplyr::mutate(
      next_start = dplyr::lead(start_time),
      segment_after = dplyr::lead(segment),
      gap_duration = next_start - end_time
    ) |>
    dplyr::filter(!is.na(gap_duration), gap_duration > window, !is.na(segment_after)) |>
    dplyr::select(
      gap_start = end_time,
      gap_end = next_start,
      gap_duration,
      segment_before = segment,
      segment_after
    )
  
  return(gaps)
}

#' Find Ground Truth Segment
#'
#' Identify the most reliable segment to serve as an anchor point for reconstruction.
#'
#' @param df Data frame with segment, distance, and time columns
#' @param min_movement Minimum total movement required (default: 500)
#' @param group_id Optional tracking group ID to filter by
#' @return Integer segment ID of ground truth segment
#' @export
find_ground_truth_segment <- function(df, min_movement = 500, group_id = NULL) {
  # Filter by group if specified
  if (!is.null(group_id)) {
    df <- df |> filter(tracking_group == group_id)
    if (nrow(df) == 0) {
      stop(sprintf("No data found for tracking_group %d", group_id))
    }
  }
  
  segment_summary <- df |> 
    dplyr::group_by(segment) |> 
    dplyr::summarise(
      n_points = dplyr::n(),
      total_movement = sum(distance, na.rm = TRUE),
      .groups = "drop"
    ) |> 
    dplyr::filter(total_movement > min_movement) |> 
    dplyr::arrange(dplyr::desc(n_points))
  
  if (nrow(segment_summary) == 0) {
    # If no segment meets movement threshold, just take the longest one
    segment_summary <- df |>
      dplyr::group_by(segment) |>
      dplyr::summarise(n_points = dplyr::n(), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(n_points))
  }
  
  # Select the segment with the longest duration and sufficient movement
  ground_truth_segment_id <- segment_summary$segment[1]
  return(ground_truth_segment_id)
}

#' Identify Isolated Tracking Groups
#'
#' Partition trajectory segments into isolated tracking groups based on large temporal gaps.
#'
#' @param df_all Data frame with segment and time columns
#' @param window Time window threshold
#' @return Data frame with tracking_group column added
#' @export
identify_isolated_groups <- function(df_all, window) {
  # Get gaps
  gaps <- identify_tracking_gaps(df_all, window)
  
  if (nrow(gaps) == 0) {
    # No large gaps found - everything is one group
    df_all <- df_all |>
      dplyr::mutate(tracking_group = 1L)
    return(df_all)
  }
  
  # Get segment ranges
  segment_ranges <- df_all |>
    dplyr::group_by(segment) |>
    dplyr::summarise(
      start_time = min(time),
      end_time = max(time),
      .groups = "drop"
    ) |>
    dplyr::arrange(start_time)
  
  # Assign groups based on gaps
  segment_ranges <- segment_ranges |>
    dplyr::mutate(tracking_group = 1L)
  
  current_group <- 1L
  for (i in seq_len(nrow(gaps))) {
    current_group <- current_group + 1L
    segment_after <- gaps$segment_after[i]
    segment_ranges <- segment_ranges |>
      dplyr::mutate(tracking_group = ifelse(segment >= segment_after, current_group, tracking_group))
  }
  
  # Join back to original data
  df_all <- df_all |>
    dplyr::left_join(segment_ranges |> select(segment, tracking_group), by = "segment")
  
  return(df_all)
}

# 3. UNCERTAINTY QUANTIFICATION ==================================

#' Compute Uncertainty Model (Ornstein-Uhlenbeck Process)
#'
#' Quantify how positional uncertainty grows as a function of time gap.
#'
#' @param df_all Data frame with segment, time, x, y, speed, direction columns
#' @param window Maximum time window to consider
#' @param n_boot Number of bootstrap samples (default: 10000)
#' @return Function that predicts uncertainty for a given time gap
#' @export
compute_uncertainty_model <- function(df_all, window, n_boot = 10000) {
  
  # Step 1: Compute empirical uncertainty through bootstrap sampling
  # Filter segments with sufficient points upfront
  segment_counts <- df_all |>
    dplyr::group_by(segment) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n >= 20)
  
  valid_segments <- segment_counts$segment
  
  if(length(valid_segments) == 0) {
    warning("No valid segments for uncertainty computation, returning constant uncertainty")
    return(function(time_gap) rep(100, length(time_gap)))
  }
  
  # Pre-filter and prepare data once
  seg_data <- df_all |>
    dplyr::filter(segment %in% valid_segments) |>
    dplyr::arrange(segment, time) |>
    dplyr::select(segment, time, x, y, speed, direction)
  
  # Pre-allocate result list
  results <- vector("list", length(valid_segments))
  
  # Process each segment in vectorized batches
  for(i in seq_along(valid_segments)) {
    seg <- valid_segments[i]
    seg_subset <- seg_data |> filter(segment == seg)
    n_points <- nrow(seg_subset)
    
    # Generate all random pairs at once (vectorized)
    idx1 <- sample(1:n_points, n_boot, replace = TRUE)
    idx2 <- sample(1:n_points, n_boot, replace = TRUE)
    
    # Ensure time ordering and different indices
    valid <- idx1 != idx2
    idx_min <- pmin(idx1[valid], idx2[valid])
    idx_max <- pmax(idx1[valid], idx2[valid])
    
    # Extract values vectorized
    point1 <- seg_subset[idx_min, ]
    point2 <- seg_subset[idx_max, ]
    
    # Compute all metrics vectorized
    time_gap <- as.numeric(point2$time - point1$time)
    
    # Filter valid time gaps (use window parameter)
    valid_gaps <- time_gap > 0 & time_gap <= window
    
    if(sum(valid_gaps) > 0) {
      # Vectorized deviation calculations
      dx <- point2$x[valid_gaps] - point1$x[valid_gaps]
      dy <- point2$y[valid_gaps] - point1$y[valid_gaps]
      
      results[[i]] <- tibble(
        segment = seg,
        time_gap = time_gap[valid_gaps],
        positional_deviation = sqrt(dx^2 + dy^2),
        speed_deviation = abs(point2$speed[valid_gaps] - point1$speed[valid_gaps]),
        directional_deviation = abs(atan2(sin(point2$direction[valid_gaps] - point1$direction[valid_gaps]), 
                                          cos(point2$direction[valid_gaps] - point1$direction[valid_gaps]))),
        initial_speed = point1$speed[valid_gaps],
        initial_direction = point1$direction[valid_gaps]
      )
    }
  }
  
  empirical_uncertainty <- bind_rows(results)
  
  if(nrow(empirical_uncertainty) == 0) {
    warning("No empirical uncertainty data, returning constant uncertainty")
    return(function(time_gap) rep(100, length(time_gap)))
  }
  
  # Step 2: Summarise uncertainty by time gap
  uncertainty_summary <- empirical_uncertainty |>
    dplyr::group_by(time_gap) |>
    dplyr::summarise(
      mean_positional = mean(positional_deviation, na.rm = TRUE),
      sd_positional = sd(positional_deviation, na.rm = TRUE),
      mean_speed = mean(speed_deviation, na.rm = TRUE),
      mean_direction = mean(directional_deviation, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::filter(n > 10)
  
  # Extract time gaps and observed standard deviations
  time_gap_vec <- uncertainty_summary$time_gap
  observed_sd <- uncertainty_summary$sd_positional
  
  # Remove any NA or non-finite values
  valid_idx <- !is.na(observed_sd) & is.finite(observed_sd)
  time_gap_vec <- time_gap_vec[valid_idx]
  observed_sd <- observed_sd[valid_idx]
  
  if(length(time_gap_vec) < 3) {
    warning("Insufficient data for model fitting, returning constant uncertainty")
    return(function(time_gap) rep(mean(observed_sd, na.rm = TRUE), length(time_gap)))
  }
  
  # Step 3: Fit Ornstein-Uhlenbeck saturation model
  ou_model <- tryCatch({
    nls(observed_sd ~ a * sqrt((1 - exp(-2 * beta * time_gap_vec)) / (2 * beta)),
        start = list(a = max(observed_sd), beta = 0.3))
  }, error = function(e) {
    # Fallback to simplified saturation form if full OU doesn't converge
    tryCatch({
      nls(observed_sd ~ a * (1 - exp(-b * time_gap_vec)),
          start = list(a = max(observed_sd), b = 0.5))
    }, error = function(e2) {
      warning("Model fitting failed, returning linear approximation")
      return(NULL)
    })
  })
  
  # Step 4: Return prediction function
  if(is.null(ou_model)) {
    # Fallback: linear interpolation/extrapolation
    return(function(time_gap) {
      approx(time_gap_vec, observed_sd, xout = time_gap, rule = 2)$y
    })
  } else {
    # Return OU model prediction function
    return(function(time_gap) {
      # Create a data frame with the correct variable name used in model fitting
      pred_data <- data.frame(time_gap_vec = time_gap)
      uncertainty <- predict(ou_model, newdata = pred_data)
      pmax(uncertainty, 0)  # Ensure non-negative
    })
  }
}

#' Calculate Log-Likelihood for Candidate Ranking
#'
#' Score candidate segments based on spatial fit quality penalized by uncertainty.
#'
#' @param residual_gap Spatial distance between prediction and observation
#' @param sigma_t Time-dependent uncertainty
#' @return Numeric log-likelihood score
#' @export
calculate_log_liklihood <- function(residual_gap, sigma_t) {
  -sigma_t - 0.5 * residual_gap
}

# 4. VISUALISATION UTILITIES ====================================

#' Generate Diagnostic Time-Series Plot
#'
#' Create time-series diagnostic plots showing spline fits, uncertainty envelopes,
#' and candidate segments.
#'
#' @param data_for_spline Data frame with observed data used for spline fitting
#' @param plot_df Data frame with spline predictions and uncertainty
#' @param features_df Data frame with candidate segment features
#' @param current_seg_number Current segment number
#' @param ground_truth_segment Ground truth segment ID
#' @param direction Either "forward" or "backward"
#' @return ggplot object
#' @keywords internal
generate_diagnostic_plot <- function(data_for_spline, plot_df, features_df,
                                      current_seg_number, ground_truth_segment, direction = "forward") {
  
  coord_col <- if (direction == "forward") "start" else "end"
  
  y_min <- min(c(data_for_spline$x + data_for_spline$y,
                 features_df[[paste0(coord_col, "_x")]] + features_df[[paste0(coord_col, "_y")]]),
               na.rm = TRUE) - 100
  y_max <- max(c(data_for_spline$x + data_for_spline$y,
                 plot_df$fit_xy,
                 features_df[[paste0(coord_col, "_x")]] + features_df[[paste0(coord_col, "_y")]]),
               na.rm = TRUE) + 100
  
  time_ref <- if (direction == "forward") max(data_for_spline$time) else min(data_for_spline$time)
  plot_df_projection <- plot_df |>
    dplyr::mutate(is_projection = if (direction == "forward") time > time_ref else time < time_ref)
  
  ggplot2::ggplot() +
    ggplot2::geom_point(data = data_for_spline, ggplot2::aes(x = time, y = x + y),
               colour = "darkblue", alpha = 0.6) +
    ggplot2::geom_line(data = plot_df_projection, ggplot2::aes(x = time, y = fit_xy),
              colour = "purple", linewidth = 0.5) +
    ggplot2::geom_ribbon(data = dplyr::filter(plot_df_projection, is_projection),
                ggplot2::aes(x = time, ymin = fit_xy - uncertainty, ymax = fit_xy + uncertainty),
                fill = "purple", alpha = 0.2) +
    ggplot2::geom_point(data = features_df,
               ggplot2::aes(x = .data[[paste0(coord_col, "_time")]], 
                   y = pred_x + pred_y, colour = factor(segment)),
               size = 2, shape = 16, alpha = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(features_df, rank == 1),
               ggplot2::aes(x = .data[[paste0(coord_col, "_time")]], 
                   y = .data[[paste0(coord_col, "_x")]] + .data[[paste0(coord_col, "_y")]]),
               size = 6, shape = 1, colour = "black") +
    ggplot2::geom_point(data = dplyr::filter(features_df, !feasible),
               ggplot2::aes(x = .data[[paste0(coord_col, "_time")]], 
                   y = .data[[paste0(coord_col, "_x")]] + .data[[paste0(coord_col, "_y")]]),
               size = 6, shape = 4, colour = "black") +
    ggplot2::geom_point(data = features_df,
               ggplot2::aes(x = .data[[paste0(coord_col, "_time")]], 
                   y = .data[[paste0(coord_col, "_x")]] + .data[[paste0(coord_col, "_y")]], 
                   colour = factor(segment),
                   text = paste0("Segment: ", segment,
                                 "\nRank: ", ifelse(is.na(rank), "Infeasible", rank),
                                 "\nlog_liklihood: ", round(log_liklihood, 6),
                                 "\nΔt: ", round(delta_t, 2), "s",
                                 "\nResidual: ", round(residual_gap, 2),
                                 "\nσ(t): ", round(sigma_t, 2),
                                 "\nSpatial Gap: ", round(spatial_gap, 2),
                                 "\nTheoretical Vel: ", round(theoretical_velocity, 2),
                                 "\nTheoretical Accel: ", round(theoretical_acceleration, 2),
                                 "\nPos (pred): (", round(pred_x, 2), ", ", round(pred_y, 2), ")",
                                 "\nPos (", coord_col, "): (", round(.data[[paste0(coord_col, "_x")]], 2), 
                                 ", ", round(.data[[paste0(coord_col, "_y")]], 2), ")",
                                 "\nFeasible: ", ifelse(feasible, "Yes", "No"))),
               size = 3, shape = 18) +
    ggplot2::geom_text(data = features_df,
              ggplot2::aes(x = .data[[paste0(coord_col, "_time")]], 
                  y = .data[[paste0(coord_col, "_x")]] + .data[[paste0(coord_col, "_y")]], 
                  label = segment, colour = factor(segment)),
              vjust = -1.2, show.legend = FALSE, size = 3.2) +
    ggplot2::labs(x = "Time (s)", y = "composite spatial position (X + Y)",
         title = paste0("Uncertainty-weighted ", toupper(direction), " - Segment: ", current_seg_number)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max))
}

#' Generate Spatial 2D Plot
#'
#' Create 2D spatial plots showing trajectories in physical coordinates.
#'
#' @param data_for_spline Data frame with observed positions used for spline fitting
#' @param plot_df_spatial Data frame with spatial spline predictions
#' @param features_df Data frame with candidate segment features
#' @param current_seg_number Current segment number
#' @param ground_truth_segment Ground truth segment ID
#' @param direction Either "forward" or "backward"
#' @return ggplot object
#' @keywords internal
generate_spatial_plot <- function(data_for_spline, plot_df_spatial, features_df,
                                   current_seg_number, ground_truth_segment, direction = "forward") {
  
  coord_col <- if (direction == "forward") "start" else "end"
  
  time_ref <- if (direction == "forward") max(data_for_spline$time) else min(data_for_spline$time)
  plot_df_spatial_projection <- plot_df_spatial |>
    dplyr::mutate(is_projection = if (direction == "forward") time > time_ref else time < time_ref)
  
  ggplot2::ggplot() +
    ggplot2::geom_point(data = data_for_spline, ggplot2::aes(x = x, y = y), colour = "darkblue", alpha = 0.6) +
    ggplot2::geom_path(data = plot_df_spatial_projection, ggplot2::aes(x = fit_x, y = fit_y), 
              colour = "purple", linewidth = 0.4) +
    ggplot2::geom_ribbon(data = dplyr::filter(plot_df_spatial_projection, is_projection),
                ggplot2::aes(x = fit_x, ymin = fit_y - uncertainty, ymax = fit_y + uncertainty),
                fill = "purple", alpha = 0.15) +
    ggplot2::geom_point(data = features_df,
               ggplot2::aes(x = pred_x, y = pred_y, colour = factor(segment)),
               size = 2, shape = 16, alpha = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(features_df, rank == 1), 
               ggplot2::aes(x = .data[[paste0(coord_col, "_x")]], y = .data[[paste0(coord_col, "_y")]]), 
               size = 6, shape = 1, colour = "black") +
    ggplot2::geom_point(data = dplyr::filter(features_df, !feasible), 
               ggplot2::aes(x = .data[[paste0(coord_col, "_x")]], y = .data[[paste0(coord_col, "_y")]]), 
               size = 6, shape = 4, colour = "black") +
    ggplot2::geom_point(data = features_df, 
               ggplot2::aes(x = .data[[paste0(coord_col, "_x")]], y = .data[[paste0(coord_col, "_y")]], 
                   colour = factor(segment),
                   text = paste0("Segment: ", segment,
                                 "\nRank: ", ifelse(is.na(rank), "Infeasible", rank),
                                 "\nlog_liklihood: ", round(log_liklihood, 6),
                                 "\nΔt: ", round(delta_t, 2), "s",
                                 "\nResidual: ", round(residual_gap, 2),
                                 "\nσ(t): ", round(sigma_t, 2),
                                 "\nSpatial Gap: ", round(spatial_gap, 2),
                                 "\nTheoretical Vel: ", round(theoretical_velocity, 2),
                                 "\nTheoretical Accel: ", round(theoretical_acceleration, 2),
                                 "\nPos (pred): (", round(pred_x, 2), ", ", round(pred_y, 2), ")",
                                 "\nPos (", coord_col, "): (", round(.data[[paste0(coord_col, "_x")]], 2), 
                                 ", ", round(.data[[paste0(coord_col, "_y")]], 2), ")",
                                 "\nFeasible: ", ifelse(feasible, "Yes", "No"))), 
               size = 3, shape = 18) +
    ggplot2::geom_text(data = features_df, 
              ggplot2::aes(x = .data[[paste0(coord_col, "_x")]], y = .data[[paste0(coord_col, "_y")]], 
                  label = segment, colour = factor(segment)), 
              vjust = -1.2, show.legend = FALSE, size = 3.2) +
    ggplot2::labs(x = "x", y = "y", 
         title = paste0("Uncertainty-weighted ", toupper(direction), " - Segment: ", current_seg_number, 
                        "\n", abs(current_seg_number - ground_truth_segment), " segments from ground truth")) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_fixed(ratio = 1,
                ylim = c(max(c(data_for_spline$y, features_df[[paste0(coord_col, "_y")]]), na.rm = TRUE),
                         min(c(data_for_spline$y, features_df[[paste0(coord_col, "_y")]]), na.rm = TRUE)),
                xlim = c(min(c(data_for_spline$x, features_df[[paste0(coord_col, "_x")]]), na.rm = TRUE), 
                         max(c(data_for_spline$x, features_df[[paste0(coord_col, "_x")]]), na.rm = TRUE)))
}

# 5. CANDIDATE SELECTION AND RANKING ============================

#' Get Candidate Segments within Time Window
#'
#' Find all segments with at least one point within the search window.
#'
#' @param df Data frame with segment and time columns
#' @param time_interval Two-element vector defining time range
#' @param current_segment Optional current segment ID to exclude
#' @param min_candidates Minimum candidates required (triggers window expansion)
#' @param direction Either "forward" or "backward"
#' @return Vector of candidate segment IDs
#' @export
get_candidate_segments <- function(df, time_interval, current_segment = NULL,
                                   min_candidates = 5, direction = "forward") {
  original_interval <- time_interval
  data_min_time <- min(df$time, na.rm = TRUE)
  data_max_time <- max(df$time, na.rm = TRUE)
  
  window_size <- abs(diff(time_interval))
  max_expansion_factor <- 5
  expansion_factor <- 1
  
  repeat {
    candidate_segments <- df |>
      dplyr::filter(time >= time_interval[1], time <= time_interval[2]) |>
      dplyr::pull(segment) |>
      unique()
    
    if (!is.null(current_segment)) {
      candidate_segments <- setdiff(candidate_segments, current_segment)
    }
    
    if (length(candidate_segments) >= min_candidates) {
      break
    }
    
    if (direction == "forward") {
      at_boundary <- time_interval[2] >= data_max_time
    } else {
      at_boundary <- time_interval[1] <= data_min_time
    }
    
    if (at_boundary || expansion_factor >= max_expansion_factor) {
      break
    }
    
    expansion_factor <- expansion_factor + 0.5
    expanded_window <- window_size * expansion_factor
    
    if (direction == "forward") {
      time_interval[2] <- original_interval[1] + expanded_window
    } else {
      time_interval[1] <- original_interval[2] - expanded_window
    }
  }
  
  return(candidate_segments)
}

#' Rank Candidate Segments Using Uncertainty-Weighted Likelihood
#'
#' Core algorithm for scoring and ranking candidate segments for linking.
#'
#' @param current_track Data frame of current trajectory
#' @param df_all Full data frame with all segments
#' @param window Search window size in time units
#' @param speed_threshold_quantile Quantile for speed feasibility (default: 0.999)
#' @param diagnostic_plots Generate diagnostic plots (default: FALSE)
#' @param compute_uncertainty Uncertainty model function (optional)
#' @param direction Either "forward" or "backward"
#' @param min_candidates Minimum candidates to consider (default: 5)
#' @return List with features_df, diagnostic_plot, and spatial_plot
#' @export
rank_candidates_uncertainty <- function(current_track,
                                        df_all,
                                        window = 30,
                                        speed_threshold_quantile = 0.9999,
                                        diagnostic_plots = FALSE,
                                        compute_uncertainty = NULL,
                                        direction = "forward",
                                        min_candidates = 5) {
  diagnostic_plot <- NULL
  spatial_plot <- NULL
  
  if (is.null(compute_uncertainty)) {
    compute_uncertainty <- compute_uncertainty_model(df_all, window = window)
  }
  
  if (direction == "forward") {
    data_for_spline <- current_track |> dplyr::arrange(time) |> dplyr::slice_tail(n = 75)
    current_seg_number <- (current_track |> dplyr::arrange(time) |> dplyr::slice_tail(n = 1))$segment
  } else {
    data_for_spline <- current_track |> dplyr::arrange(time) |> dplyr::slice_head(n = 75)
    current_seg_number <- (current_track |> dplyr::arrange(time) |> dplyr::slice_head(n = 1))$segment
  }
  
  if (nrow(data_for_spline) < 15) {
    message("Warning: Not enough data points (", nrow(data_for_spline), ") for spline fitting")
    return(list(features_df = tibble::tibble(), diagnostic_plot = NULL, spatial_plot = NULL))
  }
  
  mv_gam <- tryCatch({
    mgcv::gam(list(x ~ s(time, k = min(10, floor(nrow(data_for_spline) / 3))),
                   y ~ s(time, k = min(10, floor(nrow(data_for_spline) / 3)))),
              family = mgcv::mvn(d = 2), data = data_for_spline)
  }, error = function(e) {
    message("Error fitting multivariate GAM: ", e$message)
    return(NULL)
  })
  
  if (is.null(mv_gam)) {
    return(list(features_df = tibble::tibble(), diagnostic_plot = NULL, spatial_plot = NULL))
  }
  
  if (direction == "forward") {
    time_interval <- c(max(data_for_spline$time), max(data_for_spline$time) + window)
    new_time <- seq(from = min(data_for_spline$time), to = max(time_interval), length.out = 200)
    time_ref <- max(data_for_spline$time)
  } else {
    time_interval <- c(min(data_for_spline$time) - window, min(data_for_spline$time))
    new_time <- seq(from = min(time_interval), to = max(data_for_spline$time), length.out = 200)
    time_ref <- min(data_for_spline$time)
  }
  
  xy_pred <- predict(mv_gam, newdata = tibble::tibble(time = new_time), type = "response")
  
  candidate_segments <- get_candidate_segments(df_all, time_interval,
                                                unique(data_for_spline$segment),
                                                min_candidates = min_candidates,
                                                direction = direction)
  
  if (length(candidate_segments) == 0L) {
    return(list(features_df = tibble::tibble(), diagnostic_plot = NULL, spatial_plot = NULL))
  }
  
  if (direction == "forward") {
    candidate_coords <- df_all |>
      dplyr::filter(segment %in% candidate_segments) |>
      dplyr::group_by(segment) |>
      dplyr::arrange(time, .by_group = TRUE) |>
      dplyr::summarise(start_time = dplyr::first(time), start_x = dplyr::first(x), start_y = dplyr::first(y), .groups = "drop")
    
    current_point <- data_for_spline |>
      dplyr::arrange(time) |>
      dplyr::filter(time >= (max(time) - 1)) |>
      dplyr::select(time, x, y, speed, direction) |>
      dplyr::mutate(acceleration = speed - dplyr::lag(speed)) |>
      dplyr::summarise(time = max(time), x = dplyr::last(x), y = dplyr::last(y),
                speed = mean(speed), acceleration = mean(acceleration, na.rm = TRUE),
                direction = mean(direction))
  } else {
    candidate_coords <- df_all |>
      dplyr::filter(segment %in% candidate_segments) |>
      dplyr::group_by(segment) |>
      dplyr::arrange(time, .by_group = TRUE) |>
      dplyr::summarise(end_time = dplyr::last(time), end_x = dplyr::last(x), end_y = dplyr::last(y), .groups = "drop")
    
    current_point <- data_for_spline |>
      dplyr::arrange(time) |>
      dplyr::filter(time <= (min(time) + 1)) |>
      dplyr::select(time, x, y, speed, direction) |>
      dplyr::mutate(acceleration = speed - dplyr::lag(speed)) |>
      dplyr::summarise(time = min(time), x = dplyr::first(x), y = dplyr::first(y),
                speed = mean(speed), acceleration = mean(acceleration, na.rm = TRUE),
                direction = mean(direction))
  }
  
  coord_time <- if (direction == "forward") "start_time" else "end_time"
  coord_x <- if (direction == "forward") "start_x" else "end_x"
  coord_y <- if (direction == "forward") "start_y" else "end_y"
  
  features_df <- candidate_coords |>
    dplyr::rowwise() |>
    dplyr::mutate(
      delta_t = as.numeric(.data[[coord_time]] - current_point$time),
      xy_pred_candidate = list(predict(mv_gam, newdata = data.frame(time = .data[[coord_time]]), type = "response")),
      pred_x = xy_pred_candidate[[1]][1],
      pred_y = xy_pred_candidate[[2]][1],
      residuals_x = .data[[coord_x]] - pred_x,
      residuals_y = .data[[coord_y]] - pred_y,
      residual_gap = sqrt(residuals_x^2 + residuals_y^2),
      delta_x = .data[[coord_x]] - current_point$x,
      delta_y = .data[[coord_y]] - current_point$y,
      spatial_gap = sqrt(delta_x^2 + delta_y^2),
      theoretical_velocity = spatial_gap / abs(delta_t),
      theoretical_acceleration = (theoretical_velocity - current_point$speed) / abs(delta_t),
      sigma_t = compute_uncertainty(abs(delta_t)),
      log_liklihood = calculate_log_liklihood(residual_gap, sigma_t)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-xy_pred_candidate) |>
    dplyr::mutate(
      thr = get_speed_threshold(df_all, prob = speed_threshold_quantile),
      feasible = theoretical_velocity <= thr
    ) |>
    dplyr::arrange(dplyr::desc(feasible), dplyr::desc(log_liklihood)) |>
    dplyr::mutate(
      rank = dplyr::if_else(feasible, dplyr::row_number(), NA_integer_)
    )
  
  if (nrow(features_df) > 0 && all(is.na(features_df$rank))) {
    features_df <- features_df |>
      dplyr::arrange(dplyr::desc(log_liklihood)) |>
      dplyr::mutate(rank = replace(rep(NA_integer_, dplyr::n()), 1, 1))
  }
  
  if (diagnostic_plots) {
    plot_df <- tibble::tibble(
      time = new_time,
      fit_x = as.numeric(xy_pred[, 1]),
      fit_y = as.numeric(xy_pred[, 2]),
      fit_xy = fit_x + fit_y,
      time_gap = abs(time - time_ref),
      uncertainty = ifelse(abs(time - time_ref) > 0, compute_uncertainty(abs(time - time_ref)), 0)
    )
    
    diagnostic_plot <- generate_diagnostic_plot(data_for_spline, plot_df, features_df, 
                                                 current_seg_number, min(df_all$segment), direction)
    
    plot_df_spatial <- tibble::tibble(
      fit_x = as.numeric(xy_pred[, 1]),
      fit_y = as.numeric(xy_pred[, 2]),
      time = new_time,
      time_gap = abs(time - time_ref),
      uncertainty = ifelse(abs(time - time_ref) > 0, compute_uncertainty(abs(time - time_ref)), 0)
    )
    
    spatial_plot <- generate_spatial_plot(data_for_spline, plot_df_spatial, features_df,
                                           current_seg_number, min(df_all$segment), direction)
  }
  
  list(features_df = features_df, diagnostic_plot = diagnostic_plot, spatial_plot = spatial_plot)
}

# 5. BIDIRECTIONAL TRAJECTORY RECONSTRUCTION ====================

#' Propagate Forward from Ground Truth
#'
#' Iteratively link segments forward in time from the ground truth anchor.
#'
#' @param df_all Full data frame with all segments
#' @param ground_truth_segment Ground truth segment ID
#' @param window Search window size
#' @param speed_threshold_quantile Speed feasibility quantile (default: 0.999)
#' @param max_steps Maximum iteration steps (default: 2000)
#' @param save_plots Save diagnostic plots (default: FALSE)
#' @param min_candidates Minimum candidates to consider (default: 5)
#' @return List with track, decisions, plots, and n_groups
#' @export
propagate_forward_uncertainty <- function(df_all,
                                          ground_truth_segment,
                                          window = 30,
                                          speed_threshold_quantile = 0.9999,
                                          max_steps = 2000,
                                          save_plots = FALSE,
                                          min_candidates = 5) {

  df_all <- identify_isolated_groups(df_all, window = window)
  
  n_groups <- max(df_all$tracking_group, na.rm = TRUE)
  
  if (n_groups > 1) {
    cat(sprintf("Found %d isolated tracking groups (gaps > %s seconds)\n", n_groups, window))
  }
  
  start_group <- df_all |>
    dplyr::filter(segment == ground_truth_segment) |>
    dplyr::pull(tracking_group) |>
    unique() |>
    dplyr::first()
  
  all_tracks <- list()
  all_decisions <- list()
  all_plots <- list()
  
  for (group_id in seq_len(n_groups)) {
    group_data <- df_all |> dplyr::filter(tracking_group == group_id)
    
    if (group_id == start_group) {
      group_ground_truth_segment <- ground_truth_segment
    } else {
      group_ground_truth_segment <- tryCatch({
        find_ground_truth_segment(df_all, group_id = group_id)
      }, error = function(e) {
        return(NULL)
      })
    }
    
    if (is.null(group_ground_truth_segment)) {
      cat(sprintf("  Skipping group %d (no valid ground truth)\n", group_id))
      next
    }
    
    cat(sprintf("  Processing group %d: starting from segment %d\n", group_id, group_ground_truth_segment))
    
    compute_uncertainty <- compute_uncertainty_model(group_data, window = window)
    
    current_track <- group_data |> dplyr::filter(segment == group_ground_truth_segment)
    used_segments <- unique(current_track$segment)
    
    history_rows <- list()
    step_plots <- list()
    steps_completed <- 0
    
    for (step in seq_len(max_steps)) {
      res <- rank_candidates_uncertainty(
        current_track,
        df_all = group_data,
        window = window,
        speed_threshold_quantile = speed_threshold_quantile,
        diagnostic_plots = save_plots,
        compute_uncertainty = compute_uncertainty,
        direction = "forward",
        min_candidates = min_candidates
      )
      
      feats <- res$features_df
      
      if (nrow(feats) == 0L) break
      
      feats <- feats |> dplyr::filter(!segment %in% used_segments)
      
      if (nrow(feats) == 0L) break
      
      winner <- feats |>
        dplyr::arrange(is.na(rank), rank, dplyr::desc(log_liklihood)) |>
        dplyr::slice(1)
      
      history_rows[[length(history_rows) + 1]] <-
        feats |>
        dplyr::mutate(step = step,
               chosen = segment == winner$segment,
               tracking_group = group_id)
      
      if (save_plots) {
        step_plots[[length(step_plots) + 1]] <- list(
          diagnostic_plot = res$diagnostic_plot, 
          spatial_plot = res$spatial_plot,
          group_id = group_id
        )
      }
      
      seg_to_add <- winner$segment
      current_track <- dplyr::bind_rows(current_track, group_data |> dplyr::filter(segment == seg_to_add))
      used_segments <- c(used_segments, seg_to_add)
      steps_completed <- step
    }
    
    if (steps_completed == max_steps) {
      warning(sprintf("Forward propagation for group %d reached max_steps (%d). Consider increasing max_steps for more complete reconstruction.", 
                      group_id, max_steps), call. = FALSE)
    }
    
    all_tracks[[group_id]] <- current_track
    if (length(history_rows)) {
      all_decisions[[group_id]] <- dplyr::bind_rows(history_rows)
    }
    all_plots <- c(all_plots, step_plots)
  }
  
  combined_track <- dplyr::bind_rows(all_tracks) |> dplyr::arrange(time)
  combined_decisions <- if (length(all_decisions)) {
    dplyr::bind_rows(all_decisions) |> dplyr::arrange(tracking_group, step, dplyr::desc(chosen), dplyr::desc(log_liklihood))
  } else {
    tibble::tibble()
  }

  list(
    track = combined_track,
    decisions = combined_decisions,
    plots = all_plots,
    n_groups = n_groups
  )
}

#' Propagate Backward from Ground Truth
#'
#' Iteratively link segments backward in time from the ground truth anchor.
#'
#' @param df_all Full data frame with all segments
#' @param ground_truth_segment Ground truth segment ID
#' @param window Search window size
#' @param speed_threshold_quantile Speed feasibility quantile (default: 0.999)
#' @param max_steps Maximum iteration steps (default: 50)
#' @param collect_plots Collect diagnostic plots (default: TRUE)
#' @param min_candidates Minimum candidates to consider (default: 5)
#' @return List with track, decisions, plots, and n_groups
#' @export
propagate_backwards_uncertainty <- function(df_all,
                                            ground_truth_segment,
                                            window = 30,
                                            speed_threshold_quantile = 0.9999,
                                            max_steps = 2000,
                                            collect_plots = TRUE,
                                            min_candidates = 5) {
  
  df_all <- identify_isolated_groups(df_all, window = window)
  
  n_groups <- max(df_all$tracking_group, na.rm = TRUE)
  
  if (n_groups > 1) {
    cat(sprintf("Found %d isolated tracking groups (gaps > %s seconds)\n", n_groups, window))
  }
  
  start_group <- df_all |>
    dplyr::filter(segment == ground_truth_segment) |>
    dplyr::pull(tracking_group) |>
    unique() |>
    dplyr::first()
  
  all_tracks <- list()
  all_decisions <- list()
  all_plots <- list()
  
  for (group_id in seq_len(n_groups)) {
    group_data <- df_all |> dplyr::filter(tracking_group == group_id)
    
    if (group_id == start_group) {
      group_ground_truth_segment <- ground_truth_segment
    } else {
      group_ground_truth_segment <- tryCatch({
        find_ground_truth_segment(df_all, group_id = group_id)
      }, error = function(e) {
        return(NULL)
      })
    }
    
    if (is.null(group_ground_truth_segment)) {
      cat(sprintf("  Skipping group %d (no valid ground truth)\n", group_id))
      next
    }
    
    cat(sprintf("  Processing group %d backwards: starting from segment %d\n", group_id, group_ground_truth_segment))
    
    compute_uncertainty <- compute_uncertainty_model(group_data, window = window)
    
    current_track <- group_data |> dplyr::filter(segment == group_ground_truth_segment)
    decision_log <- list()
    per_step_plots <- list()
    step_counter <- 0
    
    repeat {
      step_counter <- step_counter + 1
      
      if (step_counter > max_steps) {
        warning(sprintf("Backward propagation for group %d reached max_steps (%d). Consider increasing max_steps for more complete reconstruction.", 
                        group_id, max_steps), call. = FALSE)
        break
      }
      
      res <- rank_candidates_uncertainty(
        current_track = current_track,
        df_all = group_data,
        window = window,
        speed_threshold_quantile = speed_threshold_quantile,
        diagnostic_plots = collect_plots,
        compute_uncertainty = compute_uncertainty,
        direction = "backward",
        min_candidates = min_candidates
      )
      
      candidates <- res$features_df
      
      if (nrow(candidates) == 0) {
        message("No more candidates at step ", step_counter, " for group ", group_id)
        break
      }
      
      winner <- candidates |>
        dplyr::filter(!is.na(rank)) |>
        dplyr::arrange(rank) |>
        dplyr::slice(1)
      
      if (nrow(winner) == 0) {
        message("No feasible candidates at step ", step_counter, " for group ", group_id)
        break
      }
      
      winner_segment <- winner$segment
      
      decision_log[[step_counter]] <- tibble::tibble(
        step = step_counter,
        chosen_segment = winner_segment,
        rank = winner$rank,
        log_liklihood = winner$log_liklihood,
        residual_gap = winner$residual_gap,
        delta_t = winner$delta_t,
        sigma_t = winner$sigma_t,
        theoretical_velocity = winner$theoretical_velocity,
        tracking_group = group_id
      )
      
      if (collect_plots) {
        per_step_plots[[length(per_step_plots) + 1]] <- list(
          diagnostic_plot = res$diagnostic_plot,
          spatial_plot = res$spatial_plot,
          group_id = group_id
        )
      }
      
      new_segment_data <- group_data |> dplyr::filter(segment == winner_segment)
      current_track <- dplyr::bind_rows(new_segment_data, current_track)
      
      message("  Step ", step_counter, ": added segment ", winner_segment, 
              " (prob = ", round(winner$log_liklihood, 3), ")")
    }
    
    all_tracks[[group_id]] <- current_track
    if (length(decision_log)) {
      all_decisions[[group_id]] <- dplyr::bind_rows(decision_log)
    }
    all_plots <- c(all_plots, per_step_plots)
  }
  
  combined_track <- dplyr::bind_rows(all_tracks) |> dplyr::arrange(time)
  combined_decisions <- if (length(all_decisions)) {
    dplyr::bind_rows(all_decisions) |> dplyr::arrange(tracking_group, step)
  } else {
    tibble::tibble()
  }
  
  list(
    track = combined_track,
    decisions = combined_decisions,
    plots = all_plots,
    n_groups = n_groups
  )
}

# 6. POST-PROCESSING ============================================

#' Smooth Reconstructed Trajectory
#'
#' Apply Savitzky-Golay or moving average filter to reduce noise.
#'
#' @param df Data frame with x_reconstructed and y_reconstructed columns
#' @param method Either "savitzky_golay" or "moving_average"
#' @param n Window size (default: "not_specified" for auto-calculation)
#' @param p Polynomial order for Savitzky-Golay (default: 3)
#' @return Data frame with x_smooth and y_smooth columns added
#' @export
smooth_track <- function(df, method = "savitzky_golay", n = "not_specified", p = 3) {
  if (n == "not_specified") {
    time_diffs <- diff(df$time[!is.na(df$time)])
    median_dt <- stats::median(time_diffs, na.rm = TRUE)
    fps <- 1 / median_dt
    n <- round(fps / 2)
    if (n %% 2 == 0) n <- n + 1
  }
  
  valid_data <- df |> dplyr::filter(!is.na(x_reconstructed), !is.na(y_reconstructed))
  if (nrow(valid_data) < n) {
    warning("Not enough reconstructed data points to smooth")
    return(df |> dplyr::mutate(x_smooth = NA_real_, y_smooth = NA_real_))
  }
  
  df_smoothed <- df |>
    dplyr::mutate(
      x_filled = zoo::na.approx(x_reconstructed, na.rm = FALSE),
      y_filled = zoo::na.approx(y_reconstructed, na.rm = FALSE)
    ) |>
    dplyr::filter(!is.na(x_filled), !is.na(y_filled))
  
  if (method == "savitzky_golay") {
    df_smoothed <- df_smoothed |>
      dplyr::mutate(
        x_smooth = signal::sgolayfilt(x_filled, p = p, n = n),
        y_smooth = signal::sgolayfilt(y_filled, p = p, n = n)
      )
  } else if (method == "moving_average") {
    df_smoothed <- df_smoothed |>
      dplyr::mutate(
        x_smooth = zoo::rollmean(x_filled, k = n, fill = NA, align = "center"),
        y_smooth = zoo::rollmean(y_filled, k = n, fill = NA, align = "center")
      )
  } else {
    stop(sprintf("Unknown smoothing method '%s'. Use 'savitzky_golay' or 'moving_average'", method))
  }
  
  df <- dplyr::left_join(df, df_smoothed |> dplyr::select(time, x_smooth, y_smooth), by = "time")
  return(df)
}

# 7. MULTI-INDIVIDUAL PIPELINE ==================================

#' Process Multiple Tracked Individuals
#'
#' Complete reconstruction pipeline for multi-individual tracking data.
#'
#' @param df_raw Data frame with columns: time, x1, y1, x2, y2, etc.
#' @param window Search window size (default: 20)
#' @param speed_threshold_quantile Speed feasibility quantile (default: 0.999)
#' @param smooth_method Smoothing method (default: "savitzky_golay")
#' @param smooth_window Smoothing window size (default: "not_specified")
#' @param min_movement Minimum movement for ground truth (default: 500)
#' @param diagnostic_plots Generate diagnostic plots (default: FALSE)
#' @param min_candidates Minimum candidates to consider (default: 5)
#' @return List with tracks data frame and optionally plots list
#' @export
process_multi_individual <- function(df_raw,
                                     window = 20,
                                     speed_threshold_quantile = 0.9999,
                                     smooth_method = "savitzky_golay",
                                     smooth_window = "not_specified",
                                     min_movement = 500,
                                     diagnostic_plots = FALSE,
                                     min_candidates = 5) {
  x_cols <- grep("^x\\d+$", names(df_raw), value = TRUE)
  n_individuals <- length(x_cols)
  
  if (n_individuals == 0) {
    stop("No individual columns found. Expected columns like x1, y1, x2, y2, etc.")
  }
  
  cat(sprintf("Processing %d individual(s)...\n", n_individuals))
  
  all_tracks <- list()
  all_plots <- list()
  
  for (i in 1:n_individuals) {
    cat(sprintf("\n=== Processing Individual %d/%d ===\n", i, n_individuals))
    
    x_col <- paste0("x", i)
    y_col <- paste0("y", i)
    
    df_individual <- df_raw |>
      dplyr::select(time, x = dplyr::all_of(x_col), y = dplyr::all_of(y_col))
    
    cat("  Step 1/6: Detecting jumps and segmenting...\n")
    df_jump <- detect_jumps(df_individual)
    df_seg <- assign_segments(df_jump)
    
    gaps <- identify_tracking_gaps(df_seg, window = window)
    if (nrow(gaps) > 0) {
      cat(sprintf("  Found %d large tracking gaps (> %s seconds)\n", nrow(gaps), window))
    }
    
    cat("  Step 2/6: Finding ground truth segment...\n")
    ground_truth_segment <- tryCatch({
      find_ground_truth_segment(df_seg, min_movement = min_movement)
    }, error = function(e) {
      warning(sprintf("Individual %d: %s", i, e$message))
      return(NULL)
    })
    
    if (is.null(ground_truth_segment)) {
      next
    }
    
    cat(sprintf("  Ground truth segment: %d\n", ground_truth_segment))
    
    cat("  Step 3/6: Forward propagation (gap-aware)...\n")
    res_forward <- suppressWarnings(propagate_forward_uncertainty(
      df_all = df_seg,
      ground_truth_segment = ground_truth_segment,
      window = window,
      speed_threshold_quantile = speed_threshold_quantile,
      max_steps = 2000,
      save_plots = diagnostic_plots,
      min_candidates = min_candidates
    ))
    final_track_forward <- res_forward$track
    
    cat("  Step 4/6: Backward propagation (gap-aware)...\n")
    res_backward <- suppressWarnings(suppressMessages(propagate_backwards_uncertainty(
      df_all = df_seg,
      ground_truth_segment = ground_truth_segment,
      window = window,
      speed_threshold_quantile = speed_threshold_quantile,
      max_steps = 2000,
      collect_plots = diagnostic_plots,
      min_candidates = min_candidates
    )))
    final_track_backward <- res_backward$track
    
    if (diagnostic_plots) {
      all_plots[[i]] <- list(
        backward = res_backward$plots,
        forward = res_forward$plots
      )
    }
    
    cat("  Step 5/6: Combining tracks...\n")
    complete_track <- dplyr::bind_rows(
      final_track_backward |> dplyr::filter(segment != ground_truth_segment),
      df_seg |> dplyr::filter(segment == ground_truth_segment),
      final_track_forward |> dplyr::filter(segment != ground_truth_segment)
    ) |>
      dplyr::arrange(time)
    
    output_track <- df_individual |>
      dplyr::select(time, x, y) |>
      dplyr::rename(x_original = x, y_original = y) |>
      dplyr::left_join(complete_track |> dplyr::select(time, x, y, segment), 
                by = "time",
                suffix = c("_original", "_reconstructed")) |>
      dplyr::rename(x_reconstructed = x, y_reconstructed = y)
    
    cat(sprintf("  Step 6/6: Smoothing (%s)...\n", smooth_method))
    output_track <- tryCatch({
      smooth_track(output_track, method = smooth_method, n = smooth_window)
    }, error = function(e) {
      warning(sprintf("Individual %d smoothing failed: %s", i, e$message))
      output_track |> dplyr::mutate(x_smooth = NA_real_, y_smooth = NA_real_)
    })
    
    output_track <- output_track |>
      dplyr::mutate(individual_number = i)
    
    all_tracks[[i]] <- output_track
    
    cat(sprintf("  ✓ Individual %d complete (%d points)\n", i, nrow(output_track)))
  }
  
  if (length(all_tracks) == 0) {
    stop("No individuals were successfully processed")
  }
  
  cat(sprintf("\n=== Combining %d individual(s) ===\n", length(all_tracks)))
  combined_tracks <- dplyr::bind_rows(all_tracks)
  
  cat(sprintf("✓ Multi-individual processing complete: %d total points\n", nrow(combined_tracks)))
  
  result <- list(tracks = combined_tracks)
  
  if (diagnostic_plots) {
    result$plots <- all_plots
  }
  
  return(result)
}

#' Smooth Multiple Tracked Individuals Without Reconstruction
#'
#' Applies smoothing to multi-individual tracking data without performing
#' trajectory reconstruction. This is useful when the raw tracking data is
#' already clean and only needs smoothing, or when you want to compare
#' smoothed raw data against reconstructed data.
#'
#' @param df_raw Data frame with columns: time, x1, y1, x2, y2, etc.
#' @param smooth_method Smoothing method: "savitzky_golay" or "moving_average" (default: "savitzky_golay")
#' @param smooth_window Smoothing window size. If "not_specified", automatically calculated from frame rate.
#' @return List with tracks data frame compatible with overlay_track_on_video()
#' @export
#' @examples
#' \dontrun{
#' # Load tracking data
#' df <- read.csv("tracking_data.csv")
#' 
#' # Smooth without reconstruction
#' result <- smooth_multi_individual(df)
#' 
#' # Create video overlay
#' overlay_track_on_video(result$tracks, "video.mp4")
#' }
smooth_multi_individual <- function(df_raw,
                                    smooth_method = "savitzky_golay",
                                    smooth_window = "not_specified") {
  # Ensure df_raw is a data frame
  if (!is.data.frame(df_raw)) {
    stop("df_raw must be a data frame")
  }
  
  # Check for time column
  if (!"time" %in% names(df_raw)) {
    stop("df_raw must have a 'time' column")
  }
  
  x_cols <- grep("^x\\d+$", names(df_raw), value = TRUE)
  n_individuals <- length(x_cols)
  
  if (n_individuals == 0) {
    stop("No individual columns found. Expected columns like x1, y1, x2, y2, etc.")
  }
  
  cat(sprintf("Smoothing %d individual(s) (no reconstruction)...\n", n_individuals))
  
  all_tracks <- list()
  
  for (i in 1:n_individuals) {
    cat(sprintf("\n=== Smoothing Individual %d/%d ===\n", i, n_individuals))
    
    x_col <- paste0("x", i)
    y_col <- paste0("y", i)
    
    # Check that columns exist
    if (!x_col %in% names(df_raw) || !y_col %in% names(df_raw)) {
      warning(sprintf("Individual %d: Missing columns %s or %s, skipping", i, x_col, y_col))
      next
    }
    
    # Extract individual data
    df_individual <- data.frame(
      time = df_raw$time,
      x = df_raw[[x_col]],
      y = df_raw[[y_col]]
    )
    
    # Create output structure with original data as "reconstructed" for compatibility
    output_track <- data.frame(
      time = df_individual$time,
      x_original = df_individual$x,
      y_original = df_individual$y,
      x_reconstructed = df_individual$x,
      y_reconstructed = df_individual$y,
      segment = 1L
    )
    
    cat(sprintf("  Smoothing (%s)...\n", smooth_method))
    
    # Calculate smooth window if not specified
    n <- smooth_window
    if (identical(n, "not_specified")) {
      time_vals <- output_track$time[!is.na(output_track$time)]
      if (length(time_vals) > 1) {
        time_diffs <- diff(time_vals)
        median_dt <- stats::median(time_diffs, na.rm = TRUE)
        if (!is.na(median_dt) && median_dt > 0) {
          fps <- 1 / median_dt
          n <- round(fps / 2)
          if (n %% 2 == 0) n <- n + 1
          if (n < 3) n <- 3  # Minimum window size
        } else {
          n <- 15  # Default fallback
        }
      } else {
        n <- 15  # Default fallback
      }
    }
    
    # Count valid data points
    valid_mask <- !is.na(output_track$x_reconstructed) & !is.na(output_track$y_reconstructed)
    n_valid <- sum(valid_mask)
    
    if (n_valid < n) {
      warning(sprintf("Individual %d: Not enough data points to smooth (%d valid, need %d)", i, n_valid, n))
      output_track$x_smooth <- NA_real_
      output_track$y_smooth <- NA_real_
    } else {
      # Interpolate missing values
      x_filled <- zoo::na.approx(output_track$x_reconstructed, na.rm = FALSE)
      y_filled <- zoo::na.approx(output_track$y_reconstructed, na.rm = FALSE)
      
      # Create smoothed data frame with filled values
      filled_mask <- !is.na(x_filled) & !is.na(y_filled)
      
      if (sum(filled_mask) < n) {
        warning(sprintf("Individual %d: Not enough filled data points to smooth", i))
        output_track$x_smooth <- NA_real_
        output_track$y_smooth <- NA_real_
      } else {
        # Extract filled values
        times_filled <- output_track$time[filled_mask]
        x_filled_vals <- x_filled[filled_mask]
        y_filled_vals <- y_filled[filled_mask]
        
        # Apply smoothing
        if (smooth_method == "savitzky_golay") {
          x_smooth_vals <- signal::sgolayfilt(x_filled_vals, p = 3, n = n)
          y_smooth_vals <- signal::sgolayfilt(y_filled_vals, p = 3, n = n)
        } else if (smooth_method == "moving_average") {
          x_smooth_vals <- zoo::rollmean(x_filled_vals, k = n, fill = NA, align = "center")
          y_smooth_vals <- zoo::rollmean(y_filled_vals, k = n, fill = NA, align = "center")
        } else {
          stop(sprintf("Unknown smoothing method '%s'. Use 'savitzky_golay' or 'moving_average'", smooth_method))
        }
        
        # Create lookup for smoothed values
        smooth_df <- data.frame(
          time = times_filled,
          x_smooth = x_smooth_vals,
          y_smooth = y_smooth_vals
        )
        
        # Merge back
        output_track <- merge(output_track, smooth_df, by = "time", all.x = TRUE)
      }
    }
    
    output_track$individual_number <- i
    
    all_tracks[[i]] <- output_track
    
    cat(sprintf("  ✓ Individual %d complete (%d points)\n", i, nrow(output_track)))
  }
  
  # Filter out NULL entries (skipped individuals)
  all_tracks <- all_tracks[!sapply(all_tracks, is.null)]
  
  if (length(all_tracks) == 0) {
    stop("No individuals were successfully processed")
  }
  
  cat(sprintf("\n=== Combining %d individual(s) ===\n", length(all_tracks)))
  combined_tracks <- do.call(rbind, all_tracks)
  rownames(combined_tracks) <- NULL
  
  cat(sprintf("✓ Multi-individual smoothing complete: %d total points\n", nrow(combined_tracks)))
  
  result <- list(tracks = combined_tracks)
  
  return(result)
}

# 8. VIDEO OVERLAY ==============================================

#' Convert Seconds to ASS Time Format
#'
#' @param t Time in seconds
#' @return Character string in ASS format (H:MM:SS.CS)
#' @export
seconds_to_ass_time <- function(t) {
  hours <- floor(t / 3600)
  remainder <- t %% 3600
  minutes <- floor(remainder / 60)
  seconds <- floor(remainder %% 60)
  centiseconds <- round((t %% 1) * 100)
  sprintf("%d:%02d:%02d.%02d", hours, minutes, seconds, centiseconds)
}

#' Get Video Resolution Using ffprobe
#'
#' @param video_path Path to video file
#' @return List with width and height
#' @export
get_video_resolution <- function(video_path) {
  cmd_res <- sprintf(
    'ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of csv=p=0 "%s"',
    video_path
  )
  result_res <- system(cmd_res, intern = TRUE)
  if (length(result_res) != 1) stop("Failed to get video resolution")
  
  dimensions <- strsplit(result_res, ",")[[1]]
  width <- as.integer(dimensions[1])
  height <- as.integer(dimensions[2])
  
  return(list(width = width, height = height))
}

#' Create Fish Tracking Overlay on Video
#'
#' Generates a video with reconstructed tracks overlaid on original footage.
#'
#' @param tracks_df Data frame with tracking results
#' @param video_path Path to input video file
#' @param output_name Output video filename (default: "overlaid.mp4")
#' @param output_location Optional directory path for output video. 
#' @param dot_size Size of tracking dots (default: 18)
#' @param use_smoothed Use smoothed track if available (default: TRUE)
#' @return Path to output video
#' @export
overlay_track_on_video <- function(tracks_df, video_path, 
                                   output_name = "overlaid.mp4",
                                   output_location = NULL,
                                   dot_size = 18,
                                   use_smoothed = TRUE) {
  
  # Convert to Path object and validate
  video_path <- normalizePath(video_path, mustWork = TRUE)
  video_dir <- dirname(video_path)
  
  # Determine output path based on output_location parameter
  if (is.null(output_location)) {
    output_path <- file.path(video_dir, output_name)
  } else {
    # Validate output_location exists
    if (!dir.exists(output_location)) {
      stop(sprintf("Output location does not exist: %s", output_location))
    }
    output_path <- file.path(output_location, output_name)
  }
  
  ass_path <- file.path(video_dir, "overlay.ass")
  
  # Define colour palette for individuals (BBGGRR format)
  # Each individual gets a bright colour for reconstructed and dark colour for original
  colour_palette_bright <- c(
    "&H0000FF&",  # Bright Red
    "&H00FF00&",  # Bright Green
    "&HFF0000&",  # Bright Blue
    "&H00FFFF&",  # Bright Yellow
    "&HFF00FF&",  # Bright Magenta
    "&HFFFF00&",  # Bright Cyan
    "&H0080FF&",  # Bright Orange
    "&HFF0080&",  # Bright Purple
    "&H00FF80&",  # Bright Spring green
    "&H80FF00&"   # Bright Chartreuse
  )
  
  colour_palette_dark <- c(
    "&H000080&",  # Dark Red
    "&H008000&",  # Dark Green
    "&H800000&",  # Dark Blue
    "&H008080&",  # Dark Yellow/Olive
    "&H800080&",  # Dark Magenta
    "&H808000&",  # Dark Cyan/Teal
    "&H004080&",  # Dark Orange
    "&H800040&",  # Dark Purple
    "&H008040&",  # Dark Spring green
    "&H408000&"   # Dark Chartreuse
  )
  
  # Helper: Convert seconds to ASS timestamp format (H:MM:SS.CS)
  seconds_to_ass_time <- function(t) {
    hours <- floor(t / 3600)
    remainder <- t %% 3600
    minutes <- floor(remainder / 60)
    seconds <- floor(remainder %% 60)
    centiseconds <- round((t - floor(t)) * 100)
    sprintf("%d:%02d:%02d.%02d", hours, minutes, seconds, centiseconds)
  }
  
  # Get video properties using ffprobe
  get_video_properties <- function(video_path) {
    # Get resolution
    cmd_res <- sprintf(
      'ffprobe -v error -select_streams v:0 -show_entries stream=width,height -of default=noprint_wrappers=1:nokey=1 "%s"',
      video_path
    )
    result_res <- system(cmd_res, intern = TRUE)
    if (length(result_res) != 2) stop("Failed to get video resolution")
    
    # Get frame rate (r_frame_rate returns as fraction like "30/1")
    cmd_fps <- sprintf(
      'ffprobe -v error -select_streams v:0 -show_entries stream=r_frame_rate -of default=noprint_wrappers=1:nokey=1 "%s"',
      video_path
    )
    result_fps <- system(cmd_fps, intern = TRUE)
    
    # Parse frame rate (handle fraction format like "30/1" or "30000/1001")
    fps_parts <- strsplit(result_fps, "/")[[1]]
    if (length(fps_parts) == 2) {
      fps <- as.numeric(fps_parts[1]) / as.numeric(fps_parts[2])
    } else {
      fps <- as.numeric(result_fps)
    }
    
    list(
      width = as.integer(result_res[1]), 
      height = as.integer(result_res[2]),
      frame_rate = fps
    )
  }
  
  cat("Getting video properties...\n")
  video_props <- get_video_properties(video_path)
  cat(sprintf("  Resolution: %dx%d\n", video_props$width, video_props$height))
  cat(sprintf("  Frame rate: %.2f fps\n", video_props$frame_rate))
  
  # Check for individual_number column
  if (!"individual_number" %in% names(tracks_df)) {
    stop("tracks_df must have an 'individual_number' column")
  }
  
  n_individuals <- max(tracks_df$individual_number, na.rm = TRUE)
  cat(sprintf("Processing %d individual(s)...\n", n_individuals))
  
  # Determine which reconstructed track to use
  recon_x_col <- if (use_smoothed && "x_smooth" %in% names(tracks_df)) "x_smooth" else "x_reconstructed"
  recon_y_col <- if (use_smoothed && "y_smooth" %in% names(tracks_df)) "y_smooth" else "y_reconstructed"
  
  track_clean <- tracks_df |>
    dplyr::filter(!is.na(time), !is.na(individual_number)) |>
    dplyr::arrange(time, individual_number) |>
    dplyr::rename(
      recon_x = dplyr::all_of(recon_x_col),
      recon_y = dplyr::all_of(recon_y_col)
    )
  
  if (nrow(track_clean) == 0) {
    stop("No valid track data (all time values are NA)")
  }
  
  max_time <- max(track_clean$time)
  
  cat(sprintf("Generating ASS subtitle file (using %s track)...\n", 
              if(use_smoothed && "x_smooth" %in% names(tracks_df)) "smoothed" else "reconstructed"))
  
  # ASS file header
  ass_header <- c(
    "[Script Info]",
    "ScriptType: v4.00+",
    sprintf("PlayResX: %d", video_props$width),
    sprintf("PlayResY: %d", video_props$height),
    "",
    "[V4+ Styles]",
    "Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding",
    sprintf("Style: Default,Arial,%d,&H00FFFFFF,&H000000FF,&H00000000,&H00000000,0,0,0,0,100,100,0,0,1,1,0,7,0,0,0,0", dot_size),
    "",
    "[Events]",
    "Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text"
  )
  
  # **OPTIMIZED: Vectorized dialogue generation**
  cat(sprintf("Generating %d dialogue lines...\n", nrow(track_clean) * 2))
  
  # Pre-compute all time values
  frame_duration <- 1 / video_props$frame_rate
  start_times <- vapply(track_clean$time, seconds_to_ass_time, character(1))
  end_times <- vapply(track_clean$time + frame_duration, seconds_to_ass_time, character(1))
  
  # Pre-compute colours for each individual
  ind_bright <- colour_palette_bright[((track_clean$individual_number - 1) %% length(colour_palette_bright)) + 1]
  ind_dark <- colour_palette_dark[((track_clean$individual_number - 1) %% length(colour_palette_dark)) + 1]
  
  # Vectorized generation of original track lines
  original_mask <- !is.na(track_clean$x_original) & !is.na(track_clean$y_original)
  ass_events_original <- sprintf(
    "Dialogue: 0,%s,%s,Default,,0,0,0,,{\\an5\\pos(%.1f,%.1f)\\fs%d\\c%s}●",
    start_times[original_mask],
    end_times[original_mask],
    track_clean$x_original[original_mask],
    track_clean$y_original[original_mask],
    dot_size,
    ind_dark[original_mask]
  )
  
  # Vectorized generation of reconstructed track lines
  recon_mask <- !is.na(track_clean$recon_x) & !is.na(track_clean$recon_y)
  ass_events_recon <- sprintf(
    "Dialogue: 0,%s,%s,Default,,0,0,0,,{\\an5\\pos(%.1f,%.1f)\\fs%d\\c%s}●",
    start_times[recon_mask],
    end_times[recon_mask],
    track_clean$recon_x[recon_mask],
    track_clean$recon_y[recon_mask],
    dot_size,
    ind_bright[recon_mask]
  )
  
  # Combine all events
  ass_events <- c(ass_events_original, ass_events_recon)
  
  # Write ASS file (single write operation)
  cat("Writing ASS file...\n")
  writeLines(c(ass_header, ass_events), ass_path)
  cat(sprintf("ASS file created: %s (%d lines)\n", ass_path, length(ass_events)))
  
  # Run ffmpeg to overlay subtitles
  cat("Running ffmpeg to create overlaid video...\n")
  ffmpeg_cmd <- sprintf(
    'ffmpeg -y -i "%s" -vf "subtitles=%s" -c:a copy -c:v libx264 -preset fast -crf 23 -t %.2f "%s"',
    video_path, ass_path, max_time, output_path
  )
  
  result <- system(ffmpeg_cmd)
  
  # Clean up ASS file
  unlink(ass_path)
  
  if (result == 0) {
    cat(sprintf("\n✓ Overlay complete: %s\n", output_path))
    cat(sprintf("  %d individuals overlaid with unique colours\n", n_individuals))
    cat(sprintf("  Original tracks: dark shade | Reconstructed: bright shade\n"))
    return(invisible(output_path))
  } else {
    stop("ffmpeg failed to create overlaid video")
  }
}

# 9. INTERACTIVE PLOT VIEWER ====================================

#' Interactive Diagnostic Plot Viewer
#'
#' Launch a Shiny gadget to browse diagnostic plots chronologically.
#'
#' @param results Output from process_multi_individual with diagnostic_plots=TRUE
#' @param individual_number Specific individual or "all" (default: 1)
#' @return Invisible index of last viewed plot
#' @export
check_spline_fits <- function(results, individual_number = 1) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for interactive plot viewing")
  }
  if (!requireNamespace("miniUI", quietly = TRUE)) {
    stop("Package 'miniUI' is required for interactive plot viewing")
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plot viewing")
  }
  
  if (!is.list(results) || !"plots" %in% names(results)) {
    stop("Input must be the output from process_multi_individual() with diagnostic_plots = TRUE")
  }
  
  if (is.null(results$plots) || length(results$plots) == 0) {
    stop("No plots found. Make sure diagnostic_plots = TRUE when running process_multi_individual()")
  }
  
  if (identical(individual_number, "all")) {
    all_plots <- unlist(lapply(results$plots, function(ind_plots) {
      if (is.null(ind_plots)) return(list())
      c(ind_plots$backward, ind_plots$forward)
    }), recursive = FALSE)
    
    # Filter out NULL plots
    all_plots <- all_plots[!sapply(all_plots, function(p) {
      is.null(p) || (is.null(p$diagnostic_plot) && is.null(p$spatial_plot))
    })]
    
    if (length(all_plots) == 0) {
      stop("No diagnostic or spatial plots available. Plots may be NULL due to insufficient data for uncertainty modelling.")
    }
    
    plots <- all_plots
  } else {
    if (!is.numeric(individual_number) || individual_number < 1 ||
        individual_number > length(results$plots)) {
      stop(sprintf("Invalid individual_number. Must be between 1 and %d", length(results$plots)))
    }
    
    ind_plots <- results$plots[[individual_number]]
    if (is.null(ind_plots)) {
      stop(sprintf("No plots found for individual %d", individual_number))
    }
    
    plots <- c(ind_plots$backward, ind_plots$forward)
    
    # Filter out NULL plots
    plots <- plots[!sapply(plots, function(p) {
      is.null(p) || (is.null(p$diagnostic_plot) && is.null(p$spatial_plot))
    })]
    
    if (length(plots) == 0) {
      stop(sprintf("No diagnostic or spatial plots available for individual %d. Plots may be NULL due to insufficient data for uncertainty modelling.", individual_number))
    }
  }
  
  # Sort plots by segment number from title
  plots <- plots[order(sapply(plots, function(p) {
    # Try to extract segment number from spatial plot title
    if (!is.null(p$spatial_plot) && !is.null(p$spatial_plot$labels$title)) {
      title <- p$spatial_plot$labels$title
      # Extract segment number from "Segment: XX"
      seg_match <- regmatches(title, regexpr("Segment: [0-9]+", title))
      if (length(seg_match) > 0) {
        seg_num <- as.numeric(sub("Segment: ", "", seg_match))
        return(seg_num)
      }
    }
    # Try diagnostic plot title as fallback
    if (!is.null(p$diagnostic_plot) && !is.null(p$diagnostic_plot$labels$title)) {
      title <- p$diagnostic_plot$labels$title
      seg_match <- regmatches(title, regexpr("Segment: [0-9]+", title))
      if (length(seg_match) > 0) {
        seg_num <- as.numeric(sub("Segment: ", "", seg_match))
        return(seg_num)
      }
    }
    return(Inf)
  }))]
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Diagnostic Plot Viewer (use ← and → arrow keys)"),
    miniUI::miniContentPanel(
      shiny::tags$script(shiny::HTML("
        $(document).on('keydown', function(e) {
          if(e.which === 37) Shiny.setInputValue('nav', 'prev', {priority: 'event'});
          if(e.which === 39) Shiny.setInputValue('nav', 'next', {priority: 'event'});
        });
      ")),
      shiny::fluidRow(
        shiny::column(12, align = "centre",
                     shiny::actionButton("btn_prev", "← Prev"),
                     shiny::htmlOutput("idx_label", inline = TRUE),
                     shiny::actionButton("btn_next", "Next →")
        )
      ),
      shiny::fluidRow(
        shiny::column(6, shiny::uiOutput("plot_container_diag")),
        shiny::column(6, shiny::uiOutput("plot_container_spatial"))
      )
    )
  )
  
  server <- function(input, output, session) {
    n <- length(plots)
    i <- shiny::reactiveVal(1L)
    step <- function(d) i(max(1L, min(n, i() + d)))
    
    shiny::observeEvent(input$btn_prev, step(-1))
    shiny::observeEvent(input$btn_next, step(+1))
    shiny::observeEvent(input$nav, ignoreInit = TRUE, {
      if (identical(input$nav, "prev")) step(-1)
      else if (identical(input$nav, "next")) step(+1)
    })
    
    output$idx_label <- shiny::renderText(sprintf("<b>%d / %d</b>", i(), n))
    output$plot_container_diag <- shiny::renderUI(plotly::plotlyOutput("p_diag", height = "700px"))
    output$plot_container_spatial <- shiny::renderUI(plotly::plotlyOutput("p_spatial", height = "700px"))
    
    output$p_diag <- plotly::renderPlotly({
      p <- plots[[i()]]$diagnostic_plot
      if (is.null(p)) {
        plotly::plot_ly() |> plotly::layout(title = "No diagnostic plot")
      } else {
        plotly::ggplotly(p)
      }
    })
    
    output$p_spatial <- plotly::renderPlotly({
      p <- plots[[i()]]$spatial_plot
      if (is.null(p)) {
        plotly::plot_ly() |> plotly::layout(title = "No spatial plot")
      } else {
        plotly::ggplotly(p)
      }
    })
    
    shiny::observeEvent(input$done, { shiny::stopApp(invisible(i())) })
  }
  
  shiny::runGadget(
    ui, server,
    viewer = shiny::dialogViewer("Diagnostic Plot Viewer", width = 1800, height = 900)
  )
}
