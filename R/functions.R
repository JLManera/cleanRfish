#' Kalman Filter Path Connection Method
#'
#' Uses Kalman filtering to predict and connect path segments with probabilistic ranking.
#' If the current segment has fewer than the minimum required points, this function
#' will incorporate points from previously included segments to ensure robust predictions.
#'
#' @param segments_summary Summary of path segments with start/end positions and times
#' @param processed_df The processed tracking data with segment IDs
#' @param time_bias Weight applied to temporal distance in path connection (higher values
#'        penalize larger time gaps more heavily)
#' @param projection_dist Maximum number of future segments to consider in the projection
#' @param process_noise Process noise parameter (controls model uncertainty)
#' @param measurement_noise Measurement noise parameter (controls observation uncertainty)
#' @param min_points Minimum number of points to use for Kalman filtering (defaults to 100)
#'
#' @return Vector of segment IDs that form the connected path
#' @importFrom FKF fkf
connect_kalman_projection <- function(segments_summary, processed_df,
                                      time_bias = 10, projection_dist = 10,
                                      process_noise = 0.1, measurement_noise = 1.0,
                                      min_points = 100) {
  # Check if FKF package is installed and load it
  if (!requireNamespace("FKF", quietly = TRUE)) {
    stop("Package 'FKF' is needed for Kalman filtering. Please install it with install.packages('FKF')")
  }

  filtered_ids <- c(1L)
  current_id <- 1L

  # Keep track of end time for delta_t calculations
  current_end_time <- NULL

  repeat {
    # Get current segment end position
    current_end <- segments_summary |>
      dplyr::filter(segment_id == current_id) |>
      dplyr::select(end_x, end_y, end_time)

    current_end_time <- current_end$end_time

    # Get current segment data
    current_segment <- processed_df |>
      dplyr::filter(as.integer(segment_id) == current_id) |>
      dplyr::arrange(time)

    # NEW: Check if we need to incorporate points from previous segments
    current_segment_points <- nrow(current_segment)
    if (current_segment_points < min_points && length(filtered_ids) > 1) {
      # Calculate how many additional points we need
      points_needed <- min_points - current_segment_points

      # Get points from previous segments in the filtered path
      previous_segments <- processed_df |>
        dplyr::filter(as.integer(segment_id) %in% filtered_ids[-length(filtered_ids)]) |>
        dplyr::arrange(dplyr::desc(time)) |>
        dplyr::slice_head(n = points_needed)

      if (nrow(previous_segments) > 0) {
        # Combine previous points with current segment
        current_segment <- dplyr::bind_rows(previous_segments, current_segment) |>
          dplyr::arrange(time)
      }
    }

    # If there's only one point in the segment, simple projection is used
    if (nrow(current_segment) < 2) {
      med_vel <- stats::median(processed_df$velocity, na.rm = TRUE)

      project_position <- function(future_time) {
        delta_t <- future_time - current_end_time
        pred_x <- current_end$end_x + med_vel * delta_t
        pred_y <- current_end$end_y
        pred_vx <- med_vel
        pred_vy <- 0

        # Return prediction with high uncertainty
        list(
          mean = c(pred_x, pred_y, pred_vx, pred_vy),
          cov = diag(c(100, 100, 10, 10))  # High uncertainty
        )
      }
    } else {
      # Initialize and update Kalman filter with current segment

      # Extract time differences for state transition
      times <- current_segment$time
      dt_vec <- diff(times)

      # Extract positions
      positions <- current_segment |>
        dplyr::select(x, y) |>
        as.matrix()

      # Initial state estimate [x, y, vx, vy]
      # Calculate velocities from first two points
      if (nrow(current_segment) >= 2) {
        init_vx <- (current_segment$x[2] - current_segment$x[1]) / dt_vec[1]
        init_vy <- (current_segment$y[2] - current_segment$y[1]) / dt_vec[1]
      } else {
        init_vx <- 0
        init_vy <- 0
      }

      a0 <- c(current_segment$x[1], current_segment$y[1], init_vx, init_vy)

      # Initial state covariance (uncertainty)
      P0 <- diag(c(1, 1, 1, 1))

      # Process noise covariance - higher values for velocity components
      Q <- diag(c(process_noise, process_noise, process_noise * 10, process_noise * 10))

      # Measurement noise covariance
      H <- matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0), 2, 4, byrow = TRUE)  # Measurement matrix

      R <- diag(rep(measurement_noise, 2))  # Measurement noise

      # Run the Kalman filter recursively through all points
      kf_result <- list(
        xf = a0,  # Filtered state
        Pf = P0   # Filtered state covariance
      )

      for (i in 2:nrow(current_segment)) {
        dt <- times[i] - times[i-1]

        # State transition matrix
        F <- matrix(c(
          1, 0, dt, 0,
          0, 1, 0, dt,
          0, 0, 1, 0,
          0, 0, 0, 1
        ), 4, 4, byrow = TRUE)

        # Prediction step
        x_pred <- F %*% kf_result$xf
        P_pred <- F %*% kf_result$Pf %*% t(F) + Q

        # Update step
        y <- matrix(positions[i, ], ncol = 1)  # Current measurement
        y_pred <- H %*% x_pred                 # Predicted measurement

        S <- H %*% P_pred %*% t(H) + R
        K <- P_pred %*% t(H) %*% solve(S)      # Kalman gain

        # Updated state and covariance
        kf_result$xf <- x_pred + K %*% (y - y_pred)
        kf_result$Pf <- (diag(4) - K %*% H) %*% P_pred
      }

      # Final state and covariance after filtering
      final_state <- kf_result$xf
      final_covariance <- kf_result$Pf

      # Create projection function using Kalman predictions
      project_position <- function(future_time) {
        dt <- future_time - current_end_time

        # State transition for projection
        F_proj <- matrix(c(
          1, 0, dt, 0,
          0, 1, 0, dt,
          0, 0, 1, 0,
          0, 0, 0, 1
        ), 4, 4, byrow = TRUE)

        # Project state forward
        pred_state <- F_proj %*% final_state
        pred_cov <- F_proj %*% final_covariance %*% t(F_proj) + Q * (dt/0.05)  # Scale process noise with time

        list(
          mean = pred_state,
          cov = pred_cov
        )
      }
    }

    # Find potential next segments
    candidates <- segments_summary |>
      dplyr::filter(segment_id > current_id) |>
      dplyr::slice_head(n = projection_dist)

    if (nrow(candidates) == 0) break

    # Calculate Mahalanobis distances to projected positions
    candidates <- candidates |>
      dplyr::mutate(
        delta_t = start_time - current_end_time
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        # Get Kalman prediction
        prediction = list(project_position(start_time)),
        # Extract components
        pred_x = prediction$mean[1],
        pred_y = prediction$mean[2],
        pred_vx = prediction$mean[3],
        pred_vy = prediction$mean[4],
        # Get position covariance submatrix (2x2 top-left of full covariance)
        pos_cov = list(prediction$cov[1:2, 1:2]),
        # Calculate Mahalanobis distance for position only
        diff_vec = list(c(start_x - pred_x, start_y - pred_y)),
        mahalanobis_dist = sqrt(t(diff_vec) %*% solve(pos_cov) %*% diff_vec),
        # Add time penalty
        distance = mahalanobis_dist + time_bias * delta_t,
        # Calculate probability (proportional to exp(-0.5 * mahalanobis_dist^2))
        probability = exp(-0.5 * mahalanobis_dist^2)
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-prediction, -pos_cov, -diff_vec) |>
      dplyr::arrange(distance, segment_id)

    if (nrow(candidates) > 0) {
      current_id <- candidates$segment_id[1]
      filtered_ids <- c(filtered_ids, current_id)
    } else break
  }

  filtered_ids
}

#' Linear Regression Path Connection Method
#'
#' Uses linear regression to predict and connect path segments with Euclidean distance ranking.
#' This method is computationally efficient but less sophisticated than Kalman filtering.
#' The function uses at most 40 points for regression to maintain efficiency.
#'
#' @param segments_summary Summary of path segments with start/end positions and times
#' @param processed_df The processed tracking data with segment IDs
#' @param time_bias Weight applied to temporal distance in path connection
#' @param projection_dist Maximum number of future segments to consider in the projection
#' @param min_points Minimum number of points to use for regression (defaults to 20)
#'
#' @return Vector of segment IDs that form the connected path
connect_linear_regression <- function(segments_summary, processed_df,
                                      time_bias = 10, projection_dist = 10,
                                      min_points = 20) {
  filtered_ids <- c(1L)
  current_id <- 1L

  # Keep track of end time for delta_t calculations
  current_end_time <- NULL

  # Set maximum number of points to use for regression
  max_points <- 40

  repeat {
    # Get current segment end position
    current_end <- segments_summary |>
      dplyr::filter(segment_id == current_id) |>
      dplyr::select(end_x, end_y, end_time)

    current_end_time <- current_end$end_time

    # Get current segment data
    current_segment <- processed_df |>
      dplyr::filter(as.integer(segment_id) == current_id) |>
      dplyr::arrange(time)

    # Check if we need to incorporate points from previous segments
    current_segment_points <- nrow(current_segment)
    if (current_segment_points < min_points && length(filtered_ids) > 1) {
      # Calculate how many additional points we need
      points_needed <- min(min_points - current_segment_points, max_points - current_segment_points)

      if (points_needed > 0) {
        # Get points from previous segments in the filtered path
        previous_segments <- processed_df |>
          dplyr::filter(as.integer(segment_id) %in% filtered_ids[-length(filtered_ids)]) |>
          dplyr::arrange(dplyr::desc(time)) |>
          dplyr::slice_head(n = points_needed)

        if (nrow(previous_segments) > 0) {
          # Combine previous points with current segment
          current_segment <- dplyr::bind_rows(previous_segments, current_segment) |>
            dplyr::arrange(time)
        }
      }
    }

    # If there are more than max_points, use only the most recent max_points
    if (nrow(current_segment) > max_points) {
      current_segment <- current_segment |>
        dplyr::slice_tail(n = max_points)
    }

    # If there's only one point in the segment, use simple linear projection
    if (nrow(current_segment) < 2) {
      med_vel <- stats::median(processed_df$velocity, na.rm = TRUE)

      # Create prediction function for simple linear extrapolation
      project_position <- function(future_time) {
        delta_t <- future_time - current_end_time
        pred_x <- current_end$end_x + med_vel * delta_t
        pred_y <- current_end$end_y

        c(pred_x, pred_y)
      }
    } else {
      # Use linear regression on the data (limited to max_points)
      use_points <- min(nrow(current_segment), max_points)
      recent_segment <- tail(current_segment, use_points)

      # Fit linear models for x and y based on time
      lm_x <- stats::lm(x ~ time, data = recent_segment)
      lm_y <- stats::lm(y ~ time, data = recent_segment)

      # Create prediction function using linear model
      project_position <- function(future_time) {
        pred_df <- data.frame(time = future_time)
        pred_x <- stats::predict(lm_x, newdata = pred_df)
        pred_y <- stats::predict(lm_y, newdata = pred_df)

        c(pred_x, pred_y)
      }
    }

    # Find potential next segments
    candidates <- segments_summary |>
      dplyr::filter(segment_id > current_id) |>
      dplyr::slice_head(n = projection_dist)

    if (nrow(candidates) == 0) break

    # Calculate Euclidean distances to projected positions
    candidates <- candidates |>
      dplyr::mutate(
        delta_t = start_time - current_end_time
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        # Get linear prediction
        prediction = list(project_position(start_time)),
        # Extract components
        pred_x = prediction[1],
        pred_y = prediction[2],
        # Calculate Euclidean distance
        euclidean_dist = sqrt((start_x - pred_x)^2 + (start_y - pred_y)^2),
        # Add time penalty
        distance = euclidean_dist + time_bias * delta_t
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-prediction) |>
      dplyr::arrange(distance, segment_id)

    if (nrow(candidates) > 0) {
      current_id <- candidates$segment_id[1]
      filtered_ids <- c(filtered_ids, current_id)
    } else break
  }

  filtered_ids
}
#' Spline Path Connection Method
#'
#' Uses spline interpolation to predict and connect path segments.
#' This method can capture non-linear trends better than simple linear regression.
#'
#' @param segments_summary Summary of path segments with start/end positions and times
#' @param processed_df The processed tracking data with segment IDs
#' @param time_bias Weight applied to temporal distance in path connection
#' @param projection_dist Maximum number of future segments to consider in the projection
#' @param min_points Minimum number of points to use for spline fitting (defaults to 25)
#'
#' @return Vector of segment IDs that form the connected path
#' Spline Path Connection Method using mgcv
#'
#' Uses GAM spline models from mgcv package to predict and connect path segments.
#' This method can capture non-linear trends better than simple linear regression.
#'
#' @param segments_summary Summary of path segments with start/end positions and times
#' @param processed_df The processed tracking data with segment IDs
#' @param time_bias Weight applied to temporal distance in path connection
#' @param projection_dist Maximum number of future segments to consider in the projection
#' @param min_points Minimum number of points to use for spline fitting (defaults to 25)
#'
#' @return Vector of segment IDs that form the connected path
connect_spline <- function(segments_summary, processed_df,
                           time_bias = 10, projection_dist = 10,
                           min_points = 25) {
  filtered_ids <- c(1L)
  current_id <- 1L

  # Keep track of end time for delta_t calculations
  current_end_time <- NULL

  repeat {
    # Get current segment end position
    current_end <- segments_summary |>
      dplyr::filter(segment_id == current_id) |>
      dplyr::select(end_x, end_y, end_time)

    current_end_time <- current_end$end_time

    # Get current segment data
    current_segment <- processed_df |>
      dplyr::filter(as.integer(segment_id) == current_id) |>
      dplyr::arrange(time)

    # Check if we need to incorporate points from previous segments
    current_segment_points <- nrow(current_segment)
    if (current_segment_points < min_points && length(filtered_ids) > 1) {
      # Calculate how many additional points we need
      points_needed <- min_points - current_segment_points

      # Get points from previous segments in the filtered path
      previous_segments <- processed_df |>
        dplyr::filter(as.integer(segment_id) %in% filtered_ids[-length(filtered_ids)]) |>
        dplyr::arrange(dplyr::desc(time)) |>
        dplyr::slice_head(n = points_needed)

      if (nrow(previous_segments) > 0) {
        # Combine previous points with current segment
        current_segment <- dplyr::bind_rows(previous_segments, current_segment) |>
          dplyr::arrange(time)
      }
    }

    # If there are not enough points, fall back to linear regression
    if (nrow(current_segment) < 5) {
      # Fall back to linear regression if not enough points for spline
      if (nrow(current_segment) < 2) {
        med_vel <- stats::median(processed_df$velocity, na.rm = TRUE)

        # Create prediction function for simple linear extrapolation
        project_position <- function(future_time) {
          delta_t <- future_time - current_end_time
          pred_x <- current_end$end_x + med_vel * delta_t
          pred_y <- current_end$end_y

          c(pred_x, pred_y)
        }
      } else {
        # Use linear regression instead
        lm_x <- stats::lm(x ~ time, data = current_segment)
        lm_y <- stats::lm(y ~ time, data = current_segment)

        # Create prediction function using linear model
        project_position <- function(future_time) {
          pred_df <- data.frame(time = future_time)
          pred_x <- stats::predict(lm_x, newdata = pred_df)
          pred_y <- stats::predict(lm_y, newdata = pred_df)

          c(pred_x, pred_y)
        }
      }
    } else {
      # Use mgcv GAM models for spline fitting
      use_points <- min(nrow(current_segment), min_points)
      recent_segment <- tail(current_segment, use_points)

      # Calculate number of unique time points for setting k
      unique_times <- length(unique(recent_segment$time))

      # Fit GAM models with optimal smoothing for x and y based on time
      gam_x <- mgcv::gam(x ~ s(time, k = min(10, unique_times)), data = recent_segment)
      gam_y <- mgcv::gam(y ~ s(time, k = min(10, unique_times)), data = recent_segment)

      # Create prediction function using GAM models
      project_position <- function(future_time) {
        # Ensure future_time is within a reasonable extrapolation range
        max_extrap_time <- max(recent_segment$time) + 0.25 * (max(recent_segment$time) - min(recent_segment$time))

        if (future_time > max_extrap_time) {
          # Linear extrapolation beyond spline range
          last_few <- tail(recent_segment, 5)
          lm_x_extrap <- stats::lm(x ~ time, data = last_few)
          lm_y_extrap <- stats::lm(y ~ time, data = last_few)

          pred_df <- data.frame(time = future_time)
          pred_x <- stats::predict(lm_x_extrap, newdata = pred_df)
          pred_y <- stats::predict(lm_y_extrap, newdata = pred_df)
        } else {
          # Use GAM prediction
          pred_df <- data.frame(time = future_time)
          pred_x <- stats::predict(gam_x, newdata = pred_df)
          pred_y <- stats::predict(gam_y, newdata = pred_df)
        }

        c(pred_x, pred_y)
      }
    }

    # Find potential next segments
    candidates <- segments_summary |>
      dplyr::filter(segment_id > current_id) |>
      dplyr::slice_head(n = projection_dist)

    if (nrow(candidates) == 0) break

    # Calculate Euclidean distances to projected positions
    candidates <- candidates |>
      dplyr::mutate(
        delta_t = start_time - current_end_time
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        # Get spline prediction
        prediction = list(project_position(start_time)),
        # Extract components
        pred_x = prediction[1],
        pred_y = prediction[2],
        # Calculate Euclidean distance
        euclidean_dist = sqrt((start_x - pred_x)^2 + (start_y - pred_y)^2),
        # Add time penalty
        distance = euclidean_dist + time_bias * delta_t
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-prediction) |>
      dplyr::arrange(distance, segment_id)

    if (nrow(candidates) > 0) {
      current_id <- candidates$segment_id[1]
      filtered_ids <- c(filtered_ids, current_id)
    } else break
  }

  filtered_ids
}
#' Find Path using selected path connection method
#'
#' This function identifies a continuous path by connecting trajectory segments using
#' the specified connection method.
#'
#' @param df Data frame containing tracking data with time, x, and y coordinates
#' @param method Method for jump detection: either "EM" (Expectation-Maximization) or "MAD" (Median Absolute Deviation)
#' @param connection_method Method for connecting path segments: "kalman" (Kalman filtering),
#'        "linear" (Linear regression), or "spline" (Spline interpolation)
#' @param time_bias Weight applied to temporal distance in path connection (higher values
#'        penalize larger time gaps more heavily)
#' @param report_jumps Logical; if TRUE, includes jump detection results in output
#' @param include_segment_id Logical; if TRUE, includes segment identifiers in output
#' @param prob_threshold Probability threshold for classifying jumps when using EM method
#' @param projection_dist Maximum number of future segments to consider in projection
#' @param process_noise Process noise parameter for Kalman filter (controls model uncertainty)
#' @param measurement_noise Measurement noise parameter for Kalman filter (controls observation uncertainty)
#' @param min_points Minimum number of points to use for the selected connection method
#'
#' @return Data frame with identified path
#' @export
find_path <- function(df, method = c("EM", "MAD"),
                      connection_method = c("kalman", "linear", "spline"),
                      time_bias = 10, report_jumps = TRUE,
                      include_segment_id = TRUE, prob_threshold = 0.95,
                      projection_dist = 10,
                      process_noise = 0.1, measurement_noise = 1.0,
                      min_points = 100) {
  method <- match.arg(method)
  connection_method <- match.arg(connection_method)
  ids <- get_individual_ids(df)
  if (length(ids) == 0) stop("No valid individual columns found")

  processed <- lapply(ids, function(id) {
    x_col <- paste0("x", id)
    y_col <- paste0("y", id)

    ind_df <- df |>
      dplyr::select(time, x = !!x_col, y = !!y_col) |>
      dplyr::arrange(time)

    original_df <- ind_df |>
      dplyr::mutate(
        original_x = x,
        original_y = y
      )

    processed_df <- original_df |>
      dplyr::filter(!is.na(x) & !is.na(y) & !is.na(time))

    if (nrow(processed_df) == 0) return(ind_df)

    jump_df <- if (method == "MAD") {
      detect_jumps_MAD(processed_df)
    } else {
      detect_jumps_EM(processed_df, prob_threshold = prob_threshold)
    }

    processed_df <- jump_df |>
      dplyr::mutate(
        segment_id = factor(cumsum(is_jump | dplyr::row_number() == 1))
      ) |>
      dplyr::group_by(segment_id) |>
      dplyr::mutate(
        segment_duration = max(time) - min(time),
        segment_length = dplyr::n()
      ) |>
      dplyr::ungroup()

    segments_summary <- processed_df |>
      dplyr::group_by(segment_id) |>
      dplyr::summarize(
        start_x = dplyr::first(x),
        start_y = dplyr::first(y),
        start_time = dplyr::first(time),
        end_x = dplyr::last(x),
        end_y = dplyr::last(y),
        end_time = dplyr::last(time),
        .groups = 'drop'
      ) |>
      dplyr::mutate(segment_id = as.integer(as.character(segment_id))) |>
      dplyr::arrange(segment_id)

    # Important fix for the type incompatibility issue:
    processed_df <- processed_df |>
      dplyr::mutate(segment_id = as.integer(as.character(segment_id)))

    # Use the selected connection method
    filtered_ids <- switch(connection_method,
                           kalman = connect_kalman_projection(segments_summary, processed_df,
                                                              time_bias, projection_dist,
                                                              process_noise, measurement_noise,
                                                              min_points = min_points),
                           linear = connect_linear_regression(segments_summary, processed_df,
                                                              time_bias, projection_dist,
                                                              min_points = min(min_points, 20)),
                           spline = connect_spline(segments_summary, processed_df,
                                                   time_bias, projection_dist,
                                                   min_points = min(min_points, 25))
    )

    cols <- c("time", "x", "y")
    if (report_jumps) cols <- c(cols, "is_jump")
    if (include_segment_id) cols <- c(cols, "segment_id")

    processed_selected <- processed_df |>
      dplyr::filter(segment_id %in% filtered_ids) |>
      dplyr::select(dplyr::all_of(cols))

    processed_joined <- processed_selected |>
      dplyr::right_join(original_df, by = "time") |>
      dplyr::rename(
        !!paste0("original_x", id) := original_x,
        !!paste0("original_y", id) := original_y,
        !!x_col := x.x,
        !!y_col := y.x
      )

    if (report_jumps && "is_jump" %in% names(processed_joined)) {
      processed_joined <- processed_joined |>
        dplyr::rename(!!paste0("is_jump", id) := is_jump)
    }

    if (include_segment_id && "segment_id" %in% names(processed_joined)) {
      processed_joined <- processed_joined |>
        dplyr::rename(!!paste0("segment_id", id) := segment_id)
    }

    columns_to_select <- c("time",
                           paste0("original_x", id),
                           paste0("original_y", id),
                           x_col,
                           y_col)
    if (report_jumps) columns_to_select <- c(columns_to_select, paste0("is_jump", id))
    if (include_segment_id) columns_to_select <- c(columns_to_select, paste0("segment_id", id))

    processed_ind <- processed_joined |>
      dplyr::select(dplyr::all_of(columns_to_select)) |>
      dplyr::arrange(time)

    processed_ind
  })

  final_df <- purrr::reduce(processed, function(a, b) {
    dplyr::full_join(a, b, by = "time")
  }, .init = df |> dplyr::select(time) |> dplyr::distinct())

  other_cols <- setdiff(names(df), c(unlist(lapply(ids, function(id) {
    c(paste0("x", id), paste0("y", id))
  })), "time"))

  if (length(other_cols) > 0) {
    final_df <- dplyr::left_join(final_df, df |> dplyr::select(time, dplyr::all_of(other_cols)), by = "time")
  }

  final_df
}

#' Complete Processing Pipeline with Selected Connection Method
#'
#' This function combines path finding and smoothing into a single pipeline.
#' It first identifies a continuous path using jump detection and the selected connection method,
#' then applies Savitzky-Golay smoothing to the resulting trajectory.
#'
#' @param df Data frame containing tracking data with time, x, and y coordinates
#' @param method Method for jump detection: either "EM" (Expectation-Maximization) or "MAD" (Median Absolute Deviation)
#' @param connection_method Method for connecting path segments: "kalman" (Kalman filtering),
#'        "linear" (Linear regression), or "spline" (Spline interpolation)
#' @param p Polynomial order for Savitzky-Golay smoothing (default is 3)
#' @param n Filter length for Savitzky-Golay smoothing (must be odd, default is 13)
#' @param time_bias Weight applied to temporal distance in path connection
#' @param process_noise Process noise parameter for Kalman filter (controls model uncertainty)
#' @param measurement_noise Measurement noise parameter for Kalman filter (controls observation uncertainty)
#' @param min_points Minimum number of points to use for the selected connection method
#'
#' @return Data frame with identified and smoothed path
#' @export
find_smooth_path <- function(df, method = c("EM", "MAD"),
                             connection_method = c("kalman", "linear", "spline"),
                             p = 3, n = 13, time_bias = 10,
                             process_noise = 0.1, measurement_noise = 1.0,
                             min_points = 100) {
  method <- match.arg(method)
  connection_method <- match.arg(connection_method)
  df |>
    find_path(method = method, connection_method = connection_method,
              time_bias = time_bias, process_noise = process_noise,
              measurement_noise = measurement_noise, min_points = min_points) |>
    smooth_path(p = p, n = n)
}

#' Detect Jumps Using Expectation-Maximization (GMM)
#'
#' This function identifies jumps or discontinuities in tracking data by fitting
#' a Gaussian mixture model (GMM) to velocity distributions using Expectation-Maximization.
#'
#' @param df Data frame containing tracking data with time, x, and y coordinates
#' @param prob_threshold Probability threshold above which a point is classified as a jump
#' @param min_velocity_points Minimum number of velocity points required to use GMM method
#' @param mad_multiplier Multiplier for Median Absolute Deviation (used as fallback)
#'
#' @return Data frame with additional columns for jump probability and classification
#' @export
detect_jumps_EM <- function(df, prob_threshold,
                            min_velocity_points = 50,
                            mad_multiplier = 8) {
  df <- df |>
    dplyr::mutate(
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),
      velocity = sqrt(delta_x^2 + delta_y^2) / delta_t,
      velocity = ifelse(is.infinite(velocity), NA, velocity)
    )

  velocities <- df$velocity[!is.na(df$velocity)]

  df$jump_prob <- NA_real_
  df$is_jump <- FALSE

  if (length(velocities) >= min_velocity_points) {
    gmm_fit <- tryCatch({
      set.seed(123)
      mixtools::normalmixEM(velocities, k = 2, maxit = 200)
    }, error = function(e) NULL)

    if (!is.null(gmm_fit)) {
      normal_comp <- which.min(gmm_fit$mu)
      valid_idx <- which(!is.na(df$velocity))
      df$jump_prob[valid_idx] <- 1 - gmm_fit$posterior[, normal_comp]
      df$is_jump <- df$jump_prob > prob_threshold
      return(df)
    }
  }

  mad_df <- detect_jumps_MAD(df, mad_multiplier = mad_multiplier)
  if (!"jump_prob" %in% names(mad_df)) {
    mad_df$jump_prob <- NA_real_
  }
  mad_df
}

#' Detect Jumps Using Median Absolute Deviation (MAD)
#'
#' This function identifies jumps or discontinuities in tracking data using
#' the Median Absolute Deviation method, which is robust to outliers.
#'
#' @param df Data frame containing tracking data with time, x, and y coordinates
#' @param mad_multiplier Multiplier for Median Absolute Deviation to set the threshold
#'
#' @return Data frame with additional columns for velocity metrics and jump classification
#' @export
detect_jumps_MAD <- function(df) {
  df |>
    dplyr::filter(!is.na(x) & !is.na(y) & !is.na(time)) |>
    arrange(time) |>
    dplyr::mutate(
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),
      distance = sqrt(delta_x^2 + delta_y^2),
      velocity = distance / delta_t,
      med_dist = stats::median(distance, na.rm = TRUE),
      mad_dist = stats::mad(distance, na.rm = TRUE),
      dist_threshold = med_dist + 70 * mad_dist,
      is_jump = distance > dist_threshold,
      jump_prob = ifelse(is_jump, 1.0, 0.0)
    )
}

#' Smooth Individual Trajectory
#'
#' Applies Savitzky-Golay filtering to smooth individual tracking trajectories.
#'
#' @param ind_df Data frame containing time, x, and y coordinates for an individual
#' @param p Polynomial order for Savitzky-Golay smoothing (default is 3)
#' @param n Filter length for Savitzky-Golay smoothing (must be odd, default is 13)
#' @param smooth_first Logical indicating whether to smooth before NA filling (default is TRUE)
#'                    If FALSE, NA values are filled first, then smoothing is applied
#'
#' @return Data frame with additional columns for smoothed x and y coordinates
smooth_individual <- function(ind_df, p = 3, n = 13, smooth_first = FALSE) {
  if (smooth_first) {
    # Original approach: smooth first, then fill NAs
    smoothed <- ind_df |>
      dplyr::filter(!is.na(x) & !is.na(y)) |>
      dplyr::arrange(time) |>
      dplyr::mutate(
        x_smooth = signal::sgolayfilt(x, p = p, n = n),
        y_smooth = signal::sgolayfilt(y, p = p, n = n)
      ) |>
      dplyr::select(time, x_smooth, y_smooth)

    ind_df |>
      dplyr::left_join(smoothed, by = "time") |>
      dplyr::mutate(
        x_smooth = zoo::na.approx(x_smooth, na.rm = FALSE),
        y_smooth = zoo::na.approx(y_smooth, na.rm = FALSE)
      )
  } else {
    # Alternative approach: fill NAs first, then smooth
    filled <- ind_df |>
      dplyr::arrange(time) |>
      dplyr::mutate(
        x_filled = zoo::na.approx(x, na.rm = FALSE),
        y_filled = zoo::na.approx(y, na.rm = FALSE)
      )

    filled |>
      dplyr::mutate(
        x_smooth = signal::sgolayfilt(x_filled, p = p, n = n),
        y_smooth = signal::sgolayfilt(y_filled, p = p, n = n)
      )
  }
}


#' Smooth path after finding it
#'
#' This function applies Savitzky-Golay filtering to smooth all identified paths.
#'
#' @param df Data frame containing tracking data with identified paths
#' @param p Polynomial order for Savitzky-Golay smoothing (default is 3)
#' @param n Filter length for Savitzky-Golay smoothing (must be odd, default is 13)
#'
#' @return Data frame with additional columns for smoothed x and y coordinates
#' @export
smooth_path <- function(df, p = 3, n = 13) {
  ids <- get_individual_ids(df)
  if (length(ids) == 0) return(df)

  smoothed <- lapply(ids, function(id) {
    x_col <- paste0("x", id)
    y_col <- paste0("y", id)

    if (!(x_col %in% names(df) && y_col %in% names(df))) return(NULL)

    ind_df <- df |>
      dplyr::select(time, x = !!x_col, y = !!y_col)

    smoothed_ind <- smooth_individual(ind_df, p = p, n = n)

    smoothed_ind |>
      dplyr::rename(
        !!paste0("x_smooth", id) := x_smooth,
        !!paste0("y_smooth", id) := y_smooth
      ) |>
      dplyr::select(-x, -y)
  })

  smoothed <- smoothed[!sapply(smoothed, is.null)]

  if (length(smoothed) == 0) return(df)

  final_smoothed <- purrr::reduce(smoothed, function(a, b) {
    dplyr::full_join(a, b, by = "time")
  }, .init = df)

  final_smoothed
}

#' Helper function to detect individual IDs from column names
#'
#' Extracts individual identifiers from column names that follow the pattern 'x1', 'y1', etc.
#'
#' @param df Data frame containing tracking data
#' @return Character vector of detected individual IDs
get_individual_ids <- function(df) {
  x_cols <- grep("^x\\d+$", names(df), value = TRUE)
  x_ids <- unique(sub("^x", "", x_cols))

  # Validate corresponding y columns exist
  valid_ids <- x_ids[paste0("y", x_ids) %in% names(df)]

  if (length(valid_ids) > 0) {
    num_ids <- suppressWarnings(as.numeric(valid_ids))
    if (!any(is.na(num_ids))) {
      return(as.character(sort(num_ids)))
    }
    return(sort(valid_ids))
  }
  character(0)
}

#' Convert seconds to ASS time format
#'
#' @param t Numeric value representing time in seconds
#' @return Character string in ASS time format (H:MM:SS.FF)
#' @export
seconds_to_ass_time <- function(t) {
  hours <- floor(t / 3600)
  remainder <- t %% 3600
  minutes <- floor(remainder / 60)
  seconds <- floor(remainder %% 60)
  centiseconds <- floor((t %% 1) * 100)
  sprintf("%d:%02d:%02d.%02d", hours, minutes, seconds, centiseconds)
}

#' Get video resolution using ffprobe
#'
#' @param video_path Character string of the video file path
#' @return List containing width and height of the video
#' @export
get_video_resolution <- function(video_path) {
  cmd <- paste(
    'ffprobe',
    '-v error',
    '-select_streams v:0',
    '-show_entries stream=width,height',
    '-of csv=p=0',
    shQuote(video_path)
  )

  output <- system(cmd, intern = TRUE)

  if (!grepl(",", output)) {
    stop(paste("Unexpected ffprobe output format:", output))
  }

  dimensions <- strsplit(output, ",")[[1]]
  width <- as.integer(dimensions[1])
  height <- as.integer(dimensions[2])

  return(list(width = width, height = height))
}

#' Update a position trail with a new position
#'
#' @param trail List of position coordinates
#' @param new_position List containing x and y coordinates
#' @param max_length Integer specifying maximum length of the trail
#' @return Updated trail list
#' @keywords internal
update_trail <- function(trail, new_position, max_length) {
  trail <- c(trail, list(new_position))
  if (length(trail) > max_length) {
    trail <- trail[(length(trail) - max_length + 1):length(trail)]
  }
  return(trail)
}

#' Generate ASS subtitle file content for fish tracking
#'
#' @param df Data frame containing fish tracking data
#' @param width Video width in pixels
#' @param height Video height in pixels
#' @param frame_rate Video frame rate
#' @param trail_length Number of previous positions to show
#' @param start_size Starting size for most recent position
#' @param size_step Size reduction per previous position
#' @param fish_colors Vector of color codes in RRGGBB format
#' @return Character vector containing ASS subtitle content
#' @keywords internal
generate_ass_content <- function(df, width, height, frame_rate,
                                 trail_length = 10, start_size = 30,
                                 size_step = 3, fish_colors = NULL) {
  # Default colors if not provided
  if (is.null(fish_colors)) {
    fish_colors <- c(
      'FF0000', '00FF00', '0000FF', 'FFFF00', 'FF00FF', '00FFFF',
      'FF4500', '8A2BE2', '32CD32', 'FF1493', 'FFD700', '1E90FF', 'FF69B4'
    )
  }

  # Initialize ASS content
  ass_content <- c(
    "[Script Info]",
    "ScriptType: v4.00+",
    paste0("PlayResX: ", width),
    paste0("PlayResY: ", height),
    "",
    "[V4+ Styles]",
    "Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding",
    "Style: Default,Arial,20,&H00FFFFFF,&H000000FF,&H00000000,&H00000000,0,0,0,0,100,100,0,0,1,1,0,7,0,0,0,0",
    "",
    "[Events]",
    "Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text"
  )

  # Determine number of fish to track
  fish_count <- 0
  for (i in 1:50) {  # Check up to 50 fish (adjust as needed)
    x_col <- paste0("x_smooth", i)
    if (x_col %in% names(df)) {
      fish_count <- i
    } else {
      break
    }
  }

  if (fish_count == 0) {
    stop("No fish tracking data columns found in the data frame")
  }

  # Initialize position history for each fish
  fish_trails <- vector("list", fish_count)
  for (i in 1:fish_count) {
    fish_trails[[i]] <- list()
  }

  # Process each frame
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    start_time <- row$time
    end_time <- start_time + (1 / frame_rate)
    start_ass <- seconds_to_ass_time(start_time)
    end_ass <- seconds_to_ass_time(end_time)

    for (fish_num in 1:fish_count) {
      x_col <- paste0("x_smooth", fish_num)
      y_col <- paste0("y_smooth", fish_num)

      # Skip if columns don't exist
      if (!(x_col %in% names(df) && y_col %in% names(df))) {
        next
      }

      x <- row[[x_col]]
      y <- row[[y_col]]

      # Skip if NA values
      if (is.na(x) || is.na(y)) {
        next
      }

      fish_index <- fish_num

      # Update fish trail
      fish_trails[[fish_index]] <- update_trail(
        fish_trails[[fish_index]],
        list(x = x, y = y),
        trail_length
      )
      trail_positions <- fish_trails[[fish_index]]

      # Create trail points with decreasing size
      for (j in 1:length(trail_positions)) {
        trail_pos <- trail_positions[[j]]
        trail_x <- trail_pos$x
        trail_y <- trail_pos$y

        # Calculate size based on position age
        n_positions <- length(trail_positions)
        steps_back <- (n_positions - j)
        font_size <- start_size - (steps_back * size_step)
        if (font_size <= 0) next

        # Convert color to ASS format (BBGGRR)
        color_hex <- fish_colors[(fish_index - 1) %% length(fish_colors) + 1]
        rr <- substr(color_hex, 1, 2)
        gg <- substr(color_hex, 3, 4)
        bb <- substr(color_hex, 5, 6)
        ass_color <- paste0("&H", bb, gg, rr)

        # Create subtitle event
        event_line <- paste0(
          "Dialogue: 0,", start_ass, ",", end_ass, ",Default,,0,0,0,,",
          "{\\an5\\pos(", trail_x, ",", trail_y, ")\\fs", font_size, "\\c", ass_color, "}●"
        )
        ass_content <- c(ass_content, event_line)
      }
    }
  }

  return(ass_content)
}

#' Create fish tracking overlay on video
#'
#' @param df Data frame containing fish tracking data with columns 'time', 'x_smooth1', 'y_smooth1', etc.
#' @param video_path Character string of the video file path
#' @param output_path Optional character string for the output video path. If NULL, creates path based on input video
#' @param trail_length Number of previous positions to show (default: 10)
#' @param start_size Starting size for most recent position (default: 30)
#' @param size_step Size reduction per previous position (default: 3)
#' @param fish_colors Optional vector of color codes in RRGGBB format
#' @param max_time Optional maximum time in seconds to process (default: NULL for entire video)
#' @param ffmpeg_options Optional list of additional ffmpeg options
#' @return Character string of the output video path
#' @export
create_fish_tracking_overlay <- function(df, video_path, output_path = NULL,
                                         trail_length = 10, start_size = 30,
                                         size_step = 3, fish_colors = NULL,
                                         max_time = NULL, ffmpeg_options = NULL) {
  # Check if data frame has required columns
  if (!("time" %in% names(df))) {
    stop("Data frame must contain a 'time' column")
  }

  if (!any(grepl("x_smooth", names(df)))) {
    stop("Data frame must contain at least one 'x_smooth' column")
  }

  # Verify video exists
  if (!file.exists(video_path)) {
    stop(paste("Video file", video_path, "not found"))
  }

  # Get video resolution
  tryCatch({
    res <- get_video_resolution(video_path)
    width <- res$width
    height <- res$height
  }, error = function(e) {
    stop(paste("Error getting video resolution:", e$message))
  })

  # Calculate frame rate from time differences
  if (nrow(df) < 2) {
    stop("Not enough data points to calculate frame rate")
  }

  frame_rate <- 1 / (df$time[2] - df$time[1])

  # Determine max_time if not provided
  if (is.null(max_time)) {
    max_time <- max(df$time) + (1 / frame_rate)
  }

  # Generate ASS content
  ass_content <- generate_ass_content(
    df = df,
    width = width,
    height = height,
    frame_rate = frame_rate,
    trail_length = trail_length,
    start_size = start_size,
    size_step = size_step,
    fish_colors = fish_colors
  )

  # Create temporary ASS file
  ass_filename <- tempfile(fileext = ".ass")
  writeLines(ass_content, ass_filename, useBytes = TRUE)
  on.exit(unlink(ass_filename), add = TRUE)

  # Determine output path if not provided
  if (is.null(output_path)) {
    video_path_parts <- strsplit(video_path, "/")[[1]]
    video_stem <- tools::file_path_sans_ext(video_path_parts[length(video_path_parts)])
    output_dir <- paste(video_path_parts[-length(video_path_parts)], collapse = "/")
    output_path <- file.path(output_dir, paste0(video_stem, "_tracked_spline.mp4"))
  }

  # Prepare ffmpeg command
  cmd_parts <- c(
    'ffmpeg', '-y',
    '-i', shQuote(video_path),
    '-vf', paste0("subtitles=", shQuote(ass_filename)),
    '-c:a', 'copy',
    '-c:v', 'libx264',
    '-preset', 'fast',
    '-crf', '23'
  )

  # Add max_time if specified
  if (!is.null(max_time)) {
    cmd_parts <- c(cmd_parts, '-t', as.character(max_time))
  }

  # Add any additional ffmpeg options
  if (!is.null(ffmpeg_options)) {
    for (opt in names(ffmpeg_options)) {
      cmd_parts <- c(cmd_parts, opt, ffmpeg_options[[opt]])
    }
  }

  # Add output path
  cmd_parts <- c(cmd_parts, shQuote(output_path))

  # Run ffmpeg command
  cmd <- paste(cmd_parts, collapse = " ")
  system_result <- system(cmd)

  if (system_result == 0) {
    message(paste("Successfully created overlay video:", output_path))
  } else {
    warning(paste("Error processing video with exit code:", system_result))
  }

  return(output_path)
}

#' Create fish tracking overlay from file
#'
#' @param tracking_data_path Character string of the path to the tracking data CSV file
#' @param video_path Character string of the video file path
#' @param output_path Optional character string for the output video path
#' @param trail_length Number of previous positions to show (default: 10)
#' @param start_size Starting size for most recent position (default: 30)
#' @param size_step Size reduction per previous position (default: 3)
#' @param fish_colors Optional vector of color codes in RRGGBB format
#' @param max_time Optional maximum time in seconds to process
#' @param ffmpeg_options Optional list of additional ffmpeg options
#' @return Character string of the output video path
#' @export
create_fish_tracking_overlay_from_file <- function(tracking_data_path, video_path,
                                                   output_path = NULL,
                                                   trail_length = 10, start_size = 30,
                                                   size_step = 3, fish_colors = NULL,
                                                   max_time = NULL, ffmpeg_options = NULL) {
  # Load tracking data
  df <- data.table::fread(tracking_data_path)

  # Call the main function
  create_fish_tracking_overlay(
    df = df,
    video_path = video_path,
    output_path = output_path,
    trail_length = trail_length,
    start_size = start_size,
    size_step = size_step,
    fish_colors = fish_colors,
    max_time = max_time,
    ffmpeg_options = ffmpeg_options
  )
}
