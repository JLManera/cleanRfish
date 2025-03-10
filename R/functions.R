# path_functions.R

#' Identify Connected Path Segments
#'
#' Processes trajectory data to identify continuous path segments by removing jumps caused by
#' tracking artifacts or measurement errors. Maintains original dataframe structure with NAs
#' where segments are removed to preserve temporal alignment.
#'
#' @param df Dataframe containing trajectory data with required columns:
#'   - x: X-coordinate values
#'   - y: Y-coordinate values
#'   - time: Timestamp values (numeric or POSIXct)
#' @return Dataframe matching input length with columns:
#'   - time: Original timestamps
#'   - original_x: Preserved original x-coordinates
#'   - original_y: Preserved original y-coordinates
#'   - x: Processed x-coordinates (NAs at discontinuity points)
#'   - y: Processed y-coordinates (NAs at discontinuity points)
#' @export
find_path <- function(df) {
  # Validate input structure
  if (!all(c("x", "y", "time") %in% names(df))) {
    stop("Input must contain 'x', 'y', and 'time' columns")
  }

  # Preserve original data and mark NAs
  original_df <- df |>
    dplyr::arrange(time) |>
    dplyr::mutate(
      original_x = x,
      original_y = y
    )

  # Process non-NA points for segment analysis
  processed_df <- original_df |>
    dplyr::filter(!is.na(x) & !is.na(y) & !is.na(time)) |>
    dplyr::mutate(
      # Calculate differentials between consecutive points
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),

      # Calculate instantaneous velocity (distance/time)
      velocity = sqrt(.data$delta_x^2 + .data$delta_y^2) / .data$delta_t,
      velocity = ifelse(is.infinite(.data$velocity), NA, .data$velocity),

      # Calculate robust velocity statistics
      med_velocity = stats::median(.data$velocity, na.rm = TRUE),
      mad_velocity = stats::mad(.data$velocity, na.rm = TRUE),

      # Set velocity threshold using Median Absolute Deviation (MAD)
      velocity_threshold = .data$med_velocity + 8 * .data$mad_velocity,

      # Identify jumps as points exceeding velocity threshold
      is_jump = .data$velocity > .data$velocity_threshold,

      # Create segment IDs based on jump locations
      segment_id = factor(cumsum(.data$is_jump | dplyr::row_number() == 1))
    ) |>
    dplyr::group_by(.data$segment_id) |>
    dplyr::mutate(
      # Calculate segment duration and length
      segment_duration = max(time) - min(time),
      segment_length = dplyr::n()
    ) |>
    dplyr::ungroup()

  # Analyze segment connections
  segments_summary <- processed_df |>
    dplyr::group_by(.data$segment_id) |>
    dplyr::summarize(
      start_x = dplyr::first(x),
      start_y = dplyr::first(y),
      start_time = dplyr::first(time),
      end_x = dplyr::last(x),
      end_y = dplyr::last(y),
      end_time = dplyr::last(time),
      .groups = 'drop'
    ) |>
    dplyr::mutate(segment_id = as.integer(as.character(.data$segment_id))) |>
    dplyr::arrange(.data$segment_id)

  # Segment connection logic using velocity projection
  filtered_ids <- c(1L)  # Initialize with first segment
  current_id <- 1L

  repeat {
    current_end <- segments_summary |>
      dplyr::filter(.data$segment_id == current_id) |>
      dplyr::select(end_x, end_y, end_time)

    # Analyze last 10 points of current segment for velocity estimation
    current_segment <- processed_df |>
      dplyr::filter(as.integer(.data$segment_id) == current_id) |>
      dplyr::arrange(dplyr::desc(time)) |>
      dplyr::slice_head(n = 10)

    # Calculate velocity parameters using linear regression
    if (nrow(current_segment) >= 2) {
      x_model <- stats::lm(x ~ time, data = current_segment)
      y_model <- stats::lm(y ~ time, data = current_segment)
      dx_dt <- stats::coef(x_model)[2]  # x-velocity
      dy_dt <- stats::coef(y_model)[2]  # y-velocity
    } else {
      dx_dt <- processed_df$med_velocity[1]
      dy_dt <- 0
    }

    # Find best candidate segments for connection
    candidates <- segments_summary |>
      dplyr::filter(.data$segment_id > current_id) |>
      dplyr::slice_head(n = 10)

    if (nrow(candidates) == 0) break

    candidates <- candidates |>
      dplyr::mutate(
        delta_t = .data$start_time - current_end$end_time,
        # Project expected position based on current velocity
        proj_x = current_end$end_x + dx_dt * .data$delta_t,
        proj_y = current_end$end_y + dy_dt * .data$delta_t,
        # Calculate composite distance metric
        distance = sqrt((.data$start_x - .data$proj_x)^2 +
                          (.data$start_y - .data$proj_y)^2 +
                          (10 * .data$delta_t)^2)  # Time penalty factor
      ) |>
      dplyr::arrange(.data$distance, .data$segment_id)

    if (nrow(candidates) > 0) {
      current_id <- candidates$segment_id[1]
      filtered_ids <- c(filtered_ids, current_id)
    } else break
  }

  # Merge results with original data structure
  processed_df |>
    dplyr::mutate(segment_num = as.integer(.data$segment_id)) |>
    dplyr::filter(.data$segment_num %in% filtered_ids) |>
    dplyr::select(time, x, y) |>
    dplyr::right_join(original_df, by = "time") |>
    dplyr::select(time, original_x, original_y, x = x.x, y = y.x) |>
    dplyr::arrange(time)
}

#' Smooth Path Coordinates with Savitzky-Golay Filter
#'
#' Applies Savitzky-Golay smoothing to path coordinates while handling missing values
#' through linear interpolation. Maintains original data structure and temporal alignment.
#'
#' @param path_df Output dataframe from find_path() containing processed coordinates
#' @param p Polynomial order for Savitzky-Golay filter (must be < n)
#' @param n Window size (number of points) for Savitzky-Golay filter (must be odd)
#' @return Smoothed dataframe matching original structure with columns:
#'   - time: Original timestamps
#'   - original_x: Preserved original x-coordinates
#'   - original_y: Preserved original y-coordinates
#'   - x: Smoothed x-coordinates
#'   - y: Smoothed y-coordinates
#' @export
smooth_path <- function(path_df, p = 3, n = 13) {
  # Create working copy and NA mask
  working_df <- path_df |>
    dplyr::arrange(time)

  original_times <- working_df |>
    mutate(place_holder = 1) |>
    select(time, place_holder)

  smoothed_df <- working_df |>
      dplyr::filter(!is.na(x) & !is.na(y)) |>
      dplyr::mutate(
        x = signal::sgolayfilt(x, p = p, n = n),
        y = signal::sgolayfilt(y, p = p, n = n)
      )

  smoothed_df <- smoothed_df |>
      dplyr::right_join(original_times, by = "time") |>
      dplyr::arrange(time)

  smoothed_df <- smoothed_df |>
      dplyr::arrange(time) |>
      dplyr::mutate(
        x = zoo::na.approx(x, na.rm = FALSE),
        y = zoo::na.approx(y, na.rm = FALSE)
      ) |>
      dplyr::select(time, x, y)


  return(smoothed_df)
}

#' Combined Path Processing Pipeline
#'
#' Wrapper function that sequentially applies path identification and smoothing.
#'
#' @param df Input dataframe with raw trajectory data
#' @param p Polynomial order for Savitzky-Golay filter (must be < n)
#' @param n Window size (number of points) for Savitzky-Golay filter (must be odd)
#' @return Fully processed dataframe with smoothed coordinates
#' @export
find_smooth_path <- function(df, p = 3, n = 13) {
  df |>
    find_path() |>
    smooth_path(p = p, n = n)
}

