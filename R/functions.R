# path_functions.R

#' Detect Jumps Using Expectation-Maximization (GMM)
#'
#' Identifies trajectory jumps by modeling velocities with a 2-component Gaussian Mixture Model
#' (normal movements vs jumps) using the EM algorithm. Automatically falls back to MAD method
#' if insufficient data for GMM.
#'
#' @param df A dataframe containing trajectory data with columns:
#'   \itemize{
#'     \item{x - Numeric vector of x-coordinates}
#'     \item{y - Numeric vector of y-coordinates}
#'     \item{time - Numeric or POSIXct vector of timestamps}
#'   }
#' @param prob_threshold Numeric [0-1]. Probability threshold for classifying jumps
#'        (default: 0.9). Values > threshold are considered jumps.
#' @param min_velocity_points Integer. Minimum number of valid velocity measurements
#'        required for GMM fitting (default: 50). If not met, falls back to MAD method.
#' @param mad_multiplier Numeric. MAD multiplier for fallback threshold calculation
#'        (default: 8). Only used if GMM fails.
#'
#' @return A dataframe with added columns:
#'   \itemize{
#'     \item{delta_x - X displacement between consecutive points}
#'     \item{delta_y - Y displacement between consecutive points}
#'     \item{delta_t - Time difference between consecutive points}
#'     \item{velocity - Calculated instantaneous velocity}
#'     \item{jump_prob - Probability of being a jump (NA if using MAD fallback)}
#'     \item{is_jump - Logical jump indicator}
#'   }
#'
#' @examples
#' \dontrun{
#' cleaned_data <- detect_jumps_EM(raw_data, prob_threshold = 0.95)
#' }
#'
#' @references
#' Benaglia T., Chauveau D., Hunter D. R., Young D. (2009). mixtools: An R Package for
#' Analyzing Finite Mixture Models. Journal of Statistical Software, 32(6), 1-29.
#'
#' @seealso \code{\link{detect_jumps_MAD}} for MAD-based jump detection
#' @export
detect_jumps_EM <- function(df, prob_threshold = 0.9,
                            min_velocity_points = 50,
                            mad_multiplier = 8) {
  # Calculate velocities
  df <- df |>
    dplyr::mutate(
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),
      velocity = sqrt(delta_x^2 + delta_y^2) / delta_t,
      velocity = ifelse(is.infinite(velocity), NA, velocity)
    )

  velocities <- df$velocity[!is.na(df$velocity)]

  # Initialize critical columns FIRST
  df$jump_prob <- NA_real_
  df$is_jump <- FALSE

  if (length(velocities) >= min_velocity_points) {
    gmm_fit <- tryCatch({
      set.seed(123)
      mixtools::normalmixEM(velocities, k = 2, maxit = 100)
    }, error = function(e) NULL)

    if (!is.null(gmm_fit)) {
      normal_comp <- which.min(gmm_fit$mu)
      valid_idx <- which(!is.na(df$velocity))
      df$jump_prob[valid_idx] <- 1 - gmm_fit$posterior[, normal_comp]
      df$is_jump <- df$jump_prob > prob_threshold
      return(df)
    }
  }

  # Fallback to MAD with jump_prob guarantee
  mad_df <- detect_jumps_MAD(df, mad_multiplier = mad_multiplier)
  if (!"jump_prob" %in% names(mad_df)) {
    mad_df$jump_prob <- NA_real_
  }
  mad_df
}

#' Detect Jumps Using Median Absolute Deviation (MAD)
#'
#' Identifies jumps using robust statistical thresholding based on median absolute deviation
#' of velocities. Suitable for datasets with limited observations or when parametric
#' assumptions are questionable.
#'
#' @param df A dataframe (see \code{detect_jumps_EM} for required columns)
#' @param mad_multiplier Numeric. Number of MADs above median to consider as jumps
#'        (default: 8). Higher values reduce false positives but increase false negatives.
#'
#' @return A dataframe with same added columns as \code{detect_jumps_EM} (excluding jump_prob)
#'
#' @examples
#' \dontrun{
#' cleaned_data <- detect_jumps_MAD(raw_data, mad_multiplier = 6)
#' }
#'
#' @references
#' Leys C., Ley C., Klein O., Bernard P., Licata L. (2013). Detecting outliers: Do not use
#' standard deviation around the mean, use absolute deviation around the median. Journal of
#' Experimental Social Psychology, 49(4), 764-766.
#'
#' @seealso \code{\link{detect_jumps_EM}} for probabilistic jump detection
#' @export
detect_jumps_MAD <- function(df, mad_multiplier = 8) {
  df |>
    dplyr::mutate(
      delta_x = x - dplyr::lag(x),
      delta_y = y - dplyr::lag(y),
      delta_t = time - dplyr::lag(time),
      velocity = sqrt(delta_x^2 + delta_y^2) / delta_t,
      velocity = ifelse(is.infinite(velocity), NA, velocity),
      med_velocity = stats::median(velocity, na.rm = TRUE),
      mad_velocity = stats::mad(velocity, na.rm = TRUE),
      velocity_threshold = med_velocity + mad_multiplier * mad_velocity,
      is_jump = velocity > velocity_threshold
    )
}

#' Identify Continuous Path Segments in Trajectory Data
#'
#' Processes trajectory data to identify biologically plausible path segments by removing discontinuities,
#' using either Expectation-Maximization (GMM) or Median Absolute Deviation (MAD) for jump detection.
#' Maintains temporal alignment by preserving original timestamps and replacing removed points with NAs.
#'
#' @param df A dataframe containing raw trajectory data with columns:
#'   \itemize{
#'     \item{x - Numeric vector of x-coordinates (NA allowed)}
#'     \item{y - Numeric vector of y-coordinates (NA allowed)}
#'     \item{time - Numeric/POSIXct timestamp vector (NA not allowed)}
#'   }
#' @param method Jump detection algorithm:
#'   \itemize{
#'     \item{"EM" - Gaussian Mixture Model (recommended for datasets >50 observations)}
#'     \item{"MAD" - Robust median-based method (faster, suitable for small datasets)}
#'   }
#'   Default: "EM"
#'
#' @return A dataframe with same rows as input, containing:
#'   \itemize{
#'     \item{time - Original timestamps (ordered)}
#'     \item{original_x - Preserved input x-coordinates}
#'     \item{original_y - Preserved input y-coordinates}
#'     \item{x - Processed x-coordinates (NA at discontinuities)}
#'     \item{y - Processed y-coordinates (NA at discontinuities)}
#'   }
#'
#' @examples
#' \dontrun{
#' # Simulate trajectory with jumps
#' set.seed(123)
#' test_data <- data.frame(
#'   time = 1:100,
#'   x = cumsum(c(0, rnorm(99)) + rep(c(0,10), c(90,10)),
#'   y = cumsum(c(0, rnorm(99)))
#'
#' # Process using EM method
#' cleaned_path <- find_path(test_data, method = "EM")
#'
#' # Process using MAD method
#' cleaned_path_mad <- find_path(test_data, method = "MAD")
#' }
#'
#' @note
#' Key processing steps:
#' \enumerate{
#'   \item Orders data by timestamp
#'   \item Detects jumps using specified method
#'   \item Creates path segments between jumps
#'   \item Connects segments using velocity-projection heuristics
#'   \item Maintains original temporal structure with NA gaps
#' }
#' Computational complexity is O(n) for MAD method, O(n + km) for EM method where
#' k is number of GMM components and m is EM iterations. For large datasets (>1e5 points),
#' consider downsampling first.
#'
#' @seealso
#' \code{\link{detect_jumps_EM}} for GMM jump detection implementation
#' \code{\link{detect_jumps_MAD}} for MAD thresholding implementation
#' \code{\link{smooth_path}} for subsequent path smoothing
#'
#' @importFrom dplyr arrange filter mutate group_by ungroup summarize first last n
#' @importFrom stats median lm coef
#' @export
find_path <- function(df, method = c("EM", "MAD")) {
  method <- match.arg(method)

  # Input validation and preprocessing
  if (!all(c("x", "y", "time") %in% names(df))) {
    stop("Input must contain 'x', 'y', and 'time' columns")
  }

  original_df <- df |>
    dplyr::arrange(time) |>
    dplyr::mutate(
      original_x = x,
      original_y = y
    )

  processed_df <- original_df |>
    dplyr::filter(!is.na(x) & !is.na(y) & !is.na(time))

  # Apply selected jump detection
  if (method == "EM") {
    jump_df <- detect_jumps_EM(processed_df)
  } else {
    jump_df <- detect_jumps_MAD(processed_df)
  }

  # Segment identification
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

  # Segment connection analysis
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

  # Segment connection logic
  filtered_ids <- c(1L)
  current_id <- 1L

  repeat {
    current_end <- segments_summary |>
      dplyr::filter(segment_id == current_id) |>
      dplyr::select(end_x, end_y, end_time)

    current_segment <- processed_df |>
      dplyr::filter(as.integer(segment_id) == current_id) |>
      dplyr::arrange(dplyr::desc(time)) |>
      dplyr::slice_head(n = 10)

    # Velocity estimation
    if (nrow(current_segment) >= 2) {
      x_model <- stats::lm(x ~ time, data = current_segment)
      y_model <- stats::lm(y ~ time, data = current_segment)
      dx_dt <- stats::coef(x_model)[2]
      dy_dt <- stats::coef(y_model)[2]
    } else {
      dx_dt <- stats::median(processed_df$velocity, na.rm = TRUE)
      dy_dt <- 0
    }

    # Candidate selection
    candidates <- segments_summary |>
      dplyr::filter(segment_id > current_id) |>
      dplyr::slice_head(n = 10)

    if (nrow(candidates) == 0) break

    candidates <- candidates |>
      dplyr::mutate(
        delta_t = start_time - current_end$end_time,
        proj_x = current_end$end_x + dx_dt * delta_t,
        proj_y = current_end$end_y + dy_dt * delta_t,
        distance = sqrt((start_x - proj_x)^2 +
                        (start_y - proj_y)^2 +
                        (10 * delta_t)^2)
      ) |>
      dplyr::arrange(distance, segment_id)

    if (nrow(candidates) > 0) {
      current_id <- candidates$segment_id[1]
      filtered_ids <- c(filtered_ids, current_id)
    } else break
  }

  # Final assembly
  processed_df |>
    dplyr::mutate(segment_num = as.integer(segment_id)) |>
    dplyr::filter(segment_num %in% filtered_ids) |>
    dplyr::select(time, x, y) |>
    dplyr::right_join(original_df, by = "time") |>
    dplyr::select(time, original_x, original_y, x = x.x, y = y.x) |>
    dplyr::arrange(time)
}

#' Savitzky-Golay Path Smoothing
#'
#' @param path_df Output from find_path
#' @param p Polynomial order
#' @param n Filter window size
#' @return Smoothed path dataframe
#' @export
smooth_path <- function(path_df, p = 3, n = 13) {
  original_times <- path_df |>
    dplyr::mutate(placeholder = 1) |>
    dplyr::select(time, placeholder)

  path_df |>
    dplyr::filter(!is.na(x) & !is.na(y)) |>
    dplyr::mutate(
      x = signal::sgolayfilt(x, p = p, n = n),
      y = signal::sgolayfilt(y, p = p, n = n)
    ) |>
    dplyr::right_join(original_times, by = "time") |>
    dplyr::arrange(time) |>
    dplyr::mutate(
      x = zoo::na.approx(x, na.rm = FALSE),
      y = zoo::na.approx(y, na.rm = FALSE)
    ) |>
    dplyr::select(-placeholder)
}

#' Complete Processing Pipeline
#'
#' @param df Raw input data
#' @param method Jump detection method
#' @param p,n Smoothing parameters
#' @return Fully processed path
#' @export
find_smooth_path <- function(df, method = c("EM", "MAD"), p = 3, n = 13) {
  method <- match.arg(method)
  df |>
    find_path(method = method) |>
    smooth_path(p = p, n = n)
}
