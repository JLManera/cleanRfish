# path_functions.R

#' Detect Jumps Using Expectation-Maximization (GMM)
#' @export
detect_jumps_EM <- function(df, prob_threshold = 0.99,
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

# Helper function to detect individual IDs from column names
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

#' Identify Continuous Path Segments (Multi-Individual Support)
#' @export
find_path <- function(df, method = c("EM", "MAD")) {
  method <- match.arg(method)
  ids <- get_individual_ids(df)
  if (length(ids) == 0) stop("No valid individual columns found")

  processed <- lapply(ids, function(id) {
    x_col <- paste0("x", id)
    y_col <- paste0("y", id)

    ind_df <- df |>
      dplyr::select(time, x = !!x_col, y = !!y_col) |>
      dplyr::arrange(time)

    # Original processing logic
    original_df <- ind_df |>
      dplyr::mutate(
        original_x = x,
        original_y = y
      )

    processed_df <- original_df |>
      dplyr::filter(!is.na(x) & !is.na(y) & !is.na(time))

    if (nrow(processed_df) == 0) return(ind_df)

    # Apply jump detection
    jump_df <- if (method == "MAD") {
      detect_jumps_MAD(processed_df)
    } else {
      detect_jumps_EM(processed_df)
    }

    # Segment processing
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

      if (nrow(current_segment) >= 2) {
        x_model <- stats::lm(x ~ time, data = current_segment)
        y_model <- stats::lm(y ~ time, data = current_segment)
        dx_dt <- stats::coef(x_model)[2]
        dy_dt <- stats::coef(y_model)[2]
      } else {
        dx_dt <- stats::median(processed_df$velocity, na.rm = TRUE)
        dy_dt <- 0
      }

      candidates <- segments_summary |>
        dplyr::filter(segment_id > current_id) |>
        dplyr::slice_head(n = 10)

      if (nrow(candidates) == 0) break

      candidates <- candidates |>
        dplyr::mutate(
          delta_t = start_time - current_end$end_time,
          proj_x = current_end$end_x + dx_dt * delta_t,
          proj_y = current_end$end_y + dy_dt * delta_t,
          distance = sqrt((start_x - proj_x)^2 + (start_y - proj_y)^2 + (10 * delta_t)^2)
        ) |>
        dplyr::arrange(distance, segment_id)

      if (nrow(candidates) > 0) {
        current_id <- candidates$segment_id[1]
        filtered_ids <- c(filtered_ids, current_id)
      } else break
    }

    processed_ind <- processed_df |>
      dplyr::mutate(segment_num = as.integer(segment_id)) |>
      dplyr::filter(segment_num %in% filtered_ids) |>
      dplyr::select(time, x, y) |>
      dplyr::right_join(original_df, by = "time") |>
      dplyr::select(time, original_x, original_y, x = x.x, y = y.x) |>
      dplyr::arrange(time) |>
      dplyr::rename(
        !!paste0("original_x", id) := original_x,
        !!paste0("original_y", id) := original_y,
        !!x_col := x,
        !!y_col := y
      )
  })

  # Combine all individuals
  final_df <- purrr::reduce(processed, function(a, b) {
    dplyr::full_join(a, b, by = "time")
  }, .init = df |> dplyr::select(time) |> dplyr::distinct())

  # Preserve non-coordinate columns
  other_cols <- setdiff(names(df), c(unlist(lapply(ids, function(id)
    c(paste0("x", id), paste0("y", id)))), "time"))
    if (length(other_cols) > 0) {
      final_df <- dplyr::left_join(final_df, df |> dplyr::select(time, !!other_cols),
                                   by = "time")
    }

    final_df
}

#' Savitzky-Golay Smoothing (Multi-Individual Support)
#' @export
smooth_path <- function(path_df, p = 3, n = 13) {
  ids <- get_individual_ids(path_df)

  smoothed <- lapply(ids, function(id) {
    x_col <- paste0("x", id)
    y_col <- paste0("y", id)

    ind_df <- path_df |>
      dplyr::select(time, x = !!x_col, y = !!y_col)

    smoothed <- ind_df |>
      dplyr::filter(!is.na(x) & !is.na(y)) |>
      dplyr::mutate(
        x = signal::sgolayfilt(x, p = p, n = n),
        y = signal::sgolayfilt(y, p = p, n = n)
      ) |>
      dplyr::right_join(ind_df |> dplyr::select(time), by = "time") |>
      dplyr::arrange(time) |>
      dplyr::mutate(
        x = zoo::na.approx(x, na.rm = FALSE),
        y = zoo::na.approx(y, na.rm = FALSE)
      ) |>
      dplyr::rename(
        !!x_col := x,
        !!y_col := y
      )
  })

  final_df <- purrr::reduce(smoothed, function(a, b) {
    dplyr::full_join(a, b, by = "time")
  }, .init = path_df |> dplyr::select(time) |> dplyr::distinct())

  # Preserve non-coordinate columns
  other_cols <- setdiff(names(path_df), c(unlist(lapply(ids, function(id)
    c(paste0("x", id), paste0("y", id)))), "time"))
    if (length(other_cols) > 0) {
      final_df <- dplyr::left_join(final_df, path_df |> dplyr::select(time, !!other_cols),
                                   by = "time")
    }

    final_df
}

#' Complete Processing Pipeline
#' @export
find_smooth_path <- function(df, method = c("EM", "MAD"), p = 3, n = 13) {
  method <- match.arg(method)
  df |>
    find_path(method = method) |>
    smooth_path(p = p, n = n)
}
