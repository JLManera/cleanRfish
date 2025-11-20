# Tests for cleanRfish V2

test_that("detect_jumps identifies discontinuities", {
  set.seed(1)
  df <- data.frame(
    time = seq(0, 10, by = 0.1),
    x = cumsum(rnorm(101, mean = 0, sd = 1)),
    y = cumsum(rnorm(101, mean = 0, sd = 1))
  )
  
  # Add a jump at row 50
  df$x[50] <- df$x[49] + 50
  df$y[50] <- df$y[49] + 50
  
  result <- detect_jumps(df, z_thresh = 3)
  
  expect_true("jump_flag" %in% names(result))
  expect_true("speed" %in% names(result))
  expect_true(any(result$jump_flag, na.rm = TRUE))
})

test_that("assign_segments creates sequential segment IDs", {
  set.seed(2)
  df <- data.frame(
    time = seq(0, 5, by = 0.1),
    x = cumsum(rnorm(51)),
    y = cumsum(rnorm(51))
  )
  
  df <- detect_jumps(df, z_thresh = 2)
  result <- assign_segments(df)
  
  expect_true("segment" %in% names(result))
  expect_true(is.integer(result$segment))
  expect_true(all(diff(result$segment) >= 0))
})

test_that("process_multi_individual handles multi-individual data", {
  skip_if_not_installed("mgcv")
  
  set.seed(3)
  n_points <- 200
  df_multi <- data.frame(
    time = seq(0, 10, length.out = n_points),
    x1 = cumsum(rnorm(n_points, mean = 0, sd = 1)),
    y1 = cumsum(rnorm(n_points, mean = 0, sd = 1)),
    x2 = cumsum(rnorm(n_points, mean = 0, sd = 1)),
    y2 = cumsum(rnorm(n_points, mean = 0, sd = 1))
  )
  
  # Add some jumps
  df_multi$x1[100] <- df_multi$x1[99] + 30
  df_multi$x2[150] <- df_multi$x2[149] + 30
  
  result <- process_multi_individual(
    df_raw = df_multi,
    window = 5,
    min_movement = 10,
    diagnostic_plots = FALSE
  )
  
  expect_true("tracks" %in% names(result))
  expect_true("individual_number" %in% names(result$tracks))
  expect_true(all(c("x_original", "y_original", "x_reconstructed", "y_reconstructed") %in% 
                    names(result$tracks)))
  expect_equal(length(unique(result$tracks$individual_number)), 2)
})

test_that("compute_uncertainty_model returns a function", {
  set.seed(4)
  df <- data.frame(
    time = seq(0, 10, by = 0.1),
    x = cumsum(rnorm(101)),
    y = cumsum(rnorm(101)),
    speed = abs(rnorm(101, mean = 5, sd = 1)),
    direction = runif(101, 0, 2 * pi)
  )
  
  df <- detect_jumps(df)
  df <- assign_segments(df)
  
  uncertainty_fn <- compute_uncertainty_model(df, window = 5)
  
  expect_true(is.function(uncertainty_fn))
  expect_true(is.numeric(uncertainty_fn(c(0.5, 1, 2))))
  expect_equal(length(uncertainty_fn(c(0.5, 1, 2))), 3)
})

test_that("smooth_track adds smoothed columns", {
  set.seed(5)
  df <- data.frame(
    time = seq(0, 10, by = 0.1),
    x_reconstructed = cumsum(rnorm(101)),
    y_reconstructed = cumsum(rnorm(101))
  )
  
  result <- smooth_track(df, method = "savitzky_golay", n = 11, p = 3)
  
  expect_true("x_smooth" %in% names(result))
  expect_true("y_smooth" %in% names(result))
})

test_that("get_speed_threshold returns positive value", {
  set.seed(6)
  df <- data.frame(
    time = seq(0, 10, by = 0.1),
    x = cumsum(rnorm(101)),
    y = cumsum(rnorm(101))
  )
  
  df <- detect_jumps(df)
  threshold <- get_speed_threshold(df, prob = 0.999)
  
  expect_true(is.numeric(threshold))
  expect_true(threshold > 0)
  expect_length(threshold, 1)
})

test_that("identify_tracking_gaps finds large gaps", {
  set.seed(7)
  # Create continuous trajectory with manual segment assignment
  n_points <- 200
  df <- data.frame(
    time = seq(0, 20, length.out = n_points),
    x = cumsum(rnorm(n_points, mean = 0, sd = 0.5)),
    y = cumsum(rnorm(n_points, mean = 0, sd = 0.5)),
    segment = c(rep(1, 50), rep(2, 50), rep(3, 50), rep(4, 50))
  )
  
  # Manually create a large time gap between segment 2 and 3
  df$time[101:200] <- df$time[101:200] + 10
  
  gaps <- identify_tracking_gaps(df, window = 2)
  
  expect_true(is.data.frame(gaps))
  expect_true(nrow(gaps) >= 1)
  expect_true(all(c("gap_duration", "gap_start", "gap_end") %in% names(gaps)))
  if (nrow(gaps) > 0) {
    expect_true(any(gaps$gap_duration > 2))
  }
})

test_that("seconds_to_ass_time formats correctly", {
  expect_equal(seconds_to_ass_time(0), "0:00:00.00")
  expect_equal(seconds_to_ass_time(65.5), "0:01:05.50")
  expect_equal(seconds_to_ass_time(3661.25), "1:01:01.25")
})

test_that("calculate_log_likelihood returns numeric", {
  result <- calculate_log_likelihood(residual_gap = 10, sigma_t = 5)
  expect_true(is.numeric(result))
  expect_length(result, 1)
  
  # Log-likelihood should decrease with larger residual gaps
  ll1 <- calculate_log_likelihood(5, 5)
  ll2 <- calculate_log_likelihood(10, 5)
  expect_true(ll1 > ll2)
})

