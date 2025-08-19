test_that("find_path returns invariant shape and columns", {
  set.seed(1)
  df <- data.frame(time = 1:200,
                   x1 = cumsum(rnorm(200)),
                   y1 = cumsum(rnorm(200)))
  df$x1[101] <- df$x1[101] + 50
  df$y1[101] <- df$y1[101] - 40

  out <- find_path(df, connection_method = "linear")
  expect_equal(nrow(out), nrow(df))
  expect_true(all(c("time","x1","y1","original_x1","original_y1") %in% names(out)))
})

test_that("smooth_path keeps row count", {
  set.seed(2)
  df <- data.frame(time = 1:101,
                   x1 = cumsum(rnorm(101)),
                   y1 = cumsum(rnorm(101)))
  cleaned  <- find_path(df, connection_method = "linear")
  smoothed <- smooth_path(cleaned, p = 3, n = 11)
  expect_equal(nrow(smoothed), nrow(df))
  expect_true(all(c("x1","y1") %in% names(smoothed)))
})
