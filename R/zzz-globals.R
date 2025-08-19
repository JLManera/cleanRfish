#' @keywords internal
#' @importFrom utils tail
#' @importFrom data.table :=
#' @importFrom stats time
NULL

utils::globalVariables(c(
  "segment_id","end_x","end_y","end_time","time","start_time","prediction",
  "start_x","pred_x","start_y","pred_y","euclidean_dist","delta_t","distance",
  "x","y","delta_x","delta_y","velocity","med_dist","mad_dist","dist_threshold",
  "is_jump","original_x","original_y","x.x","y.x","x_filled","y_filled",
  "x_smooth","y_smooth"
))
