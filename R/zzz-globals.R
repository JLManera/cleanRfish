#' @keywords internal
#' @importFrom utils tail
#' @importFrom stats time predict
#' @importFrom dplyr filter select mutate arrange group_by summarise left_join bind_rows ungroup rename lag lead n desc all_of
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon labs coord_fixed scale_y_reverse
#' @importFrom tibble tibble
#' @importFrom zoo na.approx rollmean
#' @importFrom signal sgolayfilt
#' @importFrom mgcv gam s predict.gam
NULL

# This file handles package-level imports
# Function-level exports are defined via @export in functions.R
