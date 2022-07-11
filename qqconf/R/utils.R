#' Shorthand for two numerical comparisons
#'
#' @param x numeric value
#' @param gte lower bound
#' @param lte upper bound
#'
#' @return boolean
between <- function(x, gte, lte) {

  x >= gte & x <= lte

}
