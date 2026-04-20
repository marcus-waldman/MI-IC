#' Get the 12 Competing SEM Model Specifications
#'
#' Returns \code{lavaan} model syntax strings for the 12 models in the
#' Bollen et al. (2014) SIM1 design. M1 is the true generating model.
#'
#' @return A named list of lavaan syntax strings.
#' @export
get_sim1_models <- function() {
  models <- list()

  models$M1 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M2 <- '
    F1 =~ 1*y1 + y2 + y3
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M3 <- '
    F1 =~ 1*y1 + y2 + y3
    F2 =~ 1*y5 + y4 + y6
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M4 <- '
    F1 =~ 1*y1 + y2 + y3
    F2 =~ 1*y5 + y4
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M5 <- '
    F1 =~ 1*y1 + y2 + y3
    F2 =~ 1*y5 + y4 + y6
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    y6 ~~ y7
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M6 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    y6 ~~ y7
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M7 <- '
    F1 =~ 1*y1 + y3 + y4 + y5
    F2 =~ 1*y2 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M8 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2 + F1
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M9 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F2 ~ F1
    F3 ~ F2 + F1
    y6 ~~ y7
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M10 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  models$M11 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7 + y8 + y9
    F2 ~ F1
    F1 ~~ F1
    F2 ~~ F2
  '

  models$M12 <- '
    F1 =~ 1*y1 + y2 + y3
    F2 =~ 1*y4 + y5
    F3 =~ 1*y6 + y8
    F4 =~ 1*y7 + y9
    F2 ~ F1
    F3 ~ F2
    F4 ~ F3
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
    F4 ~~ F4
  '

  return(models)
}


#' Generate Saturated (Unstructured) Model Syntax
#'
#' Builds lavaan syntax for a fully saturated model over \code{var_names}:
#' every variable gets its own variance and every pair gets a free
#' covariance. This is the reference model for the chi-square statistics.
#'
#' @param var_names Character vector of variable names (default y1..y9).
#' @return A character string of lavaan syntax.
#' @export
get_saturated_model <- function(var_names = paste0("y", 1:9)) {
  p <- length(var_names)
  lines <- character()
  for (i in seq_len(p)) {
    rhs <- paste(var_names[i:p], collapse = " + ")
    lines <- c(lines, paste0(var_names[i], " ~~ ", rhs))
  }
  paste(lines, collapse = "\n")
}
