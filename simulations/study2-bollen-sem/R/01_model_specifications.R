# ============================================================================
# Model Specifications for SIM1 (9 indicators, cross-loadings)
# ============================================================================
# Adapted from Bollen et al. (2014) via IC-FIML/model_specifications.r
# 12 competing models: M1 is the true generating model
# ============================================================================

get_sim1_models <- function() {
  models <- list()

  # M1: True model - all cross-loadings included
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

  # M2: Drop y4 on F1
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

  # M3: Drop y4 on F1 and y7 on F2
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

  # M4: Drop all cross-loadings (simple structure)
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

  # M5: Drop cross-loadings, add error correlation y6-y7
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

  # M6: All cross-loadings + error correlation y6-y7 (overspecified)
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

  # M7: Switch loadings (different factor structure)
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

  # M8: Add direct path F1 -> F3 (overspecified)
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

  # M9: All cross-loadings + error corr + direct path (overspecified)
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

  # M10: Drop structural paths (CFA only)
  models$M10 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7
    F3 =~ 1*y8 + y9 + y6 + y7
    F1 ~~ F1
    F2 ~~ F2
    F3 ~~ F3
  '

  # M11: Two factors (merge F2+F3)
  models$M11 <- '
    F1 =~ 1*y1 + y2 + y3 + y4
    F2 =~ 1*y5 + y4 + y6 + y7 + y8 + y9
    F2 ~ F1
    F1 ~~ F1
    F2 ~~ F2
  '

  # M12: Four factors (split F2 and F3)
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
