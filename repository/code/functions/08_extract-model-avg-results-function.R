# function to extract pertinent results from model-averaged object

extract_model_avg <- function(model_avg,
                              global_model,
                              hillchao,
                              fixed,
                              full = FALSE,
                              pretty_names = NULL,   # named character vector: c("Age"="Age (yrs)", ...)
                              present_char = ".") {
  
  # ---------------------------
  # 1) Global model terms + coefficient ("variable-level") names
  # ---------------------------
  trm <- terms(global_model)
  term_labels <- attr(trm, "term.labels")
  
  canon_term <- function(x) {
    # For interactions, sort components alphabetically: "B:A" -> "A:B"
    ifelse(grepl(":", x),
           vapply(strsplit(x, ":", fixed = TRUE),
                  function(parts) paste(sort(parts), collapse = ":"),
                  character(1)),
           x)
  }
  
  canon_coef <- function(x) {
    ifelse(grepl(":", x),
           vapply(strsplit(x, ":", fixed = TRUE),
                  function(parts) paste(sort(parts), collapse = ":"),
                  character(1)),
           x)
  }
  
  # model.matrix gives the coefficient-level columns actually used for estimation
  mm <- model.matrix(global_model)
  coef_names <- colnames(mm)
  
  # drop intercept column name if present
  coef_names <- coef_names[coef_names != "(Intercept)"]
  
  # map each coefficient column to its originating term using model.matrix "assign"
  # assign = 0 for intercept, 1..length(term_labels) for terms
  assign_vec <- attr(mm, "assign")
  coef_assign <- assign_vec[colnames(mm) != "(Intercept)"]
  term_for_coef <- term_labels[coef_assign]
  
  mapping_df <- data.frame(
    term = term_for_coef,          # original global-model term (e.g., "Season", "crop.stage", "Age:Treatment")
    coef = coef_names,             # coefficient/level-specific name (e.g., "SeasonLate", "crop.stageosr")
    stringsAsFactors = FALSE
  )
  
  mapping_df$term_canon <- canon_term(mapping_df$term)
  
  # Optional: pretty term names (applied to term, not coef)
  if (!is.null(pretty_names)) {
    mapping_df$Parameter <- mapping_df$coef
    hits <- mapping_df$coef %in% names(pretty_names)
    mapping_df$Parameter[hits] <- unname(pretty_names[mapping_df$coef[hits]])
  } else {
    mapping_df$Parameter <- mapping_df$term
  }
  
  # ---------------------------
  # 2) Term weights (sw): term-level importance
  # ---------------------------
  # sw() should return a named numeric vector with names = term labels
  term_sw <- sw(model_avg)
  names(term_sw) <- canon_term(names(term_sw))
  
  mapping_df$Weight <- NA_real_
  in_sw <- mapping_df$term_canon %in% names(term_sw)
  mapping_df$Weight[in_sw] <- unname(term_sw[mapping_df$term_canon[in_sw]])
  
  
  # ---------------------------
  # 3) Coefficients + SEs + 95% CIs (coef-level), aligned to global coef_names
  # ---------------------------
  est_vec <- stats::coef(model_avg, full = full)
  names(est_vec) <- canon_coef(names(est_vec))
  
  mapping_df$Estimate <- NA_real_
  coef_canon <- canon_coef(mapping_df$coef)
  in_est <- coef_canon %in% names(est_vec)
  mapping_df$Estimate[in_est] <- unname(est_vec[coef_canon[in_est]])
  
  # extract SEs from model average summary
  se_vec <- summary(model_avg)$coefmat.subset[, "Std. Error"]
  names(se_vec) <- canon_coef(names(se_vec))
  
  mapping_df$SE <- NA_real_
  in_se <- coef_canon %in% names(se_vec)
  mapping_df$SE[in_se] <- unname(se_vec[coef_canon[in_se]])
  
  # confint columns often named "2.5 %" and "97.5 %" etc; take first two numeric cols
  
  ci_mat <- confint(model_avg, full = full)
  ci_df <- as.data.frame(ci_mat)
  ci_df$coef <- rownames(ci_df)
  ci_df$coef_canon <- canon_coef(ci_df$coef)
  
  num_cols <- which(vapply(ci_df, is.numeric, logical(1)))
  lwr_col <- num_cols[1]
  upr_col <- num_cols[2]
  
  # merge using canonical coef names
  ci_df2 <- ci_df[, c("coef_canon", names(ci_df)[lwr_col], names(ci_df)[upr_col]), drop = FALSE]
  names(ci_df2) <- c("coef_canon", "Lwr", "Upr")
  
  mapping_df$coef_canon <- coef_canon
  mapping_df <- merge(mapping_df, ci_df2, by = "coef_canon", all.x = TRUE, sort = FALSE)
  
  # Because we merged, we now have duplicate Lwr/Upr if they existed; clean up:
  # Keep merged Lwr/Upr, drop the old placeholders if present
  if ("Lwr.x" %in% names(mapping_df)) {
    mapping_df$Lwr <- mapping_df$Lwr.y
    mapping_df$Upr <- mapping_df$Upr.y
    mapping_df$Lwr.x <- NULL; mapping_df$Lwr.y <- NULL
    mapping_df$Upr.x <- NULL; mapping_df$Upr.y <- NULL
  }
  
  # ---------------------------
  # 4) Add Hill-Chao + Fixed flags; final averaged-results df
  # ---------------------------
  mapping_df$HillChao <- hillchao
  mapping_df$Fixed <- fixed
  
  avg_df <- mapping_df[, c("Parameter", "term", "HillChao", "Fixed",
                           "Weight", "Estimate", "SE", "Lwr", "Upr")]
  
  # ---------------------------
  # 5) Sub-model summary table using coefArray + msTable
  # ---------------------------
  
  ms <- model_avg[["msTable"]]
  ca <- model_avg[["coefArray"]]   # [n_models, something, n_coefs]
  
  n_models <- dim(ca)[1]
  
  # Which slice of coefArray should we use to detect presence?
  # Use the first "statistic" slice; presence = any non-NA among coef columns.
  # (This works because absent terms are stored as NA for that model.)
  ca_slice <- ca[, 1, , drop = FALSE]    # dims: [n_models, 1, n_coefs]
  ca_mat   <- ca_slice[, 1, ]            # dims: [n_models, n_coefs]
  coef_cols_in_ca <- colnames(ca_mat)
  coef_cols_in_ca_canon <- canon_coef(coef_cols_in_ca)
  
  map2 <- mapping_df[, c("term_canon", "coef_canon")]
  map2 <- unique(map2)
  
  # for each term, find its coef_canon set, then match to coefArray canonical columns
  term_labels <- attr(terms(global_model), "term.labels")
  term_labels_canon <- canon_term(term_labels)
  
  presence <- matrix("", nrow = n_models, ncol = length(term_labels),
                     dimnames = list(NULL, term_labels))
  
  for (j in seq_along(term_labels)) {
    tcanon <- term_labels_canon[j]
    coefs_t <- map2$coef_canon[map2$term_canon == tcanon]
    if (length(coefs_t) == 0) next
    
    # which columns in coefArray correspond to this term (by canonical coef names)
    cols <- coef_cols_in_ca[coef_cols_in_ca_canon %in% coefs_t]
    if (length(cols) == 0) next
    
    present_t <- rowSums(!is.na(ca_mat[, cols, drop = FALSE])) > 0
    presence[present_t, j] <- present_char
  }
  
  
  presence_df <- as.data.frame(presence, stringsAsFactors = FALSE)
  
  # ---- stats from msTable ----
  # (column names in your msTable are df, AICc, Delta_AICc (or Delta/delta), weight)
  models_df <- data.frame(
    Fixed = fixed,
    HillChao = hillchao,
    Model = seq_len(n_models),
    presence_df,
    stringsAsFactors = FALSE
  )
  
  # attach model-fit columns if they exist
  models_df$df   <- ms$df
  models_df$AICc <- ms$AICc
  models_df$deltaAICc <- ms$delta
  models_df$weight <- ms$weight
  
  # order by weight (high -> low) if present
  models_df <- models_df[order(models_df$weight, decreasing = TRUE), , drop = FALSE]
  rownames(models_df) <- NULL
  
  return(list("vars"=avg_df, "models"=models_df))
}