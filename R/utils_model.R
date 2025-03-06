

#' @export
calculate_substance_irr <- function(subst_pp_df, patient_data, calc_pvalue = FALSE, control_demog = FALSE) {
  message("Computing IRR models ", Sys.time())
  # Create binary columns using dcast
  binary_dt <- dcast(subst_pp_df, patid ~ substance, fun.aggregate = length)

  # Merge with original pp values
  model_dt <- merge(binary_dt, subst_pp_df[, .(patid, pp, outcome_age, n_ltc)], by = "patid")
  # Remove duplicate pp values (if any)
  model_dt <- unique(model_dt)

  if (control_demog) {
    model_dt <- merge(model_dt, patient_data[ , .(patid, sex, eth_group, imd_quintile)], by = "patid", all.x = TRUE, all.y = FALSE)
    #print(model_dt[1:5])
  }

  # Get substance names (excluding patid and pp columns)
  substances <- setdiff(names(model_dt), c("patid", "pp", "outcome_age", "n_ltc", "sex", "eth_group", "imd_quintile"))

  # Function to run GLM and extract IRR for one substance
  run_glm <- function(sub) {
    n_exposed <- sum(model_dt[[sub]])
    if (control_demog) {
      model_formula <- paste0("pp ~ `", sub,"` + sex + eth_group + imd_quintile")
    }
    else {
      model_formula <- paste0("pp ~ `", sub,"` + outcome_age + n_ltc")
    }
    model_formula <- as.formula(model_formula)

    # Fit negative binomial model
    tryCatch({
      if (var(model_dt$pp, na.rm = TRUE)/mean(model_dt$pp, na.rm = TRUE) == 1) {
        model <- glm(model_formula, data = model_dt, family = "poisson")
        coef <- signif(exp(coef(model)[2]), digits = 4)  # IRR for the substance
        ci <- suppressMessages(signif(exp(confint(model)[2,]), digits = 4))
        return(data.table(
          substance = sub,
          IRR = coef,
          lower_CI = ci[1],
          upper_CI = ci[2],
          p_value = NA_real_,
          n_exposed = n_exposed
        ))
      } else {
        #model <- suppressWarnings(glm.nb(model_formula, data = model_dt, control = glm.control(maxit = 50, trace = FALSE)))
        model <- glm(model_formula, data = model_dt, family = poisson)
        coef <- suppressWarnings(signif(exp(coef(model)[2]), digits = 4))  # IRR for the substance
        ci <- suppressMessages(signif(exp(confint(model)[2,]), digits = 4))
        return(data.table(
          substance = sub,
          IRR = coef,
          lower_CI = ci[1],
          upper_CI = ci[2],
          #p_value = NA_real_,
          n_exposed = n_exposed
        ))
      }
    }, error = function(e) {
      print(paste("Error", e, "with", sub, ". N = ", n_exposed))
      return(data.table(
        substance = sub,
        IRR = NA_real_,
        lower_CI = NA_real_,
        upper_CI = NA_real_,
        #p_value = NA_real_,
        n_exposed = n_exposed
      ))
    }) #, silent = TRUE)


    #if (control_demog) print(summary(model))
    # if (!inherits(model, "try-error")) {
    #     #print(class(model))
    #     #print(vif(model))
    #     # Extract coefficients and confidence intervals
    #     coef <- signif(exp(coef(model)[2]), digits = 4)  # IRR for the substance
    #     ci <- suppressMessages(signif(exp(confint(model)[2,]), digits = 4))
    #
    #     # Get p-value
    #     if (calc_pvalue) {
    #         p_value <- signif(summary(model)$coefficients[2,4], digits = 4)
    #     } else {
    #         p_value = NA_real_
    #     }
    #
    #     return(data.table$data.table(
    #         substance = sub,
    #         IRR = coef,
    #         lower_CI = ci[1],
    #         upper_CI = ci[2],
    #         p_value = p_value,
    #         n_exposed = n_exposed
    #     ))
    # } else {
    #     return(data.table$data.table(
    #         substance = sub,
    #         IRR = NA_real_,
    #         lower_CI = NA_real_,
    #         upper_CI = NA_real_,
    #         p_value = NA_real_,
    #         n_exposed = n_exposed
    #     ))
    # }
  }

  # Run GLM for each substance and combine results
  results <- rbindlist(lapply(substances, run_glm))
  # Add sample sizes
  results[, n_total := nrow(model_dt)]
  results[, pct_sample := signif(n_exposed/n_total, digits = 2)]
  # if (calc_pvalue == TRUE)
  #   results[, p_value := p.adjust(p_value, method = "fdr")]
  #results[ , significance := ifelse(p_value < 0.05, "Significant", "Non-significant")]
  # Order by IRR
  setorder(results, -IRR)
  message("Finished models ", Sys.time())
  return(results)
}

# Helper function to prepare data
prepare_substance_data <- function(subst_pp_df, patient_data, control_demog = FALSE) {
  # Create binary columns using dcast
  binary_dt <- dcast(subst_pp_df, patid ~ substance, fun.aggregate = length)

  # Merge with original pp values
  model_dt <- merge(binary_dt, subst_pp_df[, .(patid, pp, outcome_age, n_ltc)], by = "patid")
  model_dt <- unique(model_dt)

  if (control_demog) {
    model_dt <- merge(model_dt,
                      patient_data[, .(patid, sex, eth_group, imd_quintile)],
                      by = "patid", all.x = TRUE, all.y = FALSE)
  }

  # Get substance names
  substances <- setdiff(names(model_dt),
                        c("patid", "pp", "outcome_age", "n_ltc", "sex", "eth_group", "imd_quintile"))

  return(list(model_dt = model_dt, substances = substances))
}

# Function to calculate IRR for a single substance
calculate_single_substance_irr <- function(substance, model_dt, control_demog = FALSE) {
  n_exposed <- sum(model_dt[[substance]])

  if (control_demog) {
    model_formula <- paste0("pp ~ `", substance,"` + sex + eth_group + imd_quintile")
  } else {
    model_formula <- paste0("pp ~ `", substance,"` + outcome_age")
  }
  model_formula <- as.formula(model_formula)

  tryCatch({
    # Check for Poisson vs negative binomial
    if (var(model_dt$pp, na.rm = TRUE)/mean(model_dt$pp, na.rm = TRUE) == 1) {
      model <- glm(model_formula, data = model_dt, family = "poisson")
    } else {
      model <- suppressWarnings(glm.nb(model_formula, data = model_dt,
                                       control = glm.control(maxit = 50, trace = FALSE)))
    }

    coef <- suppressWarnings(signif(exp(coef(model)[2]), digits = 4))
    ci <- suppressMessages(signif(exp(confint(model)[2,]), digits = 4))
    p_value <- signif(summary(model)$coefficients[2,4], digits = 4)

    result <- data.table(
      substance = substance,
      IRR = coef,
      lower_CI = ci[1],
      upper_CI = ci[2],
      p_value = p_value,
      n_exposed = n_exposed,
      n_total = nrow(model_dt)
    )

    result[, pct_sample := n_exposed/n_total]
    result[, significance := ifelse(p_value < 0.05, "Significant", "Non-significant")]

    return(result)

  }, error = function(e) {
    return(data.table(
      substance = substance,
      IRR = NA_real_,
      lower_CI = NA_real_,
      upper_CI = NA_real_,
      p_value = NA_real_,
      n_exposed = n_exposed,
      n_total = nrow(model_dt),
      pct_sample = n_exposed/nrow(model_dt),
      significance = NA_character_
    ))
  })
}

vif <- function(mod, ...) {
  if (any(is.na(coef(mod))))
    stop ("there are aliased coefficients in the model")
  v <- vcov(mod)
  assign <- attr(model.matrix(mod), "assign")
  if (names(coef(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  }
  else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}

calculate_substance_model <- function(subst_pp_df, patient_data, predictor_vars, model_type = "logistic") {
  message("Computing models ", Sys.time())

  # Create binary columns for substances using dcast
  binary_dt <- dcast(subst_pp_df, patid ~ substance, fun.aggregate = length)
  model_dt <- merge(binary_dt, unique(subst_pp_df[, .(patid, pp, outcome_age, n_ltc)]), by = "patid")

  # Merge with demographic data
  model_dt <- merge(
    model_dt,
    patient_data[, c("patid", setdiff(predictor_vars, "outcome_age")), with = FALSE],
    by = "patid",
    all.x = TRUE,
    all.y = FALSE
  )

  # Remove any duplicate rows
  model_dt <- unique(model_dt)

  # Get substance names (excluding patid and demographic variables)
  substances <- setdiff(names(binary_dt), "patid")

  # Function to run model for one substance
  run_model <- function(substance_name) {
    n_exposed <- sum(model_dt[[substance_name]])

    for (var in predictor_vars) {
      if (is.factor(model_dt[[var]]) || is.character(model_dt[[var]])) {
        model_dt[[var]] <- as.factor(model_dt[[var]])
        # Set most frequent level as reference
        freq_table <- table(model_dt[[var]])
        ref_level <- names(freq_table)[which.max(freq_table)]
        model_dt[[var]] <- relevel(model_dt[[var]], ref = ref_level)
      }
    }

    # Check for small groups in categorical variables
    for (var in predictor_vars) {
      if (is.factor(model_dt[[var]]) || is.character(model_dt[[var]])) {
        group_counts <- table(model_dt[[var]])
        small_groups <- names(group_counts[group_counts < 10])
        if (length(small_groups) > 0) {
          warning(paste("Small groups (<10) in", var, ":", paste(small_groups, collapse=", ")))
          # Either combine small groups or exclude them
          if (length(small_groups) < length(group_counts)/2) {
            model_dt[[var]] <- ifelse(model_dt[[var]] %in% small_groups, "Other", model_dt[[var]])
          } else {
            warning(paste("Excluding", var, "due to too many small groups"))
            predictor_vars <- setdiff(predictor_vars, var)
          }
        }
      }
    }

    # Update formula after potential predictor removal
    model_formula <- paste0(
      "`", substance_name, "` ~ ",
      paste(predictor_vars, collapse = " + ")
    )
    model_formula <- as.formula(model_formula)
    model_dt[, outcome_age := outcome_age/365.25]
    tryCatch({
      if (model_type == "logistic") {
        # Check for rare events
      # Standard logistic regression with increased iterations
          model <- glm(model_formula,
          						 data = model_dt[get(substance_name)<2],
          						 family = binomial(link = "logit"))
        } else if (model_type == "poisson") {
          # Fit Poisson regression for count data
          model <- glm(model_formula, data = model_dt, family = poisson(link = "log"))
        } else {
          stop("Unsupported model type")
        }

        # Extract coefficients and confidence intervals
        coefs <- coef(model)
        ci <- suppressMessages(confint(model))

        # Create results table for all predictors
        results <- data.table(
          substance = substance_name,
          predictor = names(coefs)[-1],  # exclude intercept
          estimate = exp(coefs[-1]),     # convert to OR/IRR
          lower_CI = exp(ci[-1, 1]),
          upper_CI = exp(ci[-1, 2]),
          p_value = summary(model)$coefficients[-1, 4],
          n_exposed = n_exposed,
          n_total = nrow(model_dt)
        )

        return(results)

      }, error = function(e) {
        message(paste("Error", e, "with", substance_name, ". N = ", n_exposed))
        return(data.table(
          substance = substance_name,
          predictor = predictor_vars,
          estimate = NA_real_,
          lower_CI = NA_real_,
          upper_CI = NA_real_,
          p_value = NA_real_,
          n_exposed = n_exposed,
          n_total = nrow(model_dt)
        ))
      })
    }

    # Run model for each substance and combine results
    results <- rbindlist(lapply(substances, run_model))

    # Add percentage of sample
    results[, pct_sample := signif(n_exposed/n_total, digits = 2)]

    # Adjust p-values within each substance
    results[, adj_p_value := p.adjust(p_value, method = "fdr"), by = substance]

    # Format results
    format_model_results <- function(results) {

      # Format estimates with significance
      results[, estimate_fmt := sprintf("%.2f%s",
                                        estimate,
                                        ifelse(adj_p_value < 0.05, "*", ""))]

      # Clean predictor names and add footnote markers
      results[, predictor := gsub("^.*?group|^.*?quintile|^sex", "", predictor)]

      # Reshape to wide format
      wide_results <- dcast(results,
                            substance + n_exposed ~ predictor,
                            value.var = "estimate_fmt")

      setorder(wide_results, -n_exposed)

      return(wide_results)
    }

    results <- format_model_results(results)

    message("Finished models ", Sys.time())
    return(results)
  }