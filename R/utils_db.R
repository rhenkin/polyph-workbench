buildWhereConditions <- function(terms = NULL,
                                 start_date = NULL,
                                 eth_group = NULL,
                                 sex = NULL,
                                 imd_quintile = NULL,
                                 outcome = NULL,
                                 param_count_init = 1) {

  #sql <- "WHERE patid IN (SELECT patid FROM gold_cp)"
  sql <- "WHERE 1=1"
  params <- list()
  param_count <- param_count_init

  if (!is.null(terms)) {
    placeholders <- paste(sprintf("$%d", param_count:(param_count + length(terms) - 1)),
                          collapse = ",")
    sql <- paste0(sql, " AND term IN (", placeholders, ")")
    params <- c(params, as.list(terms))
    param_count <- param_count + length(terms)
  }

  if (!is.null(start_date)) {
    sql <- paste0(sql, sprintf(" AND eventdate >= $%d::date", param_count))
    params <- c(params, list(start_date))
    param_count <- param_count + 1
  }

  if (!is.null(eth_group)) {
    placeholders <- paste(sprintf("$%d", param_count:(param_count + length(eth_group) - 1)),
                          collapse = ",")
    sql <- paste0(sql, " AND eth_group IN (", placeholders, ")")
    params <- c(params, as.list(eth_group))
    param_count <- param_count + length(eth_group)
  }

  if (!is.null(sex)) {
    placeholders <- paste(sprintf("$%d", param_count:(param_count + length(sex) - 1)),
                          collapse = ",")
    sql <- paste0(sql, " AND gender IN (", placeholders, ")")
    params <- c(params, as.list(sex))
    param_count <- param_count + length(sex)
  }

  if (!is.null(imd_quintile)) {
    placeholders <- paste(sprintf("$%d", param_count:(param_count + length(imd_quintile) - 1)),
                          collapse = ",")
    sql <- paste0(sql, " AND imd_quintile IN (", placeholders, ")")
    params <- c(params, as.list(imd_quintile))
    param_count <- param_count + length(imd_quintile)
  }

  if (!is.null(outcome)) {
    placeholder <- sprintf("$%d", param_count)
    sql <- paste0(sql, " AND outcome = ", placeholder)
    params <- c(params, as.list(outcome))
    param_count <- param_count + 1
  }

  return(list(sql = sql, params = params))
}

# buildPrescriptionQuery <- function(terms = NULL,
#                                    start_date = NULL,
#                                    eth_group = NULL,
#                                    sex = NULL,
#                                    imd_quintile = NULL,
#                                    outcome = NULL) {
#
#   conditions <- buildWhereConditions(terms, start_date, eth_group, sex, imd_quintile, outcome)
#
#   from_sub_query <- "first_outcomes_ever"
#   if (!is.null(terms)) from_sub_query <- paste(from_sub_query, "JOIN gold_cp_ltc USING (patid)")
#   if (!is.null(eth_group) | !is.null(sex) | !is.null(imd_quintile)) from_sub_query <- paste(from_sub_query, "JOIN gold_cp_patient USING (patid)")
#
#   base_sql <- paste("SELECT patid,substance,start_date,stop_date,duration FROM gold_cp WHERE EXISTS (SELECT 1 FROM", from_sub_query)
#   query <- paste0(base_sql, " ", conditions$sql, " AND gold_cp.patid = first_outcomes_ever.patid)")
#   print(query)
#   return(list(query = query, params = conditions$params))
# }

buildAcutePrescriptionQuery <- function(terms = NULL,
																				start_date = NULL,
																				eth_group = NULL,
																				sex = NULL,
																				imd_quintile = NULL,
																				outcome = NULL) {
	conditions <- buildWhereConditions(terms, start_date, eth_group, sex, imd_quintile, outcome)

	# Start building the subquery
	subquery_tables <- "first_outcomes_ever"
	if (!is.null(terms)) subquery_tables <- paste(subquery_tables, "JOIN gold_cp_ltc USING (patid)")
	if (!is.null(eth_group) | !is.null(sex) | !is.null(imd_quintile)) subquery_tables <- paste(subquery_tables, "JOIN gold_cp_patient USING (patid)")

	# Create the optimized query using JOIN instead of EXISTS
	base_sql <- paste0("SELECT gc.patid, gc.substance, gc.start_date, gc.stop_date, gc.duration FROM gold_acute_presc gc JOIN ( SELECT DISTINCT first_outcomes_ever.patid FROM ", subquery_tables)

	query <- paste0(base_sql, " ", conditions$sql, ") fo ON gc.patid = fo.patid")
	print(query)
	return(list(query = query, params = conditions$params))
}

buildPrescriptionQuery <- function(terms = NULL,
																	 start_date = NULL,
																	 eth_group = NULL,
																	 sex = NULL,
																	 imd_quintile = NULL,
																	 outcome = NULL) {
	conditions <- buildWhereConditions(terms, start_date, eth_group, sex, imd_quintile, outcome)

	# Start building the subquery
	subquery_tables <- "first_outcomes_ever"
	if (!is.null(terms)) subquery_tables <- paste(subquery_tables, "JOIN gold_cp_ltc USING (patid)")
	if (!is.null(eth_group) | !is.null(sex) | !is.null(imd_quintile)) subquery_tables <- paste(subquery_tables, "JOIN gold_cp_patient USING (patid)")

	# Create the optimized query using JOIN instead of EXISTS
	base_sql <- paste0("SELECT gc.patid, gc.substance, gc.start_date, gc.stop_date, gc.duration FROM gold_cp gc JOIN ( SELECT DISTINCT first_outcomes_ever.patid FROM ", subquery_tables)

	query <- paste0(base_sql, " ", conditions$sql, ") fo ON gc.patid = fo.patid")
	print(query)
	return(list(query = query, params = conditions$params))
}

#' @export
buildLtcQuery <- function(terms = NULL,
                          start_date = NULL,
                          eth_group = NULL,
                          sex = NULL,
                          imd_quintile = NULL,
                          outcome = NULL) {

  conditions <- buildWhereConditions(terms, start_date, eth_group, sex, imd_quintile, outcome, param_count_init = 1)

  from_sub_query <- ""
  if (!is.null(eth_group) | !is.null(sex) | !is.null(imd_quintile)) from_sub_query <- paste(from_sub_query, "INNER JOIN gold_cp_patient gcp ON gcp.patid = fo.patid")

  #base_sql <- paste("SELECT patid,eventdate,age_days,term FROM gold_cp_ltc gl WHERE EXISTS (SELECT 1 FROM", from_sub_query)
  #base_sql <- paste("WITH first_outcome_dates AS (SELECT DISTINCT patid, MIN(eventdate) as first_outcome_date FROM first_outcomes_ever WHERE outcome = $1 GROUP BY patid) SELECT gl.patid, gl.eventdate, gl.age_days, gl.term FROM first_outcome_dates fo INNER JOIN gold_cp_ltc gl ON gl.patid = fo.patid AND gl.eventdate <= fo.first_outcome_date")
  base_sql <- "SELECT gl.patid, gl.eventdate, gl.age_days, gl.term FROM first_outcomes_ever fo JOIN gold_cp_ltc gl ON gl.patid = fo.patid AND gl.eventdate <= fo.eventdate"
  base_sql <- paste(base_sql, from_sub_query)
  #conditions$params <- c(outcome, conditions$params)
  query <- paste0(base_sql, " ", conditions$sql)
  print(query)
  return(list(query = query, params = conditions$params))
}

#' @export
buildPatientQuery <- function(terms = NULL,
                              eth_group = NULL,
                              sex = NULL,
                              imd_quintile = NULL,
                              outcome = NULL) {
  conditions <- buildWhereConditions(terms, NULL, eth_group, sex, imd_quintile, outcome)
  if (is.null(terms)) {
    base_sql <- "SELECT DISTINCT patid, dob, gender as sex, eth_group, imd_quintile FROM gold_cp_patient JOIN first_outcomes_ever USING (patid)"
  } else {
    base_sql <- "SELECT DISTINCT patid, dob, gender as sex, eth_group, imd_quintile FROM gold_cp_patient JOIN first_outcomes_ever USING (patid) JOIN gold_cp_ltc USING (patid)"
  }
  query <- paste0(base_sql, " ", conditions$sql)
  print(query)
  return(list(query = query, params = conditions$params))
}

buildOutcomeQuery <- function(outcome) {
  query <- "SELECT patid,eventdate,age_days,outcome as term FROM first_outcomes_ever WHERE outcome = $1"
  return(list(query = query, params = list(outcome)))
}

#' @export
buildOutcomeLtcOptionsQuery <-  function(terms = NULL,
                                         start_date = NULL,
                                         eth_group = NULL,
                                         sex = NULL,
                                         imd_quintile = NULL) {

  conditions <- buildWhereConditions(terms, start_date, eth_group, sex, imd_quintile)

  base_sql <- "SELECT DISTINCT term FROM gold_cp_patient JOIN gold_cp_ltc USING (patid)"

  query <- paste0(base_sql, " ", conditions$sql)
  return(list(query = query, params = conditions$params))
}

#' @export
buildLtcDatesQuery <- function(ltc_name) {
  # query <- "SELECT patid, eventdate FROM ltc WHERE term IN ($1)"
  # return(list(query = query, param = ltc_name))
  paste0("SELECT patid, eventdate, age_days FROM gold_cp_ltc WHERE term = ", ltc_name, " AND patid IN (SELECT patid FROM gold_cp)")
}