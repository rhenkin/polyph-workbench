bnf_lookup <- fread("new_bnf_lkp.csv")
ltc_chapters <- fread("chapters.tsv")

# Load data.tables
gold_patient <- fread("../data/gold_cp_patient.csv")
gold_ltc <- fread("../data/gold_ltc.csv")
gold_cp <- fread("../data/gold_cp.csv")
gold_outcomes <- fread("../data/gold_outcomes.csv")
gold_acute_presc <- fread("../data/gold_acute_presc.csv", nrows = 10)

# Set keys for optimal performance
setkey(gold_patient, patid)
setkey(gold_ltc, patid, term)
setkey(gold_cp, patid, substance)
setkey(gold_outcomes, patid, outcome, eventdate)
setkey(gold_acute_presc, patid)

gold_patient <- gold_patient[!is.na(imd_quintile)]
valid_pats <- gold_patient$patid
gold_cp <- gold_cp[patid %in% valid_pats]
gold_ltc <- gold_ltc[patid %in% valid_pats]

default_excluded_bnf <- list(
    	list(category = "Chapter", value = "Anaesthesia"),
    	list(category = "Chapter", value = "Appliances"),
    	list(category = "Chapter", value = "Dressings"),
    	#list(category = "Chapter", value = "Ear, Nose And Oropharynx"),
    	list(category = "Chapter", value = "Immunological Products & Vaccines"),
    	list(category = "Chapter", value = "Incontinence Appliances"),
    	list(category = "Chapter", value = "Other Drugs And Preparations"),
    	list(category = "Chapter", value = "Stoma Appliances")
    )