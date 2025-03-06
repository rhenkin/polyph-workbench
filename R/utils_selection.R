get_filtered_substances <- function(included_items = NULL, excluded_items = NULL, bnf_lookup) {
  
  if (is.function(included_items)) included_items <- included_items()
  if (is.function(excluded_items)) excluded_items <- excluded_items()
  
  # Helper function to get substances for a category/value pair
  get_substances_for_item <- function(item) {
    category_key <- paste0("BNF_", item$category)
    selected_name <- item$value
    
    # Get all chemical substances that belong to this category value
    matching_substances <- bnf_lookup[
                 bnf_lookup[[category_key]] == selected_name
    ][nchar(BNF_Presentation_Code)==9, BNF_Chemical_Substance]
    return(matching_substances)
  }
  
  # Start with all substances if no inclusions specified
  if (is.null(included_items) || length(included_items) == 0) {
    included_substances <- bnf_lookup[nchar(BNF_Presentation_Code)==9, BNF_Chemical_Substance]
  } else {
    # Get all substances that match any inclusion criteria
    included_substances <- character(0)
    for (item in included_items) {
      matching_substances <- get_substances_for_item(item)
      included_substances <- unique(c(included_substances, matching_substances))
    }
  }
  
  # If there are exclusions, remove them from the included set
  if (!is.null(excluded_items) && length(excluded_items) > 0) {
    excluded_substances <- character(0)
    for (item in excluded_items) {
      matching_substances <- get_substances_for_item(item)
      excluded_substances <- unique(c(excluded_substances, matching_substances))
    }
    
    # Remove excluded substances from included substances
    final_substances <- setdiff(included_substances, excluded_substances)
  } else {
    final_substances <- included_substances
  }
  
  return(final_substances)
}