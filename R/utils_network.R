create_igraph_network <- function(network_dt) {
	# Ensure we have the required columns
	if (!all(c("term1", "term2", "ochiai") %in% names(network_dt))) {
		stop("Missing required columns. Need term1, term2, and ochiai")
	}

	# Create edge list with weights
	edges_df <- data.frame(
		from = network_dt$term1,
		to = network_dt$term2,
		weight = network_dt$ochiai
	)

	# Create igraph object
	g <- graph_from_data_frame(
		d = edges_df,
		directed = FALSE
	)

	# Calculate node metrics
	V(g)$degree <- degree(g)
	V(g)$betweenness <- betweenness(g)
	V(g)$community <- as.numeric(membership(cluster_louvain(g)))

	# Print debug info
	print(paste("Created graph with", vcount(g), "vertices and", ecount(g), "edges"))

	return(g)
}


# Convert igraph object to visNetwork format
prepare_visnetwork <- function(g) {
	# Prepare nodes
	nodes <- data.frame(
		id = as.character(V(g)$name),
		label = as.character(V(g)$name),
		value = as.numeric(V(g)$degree),
		group = as.numeric(V(g)$community),
		title = sprintf("Degree: %d<br>Betweenness: %.2f",
										as.numeric(V(g)$degree),
										as.numeric(V(g)$betweenness)),
		stringsAsFactors = FALSE
	)

	# Prepare edges with explicit weight handling
	edge_df <- data.frame(
		from = as.character(get.edgelist(g)[,1]),
		to = as.character(get.edgelist(g)[,2])
	)

	# Add weights if they exist
	if ("weight" %in% edge.attributes(g)) {
		edge_df$value <- E(g)$weight
		edge_df$title <- sprintf("Weight: %.2f", E(g)$weight)
	} else {
		edge_df$value <- 1
		edge_df$title <- "Weight: 1"
	}

	# Create visualization
	visNetwork(nodes, edge_df) %>%
		visIgraphLayout() %>%
		# visPhysics(solver = "forceAtlas2Based",
		#            forceAtlas2Based = list(gravitationalConstant = -100)) %>%
		visOptions(highlightNearest = TRUE,
							 nodesIdSelection = TRUE) %>%
		visLayout(randomSeed = 123) %>%
		visInteraction(navigationButtons = TRUE)  %>%
		visEdges(smooth = FALSE)  %>%
		visNodes(scaling = list(min = 5, max = 30), font = list(size = 30))
}

calculate_disease_network <- function(dt, min_prevalence = 0.01,
																			filter_method = c("standard", "background", "prevalence"),
																			percentile = 0.9,
																			background_multiplier = 2,
																			exclude_top_n_percent = NULL) {
	# Ensure input is data.table
	setDT(dt)
	filter_method <- match.arg(filter_method)

	# Total patients
	total_patients <- uniqueN(dt$patid)

	# Calculate disease frequencies and prevalence
	disease_freq <- dt[, .(frequency = .N,
												 prevalence = .N/total_patients),
										 by = term]

	# Filter diseases by minimum prevalence
	common_diseases <- disease_freq[prevalence >= min_prevalence, term]

	# Filter the original dataset to include only common diseases
	dt_filtered <- dt[term %in% common_diseases]

	# Optionally exclude highly prevalent diseases
	if(!is.null(exclude_top_n_percent)) {
		prevalence_threshold <- quantile(disease_freq$prevalence,
																		 1 - exclude_top_n_percent/100)
		excluded_diseases <- disease_freq[prevalence > prevalence_threshold, term]
		dt_filtered <- dt_filtered[!term %in% excluded_diseases]
	}

	# Calculate co-occurrences with filtered diseases
	cooccurrence <- dt_filtered[, .(term1 = term, patid)]
	setkey(cooccurrence, patid)

	pairs <- cooccurrence[cooccurrence,
												.(term1 = i.term1, term2 = x.term1),
												on = "patid",
												allow.cartesian = TRUE
	][term1 < term2
	][, .(observed = .N), by = .(term1, term2)]

	# Calculate Ochiai and background correction
	result <- pairs[disease_freq, on = .(term1 = term)
	][, `:=` (n1 = frequency,
						prev1 = prevalence)
	][disease_freq, on = .(term2 = term)
	][, `:=` (n2 = frequency,
						prev2 = prevalence)
	]

	# Add derived calculations
	result[, `:=` (
		expected = total_patients * prev1 * prev2,
		ochiai = observed / sqrt(n1 * n2),
		observed_to_expected = observed / (total_patients * prev1 * prev2),
		relative_difference = (observed - (total_patients * prev1 * prev2)) / (total_patients * prev1 * prev2)
	)]

	# Apply filtering based on method
	if(filter_method == "standard") {
		threshold <- quantile(result$ochiai, percentile, na.rm = TRUE)
		result[, above_threshold := ochiai >= threshold]
	} else if(filter_method == "background") {
		result[, above_threshold := observed_to_expected >= background_multiplier]
	} else if(filter_method == "prevalence") {
		result[, max_prevalence := pmax(prev1, prev2)]
		result[, threshold := quantile(ochiai,
																	 pmin(0.95, 0.8 + max_prevalence),
																	 na.rm = TRUE),
					 by = .(cut(max_prevalence, breaks = 10))]
		result[, above_threshold := ochiai >= threshold]
	}

	return(result)
}