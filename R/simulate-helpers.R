# functions to support simulations of population dynamics

# function to set initial conditions based on assumed constant proportional
#   survival over all age classes
set_initial <- function(adult = 500, nsim = 1, age = 1:30, adults = 3:max(age), rate = 0.3, matrix = NULL) {
  
  # work out age frequency from steady state if matrix provided,
  #   from constant survival assumption otherwise
  if (is.null(matrix)) {
    initial_age_frequency <- dexp(age, rate = rate)
  } else {
    initial_age_frequency <- Re(eigen(matrix)$vectors[, 1])
  }
  initial_age_frequency <- initial_age_frequency / sum(initial_age_frequency)
  
  # set total abundance from yoy abundance
  total_abundance <- adult / sum(initial_age_frequency[adults])
  
  # return
  matrix(
    rpois(
      nsim * length(initial_age_frequency),
      lambda = total_abundance * initial_age_frequency
    ),
    nrow = nsim,
    byrow = TRUE
  )
  
}

# function to set simulate settings based on current actions
# effects of actions:
#   gene mixing: based on Lutz et al. 2020, 80 % increase in early life survival
#   stocking: add x fingerlings in a given year
get_settings <- function(
  actions,
  discharge,
  water_temperature,
  env_water,
  climate,
  site,
  habitat = 1.5,
  genetics = 1.8,
  stocking = 20000,
  fishing = 0.1,
  broodfish = 15,
  exotic = 0.75
) {
  
  # set up parameters based on actions
  k_factor <- ifelse("habitat_restoration" %in% actions, habitat, 1)
  
  gen_factor <- c(
    "goulburn" = 1.5,
    "king" = 1.5,
    "ovens" = 1.5,
    "mitta" = 1,
    "sevens" = 1,
    "yarra" = 1.7
  )
  gen_factor <- ifelse("gene_mixing" %in% actions, genetics, gen_factor[site])
  
  exotic_factor <- 1
  if (site == "sevens")
    exotic_factor <- ifelse("exclude_exotic" %in% actions, 1, exotic)
    
  n_stocked <- 0
  if (any(grepl("stocking", actions))) 
    n_stocked <- stocking
  
  # add environmental water if needed
  if (any(grepl("env_water", actions)))
    discharge <- env_water
  
  # and calculate flow metrics from discharge time series
  lt_median <- calculate(
    value = discharge[[paste0("value_", climate)]],
    date = discharge$date_formatted,
    rescale = NULL,
    resolution = baseline()
  )$metric
  covars <- calculate_flow_metrics(
    x = discharge[[paste0("value_", climate)]],
    date = discharge$date_formatted,
    rescale = lt_median
  )
  
  if (!is.null(water_temperature)) {
    covars$temperature_effect <- calculate_simplified_temp_metrics(
      x = water_temperature[[paste0("value_", climate)]],
      date = water_temperature$date_formatted
    )
  }

  # how many years are we dealing with?
  nyear <- nrow(covars)
  
  # set up angling removals as dynamic objects
  p_capture <- ifelse("fishing_regs" %in% actions, 0, fishing)
  p_capture <- ifelse(site == "sevens", 0, p_capture)
  p_capture <- rep(p_capture, nyear)
  if (site %in% c("goulburn", "king", "ovens"))
    p_capture[1:9] <- 0
  
  # set up broodfish removals in mitta
  n_remove <- rep(0, nyear)
  if (site == "mitta")
    n_remove <- c(broodfish, rep(broodfish * 2, 4), rep(broodfish, nyear - 5))
  
  # set time period for stocking
  stocking_end <- ifelse("stocking5" %in% actions, 5, 10)
  
  # return named list of settings
  list(
    genetic_factor = gen_factor,
    n_stocked = n_stocked,
    stocking_end = stocking_end,
    n_remove = n_remove,
    p_capture = p_capture,
    covars = covars,
    k_factor = k_factor,
    exotic_factor = exotic_factor
  )
  
}

# calculate risk of falling below any threshold prior
#   to some time point
calculate_risk <- function(sims, year = NULL, start = 2, threshold = NULL, adults = 3:30, n= 100) {
  
  # how many years are there?
  nyear <- dim(sims)[3]
  
  # default to final year if year not specified
  if (is.null(year))
    year <- nyear - 1
  
  # how many adults are there?
  adult_abundance <- apply(sims[, adults, ], c(1, 3), sum)
  
  # what is our sequence of thresholds?
  if (is.null(threshold))
    threshold <- c(0, exp(seq(0, log(max(adult_abundance)), length = n)))
  
  # calculate probability of being below threshold prior to or in year
  out <- matrix(NA, nrow = length(year), ncol = length(threshold))
  for (i in seq_along(year)) {
    out[i, ] <- sapply(
      threshold,
      function(x, y, z) mean(
        apply(y[, start:(z + 1)], 1, function(.x, .y) any(.x < .y), .y = x)
      ),
      y = adult_abundance,
      z = year[i]
    )
  }
  
  # return
  out
  
}

# function to calculate expected minimum population from year
#   start to the final year
calculate_emps <- function(sims, start = 1, adults = 3:30) {
  
  # which years do we care about?
  years <- lapply(start, function(x, end) x:end, end = dim(sims)[3])
  
  # pull out adult abundances
  adult_abundance <- apply(sims[, adults, ], c(1, 3), sum)
  
  # subset to those years and get min pop 
  sapply(years, get_emps, abund = adult_abundance)
  
}

# internal function to get EMPS from a single time period
get_emps <- function(years, abund) {
  
  # what is the minimum pop size per trajectory?
  if (length(years) > 1) {
    
    # look over multiple years
    emps <- apply(
      abund[, years],
      1,
      min
    )
    
  } else {
    
    # otherwise only one year, so min = observed
    emps <- abund[, years]
    
  }
  
  # return
  mean(emps)
  
}

# function to do all of the summaries above from a file name and path
summarise_sims <- function(
  file,
  path = "outputs/simulations/", 
  threshold = seq_len(1000),
  start = 5
) {
  
  # load single simulation  
  tmp <- qread(paste0(path, file))
  
  # calculate min pop size and risk curve
  emps <- calculate_emps(tmp, adults = 3:30, start = start)
  risk <- calculate_risk(tmp, adults = 3:30, start = start, threshold = threshold)
  
  # return
  list(
    emps = emps,
    risk = risk
  )
  
}

calculate_species_benefit <- function(
  population_benefit, action_cost, climate = "1975", neff = NULL, max_action = NULL
) {
  
  # pull out site names
  nsites <- length(population_benefit)
  sites <- names(population_benefit)
  
  # pull out results for climate
  tmp <- lapply(population_benefit, function(x) x[x$climate == climate, ])
  
  # filter down to fewer actions if required
  if (!is.null(max_action)) {
    tmp <- lapply(tmp, set_max_actions, n = max_action)
  }

  # pull out top 5 action combos for each pop
  tmp <- lapply(tmp, function(x) x[order(x$persist, x$emps, decreasing = TRUE), ])
  tmp <- lapply(tmp, function(x) x[1:5, ])
  
  # create all combinations of actions  
  n_actions <- lapply(tmp, function(x) seq_len(nrow(x)))
  idx <- do.call(expand.grid, n_actions)
  
  # pull out actions and persistence for each action combo
  species_persist <- matrix(NA, nrow = nrow(idx), ncol = nsites)
  action_combos <- vector("list", length = nsites)
  for (i in seq_len(nsites)) {
    action_combos[[i]] <- tmp[[i]][idx[, i], grepl("action", colnames(tmp[[i]]))]
    species_persist[, i] <- tmp[[i]]$persist[idx[, i]]
  }
  n_action <- sapply(action_combos, ncol)
  action_combos <- do.call(cbind, action_combos)
  species_persist_mean <- apply(
    species_persist, 1, function(x) exp(mean(log(x)))
  )
  ranked_actions <- action_combos[order(species_persist_mean, decreasing = TRUE), ]
  ranked_actions <- apply(ranked_actions, 2, function(x) as.character(x))
  
  # add informative column names for the action combos
  colnames(ranked_actions) <- rep(sites, n_action)
  
  # calculate other measures
  species_persist_arith_mean <- apply(
    species_persist, 1, mean
  )
  pr_at_least_one <- apply(species_persist, 1, function(x) 1 - prod(1 - x))
  pr_all_three <- apply(species_persist, 1, prod)
  weighted_neff <- NULL
  if (!is.null(neff))
    weighted_neff <- apply(species_persist, 1, function(x) sum(neff * x))

  # combine these
  persistence <- data.frame(
    geometric_mean = species_persist_mean,
    arith_mean = species_persist_arith_mean,
    pr_any_persist = pr_at_least_one,
    pr_all_persist = pr_all_three,
    weighted_neff = weighted_neff
  )

  # add a cost for each action and calculate cost-benefit
  total_cost <- apply(ranked_actions, 1, function(i) sum(action_cost[i]))

  # reorder everything based on best outcomes
  idx <- order(persistence$arith_mean, persistence$geometric_mean, persistence$weighted_neff, persistence$pr_any_persist, decreasing = TRUE)
  
  # return
  list(
    cost = total_cost[idx],
    benefit = persistence[idx, ],
    actions = ranked_actions[idx, ]
  )
  
}

# internal function to keep only those action combos with
#    n or fewer actions
set_max_actions <- function(x, n) {
  idy <- grepl("actions", colnames(x))
  idx <- apply(x[, idy], 1, function(.x, max_n) sum(.x == "none")) >= (sum(idy) - n)
  x[idx, ]
}
