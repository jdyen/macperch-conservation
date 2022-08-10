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
  date,
  discharge,
  water_temperature,
  env_water,
  climate,
  site,
  habitat = 1.5,
  genetics = 1.8,
  stocking = 20000,
  stocking_aspirational = 40000,
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
  # turned off boost when using top-down DD
  gen_factor <- ifelse("gene_mixing" %in% actions, genetics, gen_factor[site])
  
  exotic_factor <- 1
  if (site == "sevens")
    exotic_factor <- ifelse("exclude_exotic" %in% actions, 1, exotic)
    
  n_stocked <- 0
  if (any(grepl("stocking", actions))) 
    n_stocked <- stocking
  if (any(grepl("stocking_aspirational", actions)))
    n_stocked <- stocking_aspirational
  
  # add environmental water if needed
  if (any(grepl("env_water", actions)))
    discharge <- env_water
  
  # and calculate flow metrics from discharge time series
  lt_median <- calculate(
    value = discharge,
    date = date,
    rescale = NULL,
    resolution = baseline()
  )$metric
  covars <- calculate_flow_metrics(
    x = discharge,
    date = date,
    rescale = lt_median
  )
  
  if (!is.null(water_temperature)) {
    covars$temperature_effect <- calculate_simplified_temp_metrics(
      x = water_temperature,
      date = date,
      bounds = c(0, 2)
    )
    
  }

  # how many years are we dealing with?
  nyear <- nrow(covars)
  
  # set up angling removals as dynamic objects
  p_capture <- ifelse("fishing_regulations" %in% actions, 0, fishing)
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

# make a plot of the species summaries under two/all climate change scenarios
plot_spp_summary <- function(x, nx = 15, clim = "1975", label = "Historical", col_pal = NULL, ...) {
  
  # which values do we want?
  target1 <- paste0("emps_", clim)
  target2 <- paste0("expected_neff_", clim)
  
  # reorder by relevant CC scenario  
  x <- x[order(x[, target1], decreasing = TRUE), ]
  
  # set default colour palette
  if (is.null(col_pal)) 
    col_pal <- rev(RColorBrewer::brewer.pal(4, "Set3"))
  
  # pull out target rows with systems ordered according to MS
  tmp <- as.matrix(x[seq_len(nx), c(4, 1, 6, 3, 2, 5)])
  sys_names <- colnames(tmp)
  rownames(tmp) <- seq_len(nx)
  colnames(tmp) <- seq_len(ncol(tmp))
  
  # make a basic image plot with gridlines (green = include, red = exclude)  
  image(t(tmp[rev(seq_len(nx)), ]), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "o", col = col_pal)
  grid(nx = 6, ny = nx, col = "black", lty = 1)
  
  # add rank labels
  yvals <- seq(0, 1, length = nx)
  axis(2, at = yvals, labels = as.character(rev(seq_len(nx))), las = 1, tick = FALSE)
  axis(2, at = c(-10, 10), labels = c("", ""), las = 1, lwd.ticks = 0)
  mtext("Rank", side = 2, line = 2)
  
  # add system labels
  legend_text <- c(
    "mitta" = "Lake Dartmouth", 
    "yarra" = "Yarra River",
    "sevens" = "Seven Creeks",
    "ovens" = "Ovens River", 
    "king" = "King River",
    "goulburn" = "Goulburn River"
  )
  axis(1, at = seq(0, 1, length = ncol(tmp)), labels = legend_text[sys_names], tick = FALSE, padj = -0.5)
  axis(1, at = c(-10, 10), labels = c("", ""), las = 1, lwd.ticks = 0)
  
  # add EMPS/Ne values
  values <- cbind(x[, target1][seq_len(nx)], x[, target2][seq_len(nx)])
  text(x = 1.2, y = rev(yvals), round(values[, 1]), xpd = TRUE)
  text(x = 1.3, y = rev(yvals), round(values[, 2]), xpd = TRUE)
  text(x = c(1.2, 1.3), y = 1 + 0.11, c("EMPS", "Ne"), xpd = TRUE, font = 2)
  
  # add CC label  
  mtext(label, side = 3, adj = 0, line = 0.6, cex = 1.2)
  
  # return
  out <- NULL
  
}

# helper function to get change in EMPS between best/do-nothing scenarios
get_emps_change <- function(x) {
  clims <- unique(x$climate)
  out <- rep(NA, length(clims))
  for (i in seq_along(clims)) {
    xsub <- x[x$climate == clims[i], ]
    action_tmp <- xsub[, grepl("actions", colnames(xsub))]
    idx <- apply(action_tmp, 1, function(.x) all(.x == "none"))
    out[i] <- max(xsub$emps) - xsub$emps[idx]
  }
  out
}

# helper function to get change in pr persist between best/do-nothing scenarios
get_persist_change <- function(x) {
  clims <- unique(x$climate)
  out <- rep(NA, length(clims))
  for (i in seq_along(clims)) {
    xsub <- x[x$climate == clims[i], ]
    action_tmp <- xsub[, grepl("actions", colnames(xsub))]
    idx <- apply(action_tmp, 1, function(.x) all(.x == "none"))
    out[i] <- max(xsub$persist) - xsub$persist[idx]
  }
  out
}

# helper function to filter species results to best/do-nothing scenarios only
get_spp_filtered <- function(x) {
  clims <- unique(x$climate)
  out <- vector("list", length = length(clims))
  for (i in seq_along(clims)) {
    xsub <- x[x$climate == clims[i], ]
    action_tmp <- xsub[, grepl("actions", colnames(xsub))]
    idx <- apply(action_tmp, 1, function(.x) all(.x == "none"))
    out[[i]] <- rbind(
      xsub[which.max(xsub$emps), ],
      xsub[idx, ]
    )
  }
  out
}

# sort estimates of emps/persistence by persistence then EMPS
reorder_benefit <- function(x) {
  if ("persist_mid" %in% colnames(x)) {
    out <- x[order(x$climate, x$persist_mid, x$emps_mid, decreasing = c(FALSE, TRUE, TRUE)), ]
    out$emps_mid <- round(out$emps_mid)
  } else {
    out <- x[order(x$climate, x$persist, x$emps, decreasing = c(FALSE, TRUE, TRUE)), ]
    out$emps_mid <- round(out$emps)
  }
  out
}

# extract best and do-nothing scenarios for each population and climate
extract_best <- function(x) {
  
  clims <- unique(x$climate)
  
  out <- matrix(NA, nrow = 2 * length(clims), ncol = ncol(x))
  for (i in seq_along(clims)) {
    
    xsub <- x[x$climate == clims[i], ]
    
    do_nothing <- which(apply(xsub[, grepl("actions", colnames(xsub))], 1, function(.x) all(.x == "none")))
    best <- which(xsub$persist == max(xsub$persist))
    if (length(best) > 1) {
      best <- best[xsub$emps[best] == max(xsub$emps[best])]
    }
    
    out[2 * (i - 1) + 1:2, ] <- as.matrix(xsub[c(do_nothing, best), ])
    
  }
  
  out
  
}
