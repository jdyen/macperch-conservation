# Simulation of Macquarie perch population responses to a
#   suite of targeted interventions at sites spanning
#   the range of this species

# set a seed to ensure reproducibility
set.seed(2021-07-08)

# load some packages to read in data sets
library(qs)

# to work with data sets
library(lubridate)
library(dplyr)
library(tidyr)
library(aae.hydro)

# and to simulate population dynamics
library(aae.pop.templates)

# load some helpers
source("R/actions-helpers.R")
source("R/flow-helpers.R")
source("R/simulate-helpers.R")
source("R/sensitivity-helpers.R")

# define sites
sites <- c("goulburn", "king", "ovens", "mitta", "sevens", "yarra")
catchment <- c(
  "goulburn" = "goulburn", 
  "king" = "ovens",
  "ovens" = "ovens",
  "mitta" = "upper_murray",
  "sevens" = "goulburn",
  "yarra" = "yarra"
)

# define management actions
climates <- c("historical", "drying", "variable")
actions <- list(
  goulburn = c("gene_mixing", "stocking5", "stocking10", "stocking_aspirational", "fishing_regulations", "habitat_restoration", "env_water"),
  king = c("gene_mixing", "stocking5", "stocking10", "stocking_aspirational", "fishing_regulations", "habitat_restoration"),
  ovens = c("gene_mixing", "stocking5", "stocking10", "stocking_aspirational", "fishing_regulations", "habitat_restoration"),
  mitta = c("gene_mixing", "stocking_aspirational", "fishing_regulations"),
  sevens = c("gene_mixing", "stocking5", "stocking10", "stocking_aspirational", "fishing_regulations", "habitat_restoration", "exclude_exotic"),
  yarra = c("gene_mixing", "stocking5", "stocking10", "stocking_aspirational", "fishing_regulations", "habitat_restoration", "env_water")
)

# expand management actions
actions <- lapply(actions, expand_combn)

# don't want to keep multiple stocking actions in a single run
actions <- lapply(actions, remove_conflicts, conflict = c("stocking5", "stocking10", "stocking_aspirational"))

# flatten the actions, sites, and climates to a matrix
actions <- expand_actions(actions, climates = climates)

# load historical flows
# Goulburn R @ Trawool 405201 (temp from 1997-2021, discharge from 1908-2021)
#    Latitude: 37°05'29.8"S, Longitude: 145°12'09.0"E
# King River @ LW Hovell TG 403228 (discharge from 1969-2021) 
#    Latitude: 36°54'32.6"S, Longitude: 146°23'42.0"E
#    King River discharge alternatives: Cheshunt 403227 (discharge 1967-2021) or Edi (discharge 1979-2021)
# Ovens River @ Peechelba 403241 (temperature) and @ Wangaratta 403200 (discharge)
#    Latitude: 36°21'02.4"S, Longitude: 146°19'15.4"E
# Mitta Mitta River @ Hinnomunjie 401203
#    Latitude: 36°56'45.5"S, Longitude: 147°36'20.5"E
# Seven Creeks @ DS of Polly McQuinn Weir 405234 (discharge) and @ Kialla West 405269 (temperature)
#    PMQ: Latitude: 36°53'13.2"S, Longitude: 145°40'57.9"E
#    KW: Latitude: 36°27'24.1"S, Longtitude: 145°23'56.6"E
# Yarra River @ Millgrove 229212 (until 2004, using downloaded data in "data/" for 2004-2020)
#    Latitude: 37°45'09.3"S, Longitude: 145°39'19.7"E
#    Currently ignoring water temperature data (back-filled with Ovens data for convenience)
#    Loading discharge data from saved file due to errors in ratings table in Yarra catchment
discharge_gauges <- c("goulburn" = "405201", "king" = "403228", "ovens" = "403200", "mitta" = "401203", "sevens" = "405234", "yarra" = "229212")
water_temp_gauges <- c("goulburn" = "405201", "king" = "403241", "ovens" = "403241", "mitta" = "401203", "sevens" = "405269", "yarra" = "229212")
site_names <- c("Trawool", "LWHovell", "Peechelba", "Hinnomunjie", "Kialla", "Millgrove")
site_codes <- c(405201, 403241, 403241, 401203, 405269, 229212)
site_lat_dms <- c("37d05'29.8\"S", "36d21'02.4\"S", "36d21'02.4\"S", "36d56'45.5\"S", "36d27'24.1\"S", "37d45'09.3\"S")
site_long_dms <- c("145d12'09.0\"E", "146d19'15.4\"E", "146d19'15.4\"E", "147d36'20.5\"E", "145d23'56.6\"E", "145d39'19.7\"E")
recompile_discharge <- FALSE
if (!any(grepl("_compiled", dir("data"))) | recompile_discharge) {
  
  # download discharge and water temp from WMIS
  discharge <- lapply(
    discharge_gauges,
    fetch_hydro,
    variables = c("discharge"),
    start = "1970-01-01",
    end = "2020-12-31",
    include_missing = TRUE
  )
  water_temp <- lapply(
    water_temp_gauges,
    fetch_hydro,
    variables = c("temp"),
    start = "1970-01-01",
    end = "2020-12-31",
    include_missing = TRUE
  )
  
  # add Yarra River discharge data
  yarra_discharge <- read.csv("data/Millgrove-Daily-River-Flow-Over-The-Years.csv")
  yarra_discharge <- data.frame(
    Date = paste0(
      rep(yarra_discharge$Date, ncol(yarra_discharge) - 1L),
      "-",
      gsub("X", "", rep(colnames(yarra_discharge)[-1], each = nrow(yarra_discharge)))
    ),
    mean_discharge = as.numeric(unlist(yarra_discharge[, -1]))
  )
  yarra_discharge$date_formatted <- parse_date_time(
    yarra_discharge$Date, orders = c("dmy")
  )
  yarra_discharge <- yarra_discharge %>% filter(!is.na(date_formatted))
  discharge$yarra <- discharge$yarra %>%
    left_join(yarra_discharge %>% select("date_formatted", "mean_discharge"),
              by = "date_formatted")
  discharge$yarra <- discharge$yarra %>% mutate(
    value = ifelse(is.na(value), mean_discharge, value),
    quality_code = ifelse(is.na(value), 255, 1)
  ) %>% select(-mean_discharge)

  # remove excess variable codes if provided
  discharge <- lapply(discharge, filter_varcode, priority = "141.00")
  
  # set poor quality code observations to NA
  discharge <- lapply(discharge, filter_qc, threshold = 150)
  water_temp <- lapply(water_temp, filter_qc, threshold = 150)
  
  # which years of discharge are available?
  years_available <- lapply(discharge, check_available)
  
  # fill missing years with resampling
  discharge <- mapply(resample_discharge, x = discharge, available = years_available, SIMPLIFY = FALSE)
  
  # use rolling mean to fill small number of remaining gaps
  discharge <- lapply(
    discharge,
    fill_na_rolling,
    variable = "value",
    recursive = TRUE,
    max_iter = 20
  )
  
  # fill temperature based on nearby air temperatures, with final gaps filled
  #   with rolling means
  water_temp <- mapply(
    impute_temperature,
    data = water_temp,
    target = lapply(discharge, function(x) x$value),
    site = site_codes[seq_along(water_temp)],
    latitude = site_lat_dms[seq_along(water_temp)],
    longitude = site_long_dms[seq_along(water_temp)],
    SIMPLIFY = FALSE
  )

  qsave(discharge, file = "data/discharge_compiled.qs")
  qsave(water_temp, file = "data/water_temp_compiled.qs")
  
} else {
  
  # load pre-compiled versions
  discharge <- qread("data/discharge_compiled.qs")
  water_temp <- qread("data/water_temp_compiled.qs")
  
}

# add 2-3C to Hinnomunjie water temp data
water_temp$mitta <- water_temp$mitta %>% mutate(
  value = value + 3
)

# overwrite Yarra water temp for now (missing data causes errors, ignored below anyway)
water_temp$yarra <- water_temp$ovens

# calculate scaled versions of discharge and water temp (post-1975 and post-1997)
discharge <- lapply(discharge, rescale_discharge)
water_temp <- lapply(water_temp, rescale_discharge)

# add climate change scenarios
discharge <- mapply(
  add_climate_change_scenarios,
  x = discharge, 
  catchment = catchment, 
  MoreArgs = list(
    scenario = "rcp85", reference = 2065, type = "discharge", variable = "value"
  ),
  SIMPLIFY = FALSE
)
water_temp <- mapply(
  add_climate_change_scenarios,
  x = water_temp, 
  catchment = catchment[seq_along(water_temp)], 
  MoreArgs = list(
    scenario = "rcp85", reference = 2065, type = "water_temperature", variable = "value"
  ),
  SIMPLIFY = FALSE
)

# simple function to define years as dry, average,  or wet
#    based on daily average discharge
context_fn <- function(x, thresh = c(100, 1000)) {
  ave <- median(x)
  out <- "dry"
  if (ave > thresh[1])
    out <- "average"
  if (ave > thresh[2])
    out <- "wet"
  out
}

# work out observed and predicted transitions and use these to
#    develop some scenarios
climate_levels <- c("1975", "1997", "rcp85low", "rcp85med", "rcp85high")
trans_mat <- vector("list", length = length(climate_levels))
for (i in seq_along(trans_mat)) {
  trans_mat[[i]] <- define_transition(
    discharge$goulburn[[paste0("value_", climate_levels[i])]], 
    discharge$goulburn$date_formatted, 
    context = context_fn, 
    thresh = quantile(discharge$goulburn$value, prob = c(0.4, 0.6)),
    resolution = calendar_year
  )
}

# define transitions based (loosely) on the post-1975, RCP8.5 high and post-1997 scenarios
#    - "historical" - observed fows, nothing changes from past 50 years
#    - "drying" - more likely to enter and stay in a sequence of dry years
#    - "variable" - more likely to enter dry OR wet sequences
transition_list <- list(
  "historical" = matrix(
    c(0.62, 0.38, 0, 0.23, 0.5, 0.27, 0, 0.38, 0.62), 
    nrow = 3, 
    dimnames = list(c("dry", "average", "wet"), c("dry", "average", "wet"))
  ),
  "drying" = matrix(
    c(0.9, 0.1, 0, 0.4, 0.59, 0.01, 0.1, 0.9, 0), 
    nrow = 3, 
    dimnames = list(c("dry", "average", "wet"), c("dry", "average", "wet"))
  ),
  "variable" = matrix(
    c(0.7, 0.28, 0.02, 0.4, 0.55, 0.05, 0, 0.99, 0.01), 
    nrow = 3, 
    dimnames = list(c("dry", "average", "wet"), c("dry", "average", "wet"))
  )
)

# generate climate change scenarios
env_water <- climate_scenarios <- vector("list", length = length(discharge))
for (i in seq_along(climate_scenarios)) {
  env_water[[i]] <- climate_scenarios[[i]] <- vector("list", length = length(transition_list))
  for (j in seq_along(transition_list)) {
    climate_scenarios[[i]][[j]] <- resample_scenario(
      x = discharge[[i]]$value,
      date = discharge[[i]]$date_formatted,
      y = water_temp[[i]]$value,
      n = 50,
      rep = 20,
      transition = transition_list[[j]],
      context = context_fn,
      thresh = quantile(discharge[[i]]$value, prob = c(0.4, 0.6)),
      resolution = calendar_year
    )
    env_water[[i]][[j]] <- climate_scenarios[[i]][[j]]
    idx <- grepl("value", colnames(climate_scenarios[[i]][[j]])) & 
      !grepl("_", colnames(climate_scenarios[[i]][[j]]))
    env_water[[i]][[j]][, idx] <- apply(
      climate_scenarios[[i]][[j]][, idx],
      2,
      add_environmental_water,
      date = ymd(climate_scenarios[[i]][[j]]$date),
      system = names(discharge)[i]
    )
    
    # TODO: add catch for partial years
    # not sure what's happening here, occurs a few times total, maybe leap years?
    
    # add names for the climate scenarios
    names(env_water[[i]]) <- names(climate_scenarios[[i]]) <- climates
    
  }
}
names(env_water) <- names(climate_scenarios) <- names(discharge)

# save the discharge scenarios for plotting
qsave(climate_scenarios, file = "outputs/simulations/climate-scenarios.qs")

# add environmental water to scenarios based on SWPs
env_water_obs <- discharge
env_water_obs$goulburn <- discharge$goulburn %>% 
  mutate(
    "value_1975" = add_environmental_water(
      value_1975,
      date = date_formatted,
      system = "goulburn"
    ),
    "value_1997" = add_environmental_water(
      value_1997,
      date = date_formatted,
      system = "goulburn"
    ),
    "value_rcp85low" = add_environmental_water(
      value_rcp85low,
      date = date_formatted,
      system = "goulburn"
    ),
    "value_rcp85med" = add_environmental_water(
      value_rcp85med,
      date = date_formatted,
      system = "goulburn"
    ),
    "value_rcp85high" = add_environmental_water(
      value_rcp85high,
      date = date_formatted,
      system = "goulburn"
    )
  )
env_water_obs$yarra <- discharge$yarra %>% 
  mutate(
    "value_1975" = add_environmental_water(
      value_1975,
      date = date_formatted,
      system = "yarra"
    ),
    "value_1997" = add_environmental_water(
      value_1997,
      date = date_formatted,
      system = "yarra"
    ),
    "value_rcp85low" = add_environmental_water(
      value_rcp85low,
      date = date_formatted,
      system = "yarra"
    ),
    "value_rcp85med" = add_environmental_water(
      value_rcp85med,
      date = date_formatted,
      system = "yarra"
    ),
    "value_rcp85high" = add_environmental_water(
      value_rcp85high,
      date = date_formatted,
      system = "yarra"
    )
  )

# calculate e-water contributions
ewater_contributions <- mapply(
  calculate_ewater_contribution,
  env_water_obs,
  discharge
)

# simulate population dynamics
nsim <- 1000
qsave(actions, file = "outputs/simulations/actions_sensitivity.qs")

# set some initial conditions by site
reach_length <- c("goulburn" = 149.54, "king" = 41.07, "ovens" = 102.99 + 26.02, "mitta" = 123.46, "sevens" = 12.38, "yarra" = 100)
density_per_km <- c(12, 8.1, 12, 12, 8.1, 12)
k_by_site <- reach_length * density_per_km
mp <- macquarie_perch(
  k = 1000, 
  system = "river",
  genetic_factor = 1
)
initials <- lapply(
  k_by_site,
  set_initial,
  nsim = nsim,
  matrix = mp$matrix
)

# set up a genetics effect based on current neff
effective_popsize <- c(
  "goulburn" = 330 * k_by_site[1] / k_by_site[3],
  "king" = 330 * k_by_site[2] / k_by_site[3],
  "ovens" = 330,
  "mitta" = 307,
  "sevens" = 33,
  "yarra" = 344
)

genetics_effect <- c(
  "goulburn" = 1.8,
  "king" = 1.8,
  "ovens" = 1.8,
  "mitta" = 1.8,
  "sevens" = 1.8,
  "yarra" = 1.8
)
habitat_effect <- c(
  "goulburn" = 1.1,
  "king" = 1.43,
  "ovens" = 1.05,
  "mitta" = 1.1,
  "sevens" = 1.5,
  "yarra" = 1.1
)
stocking_rate <- c(
  "goulburn" = 30000,
  "king" = 10000,
  "ovens" = 30000,
  "mitta" = 0,
  "sevens" = 5000,
  "yarra" = 30000
)
stocking_rate_aspirational <- c(
  "goulburn" = 60000,
  "king" = 15000,
  "ovens" = 40000,
  "mitta" = 30000,
  "sevens" = 5000,
  "yarra" = 40000
)

# set covariate effects parameters
covar_parameters <- list(
  shift = 150, 
  survival_param = c(0.01, -0.1),
  recruit_param = -0.01,
  spawning_param = c(0, -0.025),
  variability_param = 0.04,
  ctf_param = 0.25,
  ctf_threshold = 3
)

# loop through all actions
nscn <- sum(grepl("_", colnames(climate_scenarios$goulburn$historical)))
if (any(grepl("emps.csv", dir("outputs/tables/")))) {
  estimated_emps <- read.csv("outputs/tables/emps.csv", row.names = 1)
} else {
  estimated_emps <- matrix(NA, nrow = actions$n_required, ncol = nscn)
}
if (any(grepl("risk-curves.qs", dir("outputs/tables/")))) {
  estimated_risk <- qread("outputs/tables/risk-curves.qs")
} else {
  estimated_risk <- vector("list", length = actions$n_required)
}

# run all scenarios for all resampled climates
for (i in seq_len(actions$n_required)) {

  # initialise risk output
  estimated_risk[[i]] <- matrix(NA, nrow = 101, ncol = nscn + 1)
    
  # run through all climate scenarios, including replicates of 
  #    each set of values
  discharge_sub <- climate_scenarios[[actions$site[i]]][[actions$climate[i]]]
  env_water_sub <- env_water[[actions$site[i]]][[actions$climate[i]]]
  discharge_idx <- grepl("value", colnames(discharge_sub)) & !grepl("_", colnames(discharge_sub))
  water_temp_idx <- grepl("_", colnames(discharge_sub))
  for (j in seq_len(sum(discharge_idx))) {
    
    # get settings based on actions
    sim_settings <- get_settings(
      actions = actions$actions[i, ],
      date = discharge_sub$date,
      discharge = discharge_sub[, discharge_idx][, j],
      water_temperature = discharge_sub[, water_temp_idx][, j],
      env_water = env_water_sub[, discharge_idx][, j],
      site = actions$site[i],
      habitat = habitat_effect[actions$site[i]],
      genetics = genetics_effect[actions$site[i]],
      stocking = stocking_rate[actions$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
      stocking_aspirational = stocking_rate_aspirational[actions$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
      fishing = 0.1,
      broodfish = 15,
      exotic = 0.5
    )
    
    # turn off temperature effects in King and Ovens
    if (actions$site[i] %in% c("king", "ovens", "yarra")) {
      sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
    }
    
    # setting carrying capacity
    carrying <- k_by_site[actions$site[i]] * sim_settings$k_factor
    
    # define base population dynamics object
    ntime <- nrow(sim_settings$covars)
    mp <- macquarie_perch(
      k = k_by_site[actions$site[i]],
      system = "river",
      genetic_factor = sim_settings$genetic_factor,
      n = list(
        c(
          rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
          rep(0, ntime - sim_settings$stocking_end)
        ),
        rep(0, ntime),
        rpois(ntime, lambda = sim_settings$n_remove)
      ),
      ntime = rep(ntime, 3), 
      start = rep(1, 3), 
      end = rep(ntime, 3), 
      add = c(TRUE, TRUE, FALSE),
      contributing_min = 0.75,
      p_capture = lapply(sim_settings$p_capture, function(x) x)
    )
    
    # account for effects of exotic predators (redfin)
    mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
    
    # account for additional boosts to early life survival under genetic mixing
    idx <- transition(mp$dynamics$matrix, dims = 1)
    mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
    mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
    
    # extract correct initials
    init_set <- initials[[actions$site[i]]]
    
    # update initials to all zeros if site is new (GLB/KNG) except for
    #   an initial stocking of 20000 fingerlings
    if (actions$site[i] %in% c("goulburn", "king", "ovens")) {
      init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
      init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
    }
    
    # simulate population dynamics
    sims <- simulate(
      mp,
      nsim = nsim,
      init = init_set,
      args = list(
        covariates = c(
          format_covariates(sim_settings$covars),
          list(
            shift = covar_parameters$shift, 
            survival_param = covar_parameters$survival_param,
            recruit_param = covar_parameters$recruit_param,
            spawning_param = covar_parameters$spawning_param,
            variability_param = covar_parameters$variability_param
          )
        )
      ),
      options = list(
        update = update_binomial_leslie,
        tidy_abundances = floor
      )
    )
    
    # calculate and store outputs
    estimated_emps[i, j] <- emps(sims, subset = 3:30, times = 11:ntime)
    if (is.na(estimated_risk[[i]][1, 1]))
      estimated_risk[[i]][, 1]  <- get_cdf(sims, subset = 3:30, times = 11:ntime)$prob
    estimated_risk[[i]][, j + 1] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
    write.csv(estimated_emps, file = "outputs/tables/emps.csv")
    qsave(estimated_risk, file = "outputs/tables/risk-curves.qs")

  }
  
}

# run scenarios to test sensitivity to the effects of genetic mixing
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  genetic_effect = seq(0.5, 2.5, length = nsens),
  site = unique(actions$site)
)
sensitivity$site <- as.character(sensitivity$site)
sensitivity_emps <- rep(NA, times = nrow(sensitivity))
sensitivity_risk <- matrix(NA, nrow = nrow(sensitivity), ncol = 101)
for (i in seq_len(nrow(sensitivity))) {
  
  # set simulation settings
  sim_settings <- get_settings(
    actions = c("gene_mixing", rep("none", 4)),
    date = discharge[[sensitivity$site[i]]]$date_formatted,
    discharge = discharge[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    water_temperature = water_temp[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    env_water = env_water_obs[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    site = sensitivity$site[i],
    habitat = habitat_effect[sensitivity$site[i]],
    genetics = sensitivity$genetic_effect[i],
    stocking = stocking_rate[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    stocking_aspirational = stocking_rate_aspirational[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = 0.1,
    broodfish = 15,
    exotic = 0.5
  )
  
  # turn off temperature effects in King and Ovens
  if (sensitivity$site[i] %in% c("king", "ovens", "yarra")) {
    sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
  }
  
  # setting carrying capacity
  carrying <- k_by_site[sensitivity$site[i]] * sim_settings$k_factor
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  mp <- macquarie_perch(
    k = k_by_site[sensitivity$site[i]],
    system = "river",
    genetic_factor = sim_settings$genetic_factor,
    n = list(
      c(
        rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
        rep(0, ntime - sim_settings$stocking_end)
      ),
      rep(0, ntime),
      rpois(ntime, lambda = sim_settings$n_remove)
    ),
    ntime = rep(ntime, 3), 
    start = rep(1, 3), 
    end = rep(ntime, 3), 
    add = c(TRUE, TRUE, FALSE),
    contributing_min = 0.75,
    p_capture = lapply(sim_settings$p_capture, function(x) x)
  )
  
  # account for effects of exotic predators (redfin)
  mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
  
  # account for additional boosts to early life survival under genetic mixing
  idx <- transition(mp$dynamics$matrix, dims = 1)
  mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
  mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
  
  # extract correct initials
  init_set <- initials[[sensitivity$site[i]]]
  
  # update initials to all zeros if site is new (GLB/KNG) except for
  #   an initial stocking of 20000 fingerlings
  if (sensitivity$site[i] %in% c("goulburn", "king", "ovens")) {
    init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
    init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
  }
  
  # simulate population dynamics
  sims <- simulate(
    mp,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = c(
        format_covariates(sim_settings$covars),
        list(
          shift = covar_parameters$shift, 
          survival_param = covar_parameters$survival_param,
          recruit_param = covar_parameters$recruit_param,
          spawning_param = covar_parameters$spawning_param,
          variability_param = covar_parameters$variability_param
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # calculate and store outputs
  sensitivity_emps[i] <- emps(sims, subset = 3:30, times = 11:ntime)
  sensitivity_risk[i, ] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
  write.csv(sensitivity_emps, file = "outputs/tables/emps-sensitivity-genetics.csv")
  qsave(sensitivity_risk, file = "outputs/tables/risk-curves-sensitivity-genetics.qs")
  
}

# run scenarios to test sensitivity to the effects of habitat rehab
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  habitat_effect = seq(1, 2, length = nsens),
  site = unique(actions$site)
)
sensitivity$site <- as.character(sensitivity$site)
sensitivity_emps <- rep(NA, times = nrow(sensitivity))
sensitivity_risk <- matrix(NA, nrow = nrow(sensitivity), ncol = 101)
for (i in seq_len(nrow(sensitivity))) {
  
  # set simulation settings
  sim_settings <- get_settings(
    actions = c("habitat_restoration", rep("none", 4)),
    date = discharge[[sensitivity$site[i]]]$date_formatted,
    discharge = discharge[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    water_temperature = water_temp[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    env_water = env_water_obs[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    site = sensitivity$site[i],
    habitat = sensitivity$habitat_effect[i],
    genetics = genetics_effect[sensitivity$site[i]],
    stocking = stocking_rate[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    stocking_aspirational = stocking_rate_aspirational[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = 0.1,
    broodfish = 15,
    exotic = 0.5
  )
  
  # turn off temperature effects in King and Ovens
  if (sensitivity$site[i] %in% c("king", "ovens", "yarra")) {
    sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
  }
  
  # setting carrying capacity
  carrying <- k_by_site[sensitivity$site[i]] * sim_settings$k_factor
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  mp <- macquarie_perch(
    k = k_by_site[sensitivity$site[i]],
    system = "river",
    genetic_factor = sim_settings$genetic_factor,
    n = list(
      c(
        rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
        rep(0, ntime - sim_settings$stocking_end)
      ),
      rep(0, ntime),
      rpois(ntime, lambda = sim_settings$n_remove)
    ),
    ntime = rep(ntime, 3), 
    start = rep(1, 3), 
    end = rep(ntime, 3), 
    add = c(TRUE, TRUE, FALSE),
    contributing_min = 0.75,
    p_capture = lapply(sim_settings$p_capture, function(x) x)
  )
  
  # account for effects of exotic predators (redfin)
  mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
  
  # account for additional boosts to early life survival under genetic mixing
  idx <- transition(mp$dynamics$matrix, dims = 1)
  mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
  mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
  
  # extract correct initials
  init_set <- initials[[sensitivity$site[i]]]
  
  # update initials to all zeros if site is new (GLB/KNG) except for
  #   an initial stocking of 20000 fingerlings
  if (sensitivity$site[i] %in% c("goulburn", "king", "ovens")) {
    init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
    init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
  }
  
  # simulate population dynamics
  sims <- simulate(
    mp,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = c(
        format_covariates(sim_settings$covars),
        list(
          shift = covar_parameters$shift, 
          survival_param = covar_parameters$survival_param,
          recruit_param = covar_parameters$recruit_param,
          spawning_param = covar_parameters$spawning_param,
          variability_param = covar_parameters$variability_param
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # calculate and store outputs
  sensitivity_emps[i] <- emps(sims, subset = 3:30, times = 11:ntime)
  sensitivity_risk[i, ] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
  write.csv(sensitivity_emps, file = "outputs/tables/emps-sensitivity-habitat.csv")
  qsave(sensitivity_risk, file = "outputs/tables/risk-curves-sensitivity-habitat.qs")
  
}

# run scenarios to test sensitivity to the effects of stocking rates
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  stocking_rate = seq(0, 100000, length = nsens),
  site = unique(actions$site)
)
sensitivity$site <- as.character(sensitivity$site)
sensitivity_emps <- rep(NA, times = nrow(sensitivity))
sensitivity_risk <- matrix(NA, nrow = nrow(sensitivity), ncol = 101)
for (i in seq_len(nrow(sensitivity))) {
  
  # set simulation settings
  sim_settings <- get_settings(
    actions = c("stocking10", rep("none", 4)),
    date = discharge[[sensitivity$site[i]]]$date_formatted,
    discharge = discharge[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    water_temperature = water_temp[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    env_water = env_water_obs[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    site = sensitivity$site[i],
    habitat = habitat_effect[sensitivity$site[i]],
    genetics = genetics_effect[sensitivity$site[i]],
    stocking = sensitivity$stocking_rate[i] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    stocking_aspirational = stocking_rate_aspirational[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = 0.1,
    broodfish = 15,
    exotic = 0.5
  )
  
  # turn off temperature effects in King and Ovens
  if (sensitivity$site[i] %in% c("king", "ovens", "yarra")) {
    sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
  }
  
  # setting carrying capacity
  carrying <- k_by_site[sensitivity$site[i]] * sim_settings$k_factor
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  mp <- macquarie_perch(
    k = k_by_site[sensitivity$site[i]],
    system = "river",
    genetic_factor = sim_settings$genetic_factor,
    n = list(
      c(
        rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
        rep(0, ntime - sim_settings$stocking_end)
      ),
      rep(0, ntime),
      rpois(ntime, lambda = sim_settings$n_remove)
    ),
    ntime = rep(ntime, 3), 
    start = rep(1, 3), 
    end = rep(ntime, 3), 
    add = c(TRUE, TRUE, FALSE),
    contributing_min = 0.75,
    p_capture = lapply(sim_settings$p_capture, function(x) x)
  )
  
  # account for effects of exotic predators (redfin)
  mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
  
  # account for additional boosts to early life survival under genetic mixing
  idx <- transition(mp$dynamics$matrix, dims = 1)
  mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
  mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
  
  # extract correct initials
  init_set <- initials[[sensitivity$site[i]]]
  
  # update initials to all zeros if site is new (GLB/KNG) except for
  #   an initial stocking of 20000 fingerlings
  if (sensitivity$site[i] %in% c("goulburn", "king", "ovens")) {
    init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
    init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
  }
  
  # simulate population dynamics
  sims <- simulate(
    mp,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = c(
        format_covariates(sim_settings$covars),
        list(
          shift = covar_parameters$shift, 
          survival_param = covar_parameters$survival_param,
          recruit_param = covar_parameters$recruit_param,
          spawning_param = covar_parameters$spawning_param,
          variability_param = covar_parameters$variability_param
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # calculate and store outputs
  sensitivity_emps[i] <- emps(sims, subset = 3:30, times = 11:ntime)
  sensitivity_risk[i, ] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
  write.csv(sensitivity_emps, file = "outputs/tables/emps-sensitivity-stocking.csv")
  qsave(sensitivity_risk, file = "outputs/tables/risk-curves-sensitivity-stocking.qs")
  
}

# run scenarios to test sensitivity to the effects of angling rates
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  fishing_rate = seq(0, 0.5, length = nsens),
  action = c("none", "fishing_regulations"),
  site = unique(actions$site)
)
sensitivity$site <- as.character(sensitivity$site)
sensitivity_emps <- rep(NA, times = nrow(sensitivity))
sensitivity_risk <- matrix(NA, nrow = nrow(sensitivity), ncol = 101)
for (i in seq_len(nrow(sensitivity))) {
  
  # set simulation settings
  sim_settings <- get_settings(
    actions = c(sensitivity$action[i], rep("none", 4)),
    date = discharge[[sensitivity$site[i]]]$date_formatted,
    discharge = discharge[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    water_temperature = water_temp[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    env_water = env_water_obs[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    site = sensitivity$site[i],
    habitat = habitat_effect[sensitivity$site[i]],
    genetics = genetics_effect[sensitivity$site[i]],
    stocking = sensitivity$stocking_rate[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    stocking_aspirational = stocking_rate_aspirational[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = sensitivity$fishing_rate[i],
    broodfish = 15,
    exotic = 0.5
  )
  
  # turn off temperature effects in King and Ovens
  if (sensitivity$site[i] %in% c("king", "ovens", "yarra")) {
    sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
  }
  
  # setting carrying capacity
  carrying <- k_by_site[sensitivity$site[i]] * sim_settings$k_factor
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  mp <- macquarie_perch(
    k = k_by_site[sensitivity$site[i]],
    system = "river",
    genetic_factor = sim_settings$genetic_factor,
    n = list(
      c(
        rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
        rep(0, ntime - sim_settings$stocking_end)
      ),
      rep(0, ntime),
      rpois(ntime, lambda = sim_settings$n_remove)
    ),
    ntime = rep(ntime, 3), 
    start = rep(1, 3), 
    end = rep(ntime, 3), 
    add = c(TRUE, TRUE, FALSE),
    contributing_min = 0.75,
    p_capture = lapply(sim_settings$p_capture, function(x) x)
  )
  
  # account for effects of exotic predators (redfin)
  mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
  
  # account for additional boosts to early life survival under genetic mixing
  idx <- transition(mp$dynamics$matrix, dims = 1)
  mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
  mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
  
  # extract correct initials
  init_set <- initials[[sensitivity$site[i]]]
  
  # update initials to all zeros if site is new (GLB/KNG) except for
  #   an initial stocking of 20000 fingerlings
  if (sensitivity$site[i] %in% c("goulburn", "king", "ovens")) {
    init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
    init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
  }
  
  # simulate population dynamics
  sims <- simulate(
    mp,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = c(
        format_covariates(sim_settings$covars),
        list(
          shift = covar_parameters$shift, 
          survival_param = covar_parameters$survival_param,
          recruit_param = covar_parameters$recruit_param,
          spawning_param = covar_parameters$spawning_param,
          variability_param = covar_parameters$variability_param
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # calculate and store outputs
  sensitivity_emps[i] <- emps(sims, subset = 3:30, times = 11:ntime)
  sensitivity_risk[i, ] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
  write.csv(sensitivity_emps, file = "outputs/tables/emps-sensitivity-angling.csv")
  qsave(sensitivity_risk, file = "outputs/tables/risk-curves-sensitivity-angling.qs")
  
}

# run scenarios to test sensitivity to the effects of exotic species
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  action = c("none", "exclude_exotic"),
  exotic = seq(0, 0.8, length = nsens),
  exclusion_effects = seq(0.8, 1, length = 6),
  site = "sevens"
)
sensitivity$site <- as.character(sensitivity$site)
sensitivity_emps <- rep(NA, times = nrow(sensitivity))
sensitivity_risk <- matrix(NA, nrow = nrow(sensitivity), ncol = 101)
for (i in seq_len(nrow(sensitivity))) {
  
  # set simulation settings
  sim_settings <- get_settings(
    actions = c(sensitivity$action[i], rep("none", 4)),
    date = discharge[[sensitivity$site[i]]]$date_formatted,
    discharge = discharge[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    water_temperature = water_temp[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    env_water = env_water_obs[[sensitivity$site[i]]][[paste0("value_", sensitivity$climate[i])]],
    site = sensitivity$site[i],
    habitat = habitat_effect[sensitivity$site[i]],
    genetics = genetics_effect[sensitivity$site[i]],
    stocking = sensitivity$stocking_rate[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    stocking_aspirational = stocking_rate_aspirational[sensitivity$site[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = 0.1,
    broodfish = 15,
    exotic = sensitivity$exotic[i]
  )
  
  # change exclusion effectiveness for exotic species
  if (sensitivity$action[i] == "exclude_exotic")
    sim_settings$exotic_factor <- sensitivity$exclusion_effects[i]
  
  # turn off temperature effects in King and Ovens
  if (sensitivity$site[i] %in% c("king", "ovens", "yarra")) {
    sim_settings$covars$temperature_effect <- rep(1, nrow(sim_settings$covars))
  }
  
  # setting carrying capacity
  carrying <- k_by_site[sensitivity$site[i]] * sim_settings$k_factor
  
  # define base population dynamics object
  ntime <- nrow(sim_settings$covars)
  mp <- macquarie_perch(
    k = k_by_site[sensitivity$site[i]],
    system = "river",
    genetic_factor = sim_settings$genetic_factor,
    n = list(
      c(
        rpois(sim_settings$stocking_end, lambda = sim_settings$n_stocked),
        rep(0, ntime - sim_settings$stocking_end)
      ),
      rep(0, ntime),
      rpois(ntime, lambda = sim_settings$n_remove)
    ),
    ntime = rep(ntime, 3), 
    start = rep(1, 3), 
    end = rep(ntime, 3), 
    add = c(TRUE, TRUE, FALSE),
    contributing_min = 0.75,
    p_capture = lapply(sim_settings$p_capture, function(x) x)
  )
  
  # account for effects of exotic predators (redfin)
  mp$dynamics$matrix[1, ] <- sim_settings$exotic_factor * mp$dynamics$matrix[1, ]
  
  # account for additional boosts to early life survival under genetic mixing
  idx <- transition(mp$dynamics$matrix, dims = 1)
  mp$dynamics$matrix[idx] <- sim_settings$genetic_factor * mp$dynamics$matrix[idx]
  mp$dynamics$matrix[idx] <- ifelse(mp$dynamics$matrix[idx] > 0.99, 0.99, mp$dynamics$matrix[idx])
  
  # extract correct initials
  init_set <- initials[[sensitivity$site[i]]]
  
  # update initials to all zeros if site is new (GLB/KNG) except for
  #   an initial stocking of 20000 fingerlings
  if (sensitivity$site[i] %in% c("goulburn", "king", "ovens")) {
    init_set <- matrix(0, nrow = nrow(init_set), ncol = ncol(init_set))
    init_set[, 1] <- rpois(n = nrow(init_set), lambda = 20000 * 0.5 * 0.13)
  }
  
  # simulate population dynamics
  sims <- simulate(
    mp,
    nsim = nsim,
    init = init_set,
    args = list(
      covariates = c(
        format_covariates(sim_settings$covars),
        list(
          shift = covar_parameters$shift, 
          survival_param = covar_parameters$survival_param,
          recruit_param = covar_parameters$recruit_param,
          spawning_param = covar_parameters$spawning_param,
          variability_param = covar_parameters$variability_param
        )
      )
    ),
    options = list(
      update = update_binomial_leslie,
      tidy_abundances = floor
    )
  )
  
  # calculate and store outputs
  sensitivity_emps[i] <- emps(sims, subset = 3:30, times = 11:ntime)
  sensitivity_risk[i, ] <- get_cdf(sims, subset = 3:30, times = 11:ntime)$value
  write.csv(sensitivity_emps, file = "outputs/tables/emps-sensitivity-exotic.csv")
  qsave(sensitivity_risk, file = "outputs/tables/risk-curves-sensitivity-exotic.qs")
  
}

# load and plot sensitivity analyses for individual parameters
file_list <- dir("outputs/tables/")
emps_param <- lapply(
  file_list[grepl("emps-sensitivity", file_list)],
  function(x) read.csv(paste0("outputs/tables/", x), row.names = 1)
)
risk_param <- lapply(
  file_list[grepl("risk-curves-sensitivity", file_list)],
  function(x) qread(paste0("outputs/tables/", x))
)

# plot sensitivity to angling effects by site
nsens <- 20
sensitivity <- expand.grid(
  climate = "1975",
  fishing_rate = seq(0, 0.5, length = nsens),
  action = c("none", "fishing_regulations"),
  site = unique(actions$site)
)
idx <- sensitivity$action == "none"
plot_sensitivity(
  emps_param[[1]]$x[idx], 
  risk = risk_param[[1]][idx, ], 
  sens = sensitivity[idx, ], 
  par = "fishing_rate",
  # secondary_group = as.integer(as.factor(sensitivity$action)),
  emps_file = "outputs/figs/FigS10-emps-sens-angling.png",
  risk_file = "outputs/figs/FigS11-risk-sens-angling.png"
)

# plot sensitivity to exotic species
sensitivity <- expand.grid(
  climate = "1975",
  action = c("none", "exclude_exotic"),
  exotic_effect = seq(0, 0.8, length = nsens),
  exclusion_effects = as.character(seq(0.8, 1, length = 6)),
  site = "sevens"
)
plot_sensitivity(
  emps_param[[2]]$x,
  risk = risk_param[[2]], 
  sens = sensitivity, 
  par = "exotic_effect",
  group = "exclusion_effects",
  panel_names = c(
    "0.8" = "80%", 
    "0.84" = "84%",
    "0.88" = "88%",
    "0.92" = "92%",
    "0.96" = "96%",
    "1" = "100%"
  ),
  secondary_group = as.integer(as.factor(sensitivity$action)),
  emps_file = "outputs/figs/FigS12-emps-sens-exotic.png",
  risk_file = "outputs/figs/FigS13-risk-sens-exotic.png"
)


# plot sensitivity to genetic mixing
sensitivity <- expand.grid(
  climate = "1975",
  genetic_effect = seq(0.5, 2.5, length = nsens),
  site = unique(actions$site)
)
plot_sensitivity(
  emps_param[[3]]$x,
  risk = risk_param[[3]],
  sens = sensitivity, 
  par = "genetic_effect",
  emps_file = "outputs/figs/FigS14-emps-sens-genetic.png",
  risk_file = "outputs/figs/FigS15-risk-sens-genetic.png"
)

# plot sensitivity to habitat effects
sensitivity <- expand.grid(
  climate = "1975",
  habitat_effect = seq(1, 2, length = nsens),
  site = unique(actions$site)
)
plot_sensitivity(
  emps_param[[4]]$x, 
  risk = risk_param[[4]], 
  sens = sensitivity, 
  par = "habitat_effect",
  emps_file = "outputs/figs/FigS16-emps-sens-habitat.png",
  risk_file = "outputs/figs/FigS17-risk-sens-habitat.png"
)

# plot sensitivity to stocking rate
sensitivity <- expand.grid(
  climate = "1975",
  stocking_rate = seq(0, 100000, length = nsens),
  site = unique(actions$site)
)
plot_sensitivity(
  emps_param[[5]]$x,
  risk = risk_param[[5]], 
  sens = sensitivity, 
  par = "stocking_rate",
  emps_file = "outputs/figs/FigS18-emps-sens-stocking.png",
  risk_file = "outputs/figs/FigS19-risk-sens-stocking.png"
)

# load main emps and risk results to check sensitivity to climate
cc_emps <- read.csv("outputs/tables/emps.csv", row.names = 1)
cc_risk <- qread("outputs/tables/risk-curves.qs")

# list all simulations so we can load and summarise each
actions <- qread("outputs/simulations/actions_sensitivity.qs")

# extract pr_persist from risk curves, based on n = 100
get_pr_persist_sens <- function(x, thresh) {
  idx <- apply(
    x[, -1],
    2,
    function(y) which.min(abs(y - thresh))
  )
  1 - x[idx, 1]
}
pr_persist <- t(mapply(
  get_pr_persist_sens,
  x = cc_risk,
  thresh = 0.25 * k_by_site[actions$site]
))
pr_persist_summary <- t(apply(pr_persist, 1, quantile, p = c(0.1, 0.5, 0.9)))
emps_summary <- t(apply(cc_emps, 1, quantile, p = c(0.1, 0.5, 0.9)))

# and split pr_persist by system
extract_pop_persist <- function(site, actions, persist, emps) {
  data.frame(
    climate = actions$climate[actions$site == site],
    actions = actions$actions[actions$site == site, ],
    persist_lower = persist[actions$site == site, 1], 
    persist_mid = persist[actions$site == site, 2], 
    persist_upper = persist[actions$site == site, 3], 
    emps_lower = emps[actions$site == site, 1], 
    emps_mid = emps[actions$site == site, 2], 
    emps_upper = emps[actions$site == site, 3]
  )
}
pop_persist <- lapply(
  sites,
  extract_pop_persist,
  actions = actions,
  persist = pr_persist_summary, 
  emps = emps_summary
)
names(pop_persist) <- sites

# population level tables
pop_outcomes <- do.call(
    rbind, lapply(
      pop_persist,
      reorder_benefit
    )
  )

# prepare plots of risk curves
sites_to_plot <- c("mitta", "yarra", "sevens", "ovens", "king", "goulburn")
action_label <- c(
  "none" = "None*",
  "gene_mixing" = "Augmented gene flow", 
  "stocking" = "Stocking", 
  "stocking5" = "Stocking (5 years)", 
  "stocking10" = "Stocking (10 years)",
  "stocking_aspirational" = "Stocking (10 years, high rate)",
  "fishing_regulations" = "Fishing regulations",
  "habitat_restoration" = "Instream habitat",
  "env_water" = "Environmental water",
  "exclude_exotic" = "Exotic species control"
)
climate_label <- c("Post-1975", "Post-1997", "RCP8.5 low", "RCP8.5 medium", "RCP8.5 high")
risk_sites <- c(
  "Lake Dartmouth", 
  "Yarra River",
  "Seven Creeks",
  "Ovens River", 
  "King River",
  "Goulburn River" 
)

# set up generic colour palette for all sites (that differ in actions)
all_actions <- unique(c(actions$actions))
all_actions <- all_actions[c(1:5, 7, 6, 8:9)]
col_pal <- RColorBrewer::brewer.pal(length(all_actions) + 1, "Set3")

# and plot it
for (i in seq_along(climates)) {
  
  # initialise plot
  file <- switch(
    climates[i],
    "historical" = "FigS20",
    "drying" = "FigS21",
    "variable" = "FigS22"
  )
  png(
    file = paste0("outputs/figs/", file, ".png"),
    units = "in",
    width = 8,
    height = (8 / 3) * 2.35,
    pointsize = 12,
    res = 600
  )
  layout(mat = cbind(c(1, 4, 7), c(2, 5, 7), c(3, 6, 7)), heights = c(rep(1, 2), 0.35))
  par(mar = c(4.1, 4.8, 2.1, 1.1))
  
  for (j in seq_along(sites_to_plot)) {
    
    # subset to relevant risk curves
    idx <- actions$site == sites_to_plot[j] &
      actions$climate == climates[i]
    risk_sub <- cc_risk[idx]
    
    # and actions
    actions_sub <- actions$actions[idx, ]
    
    # remove rows with continuous stocking
    idy <- which(actions$actions[idx, 2] == "none")
    idy <- c(idy, nrow(actions_sub))
    risk_sub <- risk_sub[idy]
    actions_sub <- actions_sub[idy, ]
    
    # set up plot
    plot(
      risk_sub[[1]][, 2],
      risk_sub[[1]][, 1], 
      ylim = c(0, 1), 
      las = 1,
      type = "n", 
      bty = "l",
      xlim = range(unlist(lapply(risk_sub, function(x) x[, 2:ncol(x)]))),
      cex.lab = 1.25,
      xlab = ifelse(j %in% c(4:6), "Threshold population size (TPS)", ""),
      ylab = ifelse(j %in% c(1, 4), "Pr(N < TPS)", "")
    )
    
    # add each curve
    col_idx <- match(actions_sub[-nrow(actions_sub), 1], all_actions)
    col_idx <- c(col_idx, length(col_pal))
    for (k in seq_along(risk_sub)) {
      
      if (sites_to_plot[j] == "sevens" & actions_sub[k, 1] == "fishing_regulations") {
        # do not plot
      } else {
        for (m in 2:ncol(risk_sub[[1]]))
          lines(risk_sub[[k]][, m], risk_sub[[1]][, 1], col = col_pal[col_idx[k]], lwd = 2.5, lty = 1)
      }
      
    }
    
    # add climate label
    mtext(risk_sites[j], side = 3, line = 0.1, adj = 0)
    
    # add legend
    if (j == 6) {
      par(mar = rep(0, 4))
      plot(seq(0, 1, length = 10) ~ seq(0, 200, length = 10),
           type = "n",
           xaxt = "n",
           yaxt = "n",
           xlab = "",
           ylab = "",
           bty = "n")
      actions_text <- action_label[all_actions]
      legend_text <- c(actions_text, "All")
      legend(
        x = 100,
        y = 0.5, 
        xjust = 0.5,
        yjust = 0.5,
        ncol = 3,
        legend = legend_text,
        lty = 1,
        lwd = 2, 
        col = col_pal, 
        bty = "n", 
        cex = 1.25,
        y.intersp = 0.75,
        x.intersp = 1,
        xpd = TRUE
      )
    }
    
  }
  
  # shut down plot
  dev.off()
  
}

# barplot of frequency of action inclusion under each climate for each pop
pop_outcomes$sites <- sapply(strsplit(rownames(pop_outcomes), "\\."), function(x) x[1])
sites_to_plot <- c("mitta", "yarra", "sevens", "ovens", "king", "goulburn")
action_freq_mid <- action_freq_lower <- action_freq_upper <- vector("list", length = length(sites_to_plot))
all_actions <- c(
  "pop_establish",
  "gene_mixing",
  "stocking5",          
  "stocking10",
  "stocking_aspirational",
  "habitat_restoration",
  "fishing_regulations",
  "env_water",
  "exclude_exotic"
)

for (i in seq_along(sites_to_plot)) {
  action_freq_mid[[i]] <- action_freq_lower[[i]] <- action_freq_upper[[i]] <- matrix(0, nrow = length(climates), ncol = 10)
  colnames(action_freq_mid[[i]]) <- colnames(action_freq_lower[[i]]) <- colnames(action_freq_upper[[i]]) <- c("n", all_actions)
  for (j in seq_along(climates)) {
    idx <- pop_outcomes$sites == sites_to_plot[i] & pop_outcomes$climate == climates[j]
    x_sub <- pop_outcomes[idx, ]
    freq_table <- get_frequency(x_sub, "emps_mid")
    action_freq_mid[[i]][j, match(names(freq_table), colnames(action_freq_mid[[i]]))] <- freq_table
    freq_table <- get_frequency(x_sub, "emps_lower")
    action_freq_lower[[i]][j, match(names(freq_table), colnames(action_freq_lower[[i]]))] <- freq_table
    freq_table <- get_frequency(x_sub, "emps_upper")
    action_freq_upper[[i]][j, match(names(freq_table), colnames(action_freq_upper[[i]]))] <- freq_table
  }
}

# remove fishing regs from Sevens if they're included
if ("fishing_regulations" %in% colnames(action_freq_mid[[3]]))
  action_freq_mid[[3]][, colnames(action_freq_mid[[3]]) == "fishing_regulations"] <- 0
if ("fishing_regulations" %in% colnames(action_freq_lower[[3]]))
  action_freq_lower[[3]][, colnames(action_freq_lower[[3]]) == "fishing_regulations"] <- 0
if ("fishing_regulations" %in% colnames(action_freq_upper[[3]]))
  action_freq_upper[[3]][, colnames(action_freq_upper[[3]]) == "fishing_regulations"] <- 0

# plot results
png(file = "outputs/figs/FigS23.png", width = 8, height = (8 / 3) * 2.35, units = "in", res = 600, pointsize = 12)
plot_frequency(x = action_freq_mid, sites = sites_to_plot, actions = all_actions)
dev.off()
png(file = "outputs/figs/FigS24.png", width = 8, height = (8 / 3) * 2.35, units = "in", res = 600, pointsize = 12)
plot_frequency(x = action_freq_lower, sites = sites_to_plot, actions = all_actions)
dev.off()
png(file = "outputs/figs/FigS25.png", width = 8, height = (8 / 3) * 2.35, units = "in", res = 600, pointsize = 12)
plot_frequency(x = action_freq_upper, sites = sites_to_plot, actions = all_actions)
dev.off()

# plot the climate scenarios
climate_scenarios <- qread("outputs/simulations/climate-scenarios.qs")
col_pal <- scales::alpha(viridisLite::viridis(20), 0.25)
for (i in seq_along(climate_scenarios)) {
  png(file = paste0("outputs/figs/FigS", 25 + i, ".png"), width = 6, height = 10, units = "in", res = 600, pointsize = 12)
  par(mfrow = c(3, 1), mar = c(3, 5, 1.8, 0.5))
  tmp <- climate_scenarios[[i]]
  plot(
    tmp$historical$value1 ~ tmp$historical$date, 
    type = "n",
    ylim = range(unlist(lapply(tmp, function(x) x[, 2:21]))),
    ylab = "", 
    xlab = "", 
    bty = "l", 
    las = 1
  )
  mtext("Year", side = 1, adj = 0.5, line = 1.9)
  mtext("Discharge (ML/d)", side = 2, adj = 0.5, line = 3.4)
  for (j in seq_len((ncol(tmp$historical) - 1) / 2))
    lines(tmp$historical[, j + 1] ~ tmp$historical$date, col = col_pal[j])
  mtext("Historical climate", side = 3, adj = 1, line = 0)
  plot(
    tmp$drying$value1 ~ tmp$drying$date, 
    type = "n",
    ylim = range(unlist(lapply(tmp, function(x) x[, 2:21]))),
    ylab = "", 
    xlab = "", 
    bty = "l", 
    las = 1
  )
  mtext("Year", side = 1, adj = 0.5, line = 1.9)
  mtext("Discharge (ML/d)", side = 2, adj = 0.5, line = 3.4)
  for (j in seq_len((ncol(tmp$historical) - 1) / 2))
    lines(tmp$drying[, j + 1] ~ tmp$drying$date, col = col_pal[j])
  mtext("Drying climate", side = 3, adj = 1, line = 0)
  plot(
    tmp$variable$value1 ~ tmp$variable$date, 
    type = "n",
    ylim = range(unlist(lapply(tmp, function(x) x[, 2:21]))),
    ylab = "", 
    xlab = "", 
    bty = "l", 
    las = 1
  )
  mtext("Year", side = 1, adj = 0.5, line = 1.9)
  mtext("Discharge (ML/d)", side = 2, adj = 0.5, line = 3.4)
  for (j in seq_len((ncol(tmp$historical) - 1) / 2))
    lines(tmp$variable[, j + 1] ~ tmp$variable$date, col = col_pal[j])
  mtext("Variable climate", side = 3, adj = 1, line = 0)
  dev.off()
}
