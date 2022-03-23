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
climates <- c("1975", "1997", "rcp85low", "rcp85med", "rcp85high")
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

# add environmental water to scenarios based on SWPs
env_water <- discharge
env_water$goulburn <- discharge$goulburn %>% 
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
env_water$yarra <- discharge$yarra %>% 
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
  env_water,
  discharge
)

# simulate population dynamics
nsim <- 1000
qsave(actions, file = "outputs/simulations/actions_list.qs")

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
for (i in seq_len(actions$n_required)) {
  
  # get settings based on actions
  sim_settings <- get_settings(
    actions = actions$actions[i, ],
    discharge = discharge[[actions$site[i]]],
    water_temperature = water_temp[[actions$site[i]]],
    env_water = env_water[[actions$site[i]]],
    climate = actions$climate[i],
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
  
  # save outputs
  qsave(sims, file = paste0("outputs/simulations/sims_", i, ".qs"))
  
}

# list all simulations so we can load and summarise each
actions <- qread("outputs/simulations/actions_list.qs")
sims_list <- paste0("sims_", seq_len(actions$n_required), ".qs")

# calculate population-level outcomes
sim_sum <- lapply(
  sims_list,
  summarise_sims,
  start = 11,
  threshold = seq_len(5000)
)

# save sim_sum to avoid recalculating each time
qsave(sim_sum, file = "outputs/simulations/sim_sum.qs")

# and re-load saved version
sim_sum <- qread("outputs/simulations/sim_sum.qs")

# extract emps and risk curves
emps <- sapply(sim_sum, function(x) x$emps)
risk_curves <- do.call(
  rbind,
  lapply(sim_sum, function(x) x$risk[1, ])
)

# extract pr_persist from risk curves, based on n = 100
pr_persist <- mapply(
  function(x, y) 1 - x$risk[1, y],
  x = sim_sum, 
  y = floor(0.25 * k_by_site[actions$site])
)

# and split pr_persist by system
extract_pop_persist <- function(site, actions, persist, emps) {
  data.frame(
    climate = actions$climate[actions$site == site],
    actions = actions$actions[actions$site == site, ],
    persist = persist[actions$site == site], 
    emps = emps[actions$site == site]
  )
}
pop_persist <- lapply(
  sites,
  extract_pop_persist,
  actions = actions,
  persist = pr_persist, 
  emps = emps
)
names(pop_persist) <- sites

# population level tables
write.csv(
  do.call(
    rbind, lapply(
      pop_persist,
      reorder_benefit
    )
  ),
  file = "outputs/tables/population_outcomes.csv"
)

# and grab best and do-nothing scenarios (for Table 3 in main text)
table3 <- lapply(pop_persist, extract_best)
table3 <- cbind(
  rep(names(pop_persist), times = sapply(table3, nrow)),
  do.call(rbind, table3)
)
write.csv(table3, file = "outputs/tables/table3.csv")

# prepare plots of trajectories
traj_sites <- c(
  "Lake Dartmouth", 
  "Yarra River",
  "Seven Creeks",
  "Ovens River", 
  "King River",
  "Goulburn River" 
)
run_lengths <- rle(actions$site)
pinch_points <- run_lengths$lengths
site_id <- run_lengths$values
start_index <- c(1, cumsum(pinch_points)[-length(pinch_points)] + 1)
end_index <- cumsum(pinch_points)
sims_to_plot <- list(
  mitta = lapply(1:5, function(x) c(start_index[site_id == "mitta"][x], end_index[site_id == "mitta"][x])),
  yarra = lapply(1:5, function(x) c(start_index[site_id == "yarra"][x], end_index[site_id == "yarra"][x])),
  sevens = lapply(1:5, function(x) c(start_index[site_id == "sevens"][x], end_index[site_id == "sevens"][x])),
  ovens = lapply(1:5, function(x) c(start_index[site_id == "ovens"][x], end_index[site_id == "ovens"][x])),
  king = lapply(1:5, function(x) c(start_index[site_id == "king"][x], end_index[site_id == "king"][x])),
  goulburn = lapply(1:5, function(x) c(start_index[site_id == "goulburn"][x], end_index[site_id == "goulburn"][x]))
)
ylim_traj <- c(3500, 3000, 350, 2700, 700, 4000)
png(file = paste0("outputs/figs/Fig2.png"),
     height = (8 / 3) * 2.35,
     width = 8,
     res = 600,
     units = "in",
     pointsize = 12)
layout(
  mat = matrix(c(seq_along(sims_to_plot), rep(length(sims_to_plot) + 1, 3)), nrow = 3, byrow = TRUE),
  heights = c(rep(1, length(sims_to_plot) / 3), 0.35)
)
par(mar = c(4.1, 4.8, 1.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(5, "Set1")
for (i in seq_along(sims_to_plot)) {
  
  sim_plot <- list(
    best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][1], ".qs")),
    baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][2], ".qs"))
  )
  xplot <- seq_len(dim(sim_plot[[1]])[3]) + 1969
  plot(subset(sim_plot[[1]], subset = 4:30)[1, 1, ] ~ xplot,
       xlim = c(1969, 2019),
       ylim = c(0, ylim_traj[i]),
       type = "n",
       xlab = "Year", ylab = "Abundance", las = 1, bty = "l")
  mtext(traj_sites[i], side = 3, line = -1, adj = 1)
  
  for (j in seq_along(sims_to_plot[[i]])) {
    
    sim_plot <- list(
      best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][1], ".qs")),
      baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][2], ".qs"))
    )
    traj_sub <- sample(seq_len(nsim), size = min(c(30, nsim)), replace = FALSE)
    
    for (w in seq_along(sim_plot)) {
      
      # plot trajectories
      for (k in traj_sub) {
        yplot <- apply(
          subset(sim_plot[[w]], subset = 4:30)[k, , ],
          2, 
          sum
        )
        lines(yplot ~ xplot, lty = w, col = scales::alpha(col_pal[j], 0.05)) 
      }
      
      # add mean lines
      ymean <- apply(
        apply(
          subset(sim_plot[[w]], subset = 4:30),
          c(1, 3),
          sum
        ),
        2,
        mean
      )
      lines(ymean ~ xplot, col = col_pal[j], lty = w)
      
    }
    
  }
  
  # add a legend
  if (i == length(sims_to_plot)) {
    x <- seq(0, 1, by = 0.1)
    plot(seq(0, 1, by = 0.1) ~ x,
         bty = "n", xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", type = "n")
    legend(x = 0.5, y = -0.2, 
           ncol = 2,
           xpd = TRUE,
           legend = c("Historical: do nothing", "Post-1997: do nothing", "RCP8.5 low: do nothing", "RCP8.5 medium: do nothing", "RCP8.5 high: do nothing",
                      "Historical: all interventions", "Post-1997: all interventions", "RCP8.5 low: all interventions", "RCP8.5 medium: all interventions", "RCP8.5 high: all interventions"),
           col = rep(col_pal, times = 2), 
           lty = rep(seq_along(sim_plot), each = 5), 
           lwd = 2, xjust = 0.5, yjust = 0.5)
  }
  
}

# close plotting device
dev.off()

# prepare plots of risk curves
sites_to_plot <- c("mitta", "yarra", "sevens", "ovens", "king", "goulburn")
action_label <- c(
  "none" = "None*",
  "gene_mixing" = "Assisted gene flow", 
  "stocking" = "Stocking", 
  "stocking5" = "Stocking (5 years)", 
  "stocking10" = "Stocking (10 years)",
  "stocking_aspirational" = "Stocking (10 years, high rate)",
  "fishing_regulations" = "Fishing regulations",
  "habitat_restoration" = "Instream habitat",
  "env_water" = "Environmental water",
  "exclude_exotic" = "Exclude exotic species"
)
climate_label <- c("Post-1975", "Post-1997", "RCP8.5 low", "RCP8.5 medium", "RCP8.5 high")
max_thresh <- c(4000, 4000, 50, 500, 250, 4000)
threshold <- lapply(seq_along(unique(actions$site)), function(x) seq_len(max_thresh[x]))
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
col_pal <- RColorBrewer::brewer.pal(length(all_actions), "Set3")

# and plot it
for (i in seq_along(climates)) {
  
  # initialise plot
  file <- switch(
    climates[i],
    "1975" = "FigS4",
    "1997" = "FigS5",
    "rcp85low" = "FigS6",
    "rcp85med" = "Fig3",
    "rcp85high" = "FigS7"
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
    risk_sub <- risk_curves[idx, threshold[[j]]]
    
    # and actions
    actions_sub <- actions$actions[idx, ]
    
    # remove rows with continuous stocking
    idy <- which(actions$actions[idx, 2] == "none")
    idy <- c(idy, nrow(actions_sub))
    risk_sub <- risk_sub[idy, ]
    actions_sub <- actions_sub[idy, ]
    
    # set up plot
    plot(
      risk_sub[1, ] ~ threshold[[j]], 
      ylim = c(0, 1), 
      las = 1,
      type = "n", 
      bty = "l",
      xlim = c(0, max_thresh[j]),
      cex.lab = 1.25,
      xlab = ifelse(j %in% c(4:6), "Threshold population size (TPS)", ""),
      ylab = ifelse(j %in% c(1, 4), "Pr(N < TPS)", "")
    )
    
    # add each curve
    col_idx <- match(actions_sub[-nrow(actions_sub), 1], all_actions)
    col_idx <- c(col_idx, length(col_pal))
    for (k in seq_len(nrow(risk_sub))) {
      
      if (sites_to_plot[j] == "sevens" & actions_sub[k, 1] == "fishing_regulations") {
        # do not plot
      } else {
        lines(risk_sub[k, ] ~ threshold[[j]], col = col_pal[col_idx[k]], lwd = 2.5, lty = 1)
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

# plot covariate effects
# calculate covariates in each system and climate
covars <- vector("list", length = length(sites))
names(covars) <- sites
average_water_temperature <- vector("list", length = length(sites))
names(average_water_temperature) <- sites
water_temp$yarra <- water_temp$ovens
for (i in seq_along(sites)) {
  
  # initialise output vector
  covars[[i]] <- vector("list", length = length(climates))
  names(covars[[i]]) <- climates
  average_water_temperature[[i]] <- vector("list", length = length(climates))
  names(average_water_temperature[[i]]) <- climates
  
  for (j in seq_along(climates)) {
    
    # get settings based on actions
    sim_settings <- get_settings(
      actions = actions$actions[1, ],
      discharge = discharge[[sites[i]]],
      water_temperature = water_temp[[sites[i]]],
      env_water = env_water[[sites[i]]],
      climate = climates[j],
      site = sites[i],
      habitat = habitat_effect[sites[i]],
      genetics = genetics_effect[sites[i]],
      stocking = stocking_rate[sites[i]] * 0.5 * 0.13,    # account for fingerling survival and sex ratio
      fishing = 0.1,
      broodfish = 15,
      exotic = 0.75
    )
    
    # calculate average temperature over the spawning period (Oct-Dec)
    average_water_temperature[[i]][[j]] <- calculate(
      value = water_temp[[sites[i]]][[paste0("value_", climates[j])]],
      date = water_temp[[sites[i]]]$date_formatted,
      rescale = NULL,
      resolution = survey(season = 10:12)
    )$metric
    
    # turn off temperature effects in King, Ovens and Yarra
    if (sites[i] %in% c("king", "ovens", "yarra")) {
      sim_settings$covars$temperature_effect <- ifelse(
        sim_settings$covars$temperature_effect < 1, 
        1,
        sim_settings$covars$temperature_effect
      )
    }
    
    # save covariates
    covars[[i]][[j]] <- sim_settings$covars
    
  }
  
}

# calculate range of covariates over all climates
min_val <- t(sapply(covars, function(x) apply(sapply(x, function(y) apply(y, 2, min)), 1, min)))
max_val <- t(sapply(covars, function(x) apply(sapply(x, function(y) apply(y, 2, max)), 1, max)))
temp_range <- sapply(average_water_temperature, function(x) range(unlist(x)))

# and plot for each river
nplot <- 100
var_names <- colnames(min_val)

# initialise plotting device
png(
  file = paste0("outputs/figs/FigS1.png"),
  width = 8,
  height = (8 / 3) * 2.35,
  units = "in",
  res = 300,
  pointsize = 14
)
layout(mat = cbind(c(1, 4, 7), c(2, 5, 7), c(3, 6, 7)), heights = c(rep(1, 2), 0.35))
par(mar = c(4.1, 4.8, 2.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(6, "RdYlGn")

# pull out formatted variable names
formatted_var_name <- gsub("_", " ", var_names)
formatted_var_name <- gsub(" effect", "", formatted_var_name)
formatted_var_name <- gsub("temperature", "water temperature", formatted_var_name)
formatted_var_name <- gsub("^min", "minimum", formatted_var_name)
formatted_var_name <- paste0(toupper(substr(formatted_var_name, 1, 1)), substr(formatted_var_name, 2, nchar(formatted_var_name)))

# work through each site and variable
ymin <- c(0.4, 0.7, 0.6, 1.0, 0, 0.5)
ymax <- c(1.05, 1.05, 1.05, 1.45, 1, 1.5)
yjitter <- seq(-0.02, 0.02, length = length(sites))
for (j in seq_len(ncol(min_val))) {
  
  # generic sequence of covariate j to set up plot window
  covar_range <-  seq(min(min_val[, j]), max(max_val[, j]), length = nplot)
  
  # truncate CTF effect to low flows only
  if (j == 5)
    covar_range <- seq(min(min_val[, j]), 5, length = nplot)
  
  # plot temperature against actual average temps
  if (j == 6)
    covar_range <- seq(min(temp_range[1, ]), max(temp_range[2, ]), length = nplot)
  
  # set up plot
  plot(
    sample(c(ymin[j], ymax[j]), size = nplot, replace = TRUE) ~ covar_range, 
    ylim = c(ymin[j], ymax[j]), 
    las = 1,
    type = "n", 
    bty = "l",
    xlab = formatted_var_name[j],
    ylab = "Effect"
  )
  
  mtext(letters[j], side = 3, adj = 0, line = 0.5)
  
  for (i in seq_along(sites)) {
    
    # fill empty data.frame with sequence for a single variable
    covar_set <-  seq(min_val[i, j], max_val[i, j], length = nplot)
    
    # truncate CTF effect to low flows only
    if (j == 5)
      covar_set <- seq(min(min_val[i, j]), 5, length = nplot)
    
    mat <- 1
    x <- data.frame(
      average_daily_flow = rep(0, nplot),
      spawning_flow = rep(0, nplot),
      spawning_variability = rep(0, nplot),
      river_height_change = rep(0, nplot),
      min_daily_flow = rep(0, nplot),
      temperature_effect = rep(0, nplot)
    )
    x[var_names[j]] <- covar_set

    # calculate covariate effect and pull out variable j only    
    effect <- extract_covar_effects(
      mat, 
      x, 
      spawning_param = covar_parameters$spawning_param, 
      variability_param = covar_parameters$variability_param,
      recruit_param = covar_parameters$recruit_param,
      shift = covar_parameters$shift, 
      survival_param = covar_parameters$survival_param,
      ctf_param = covar_parameters$ctf_param,
      ctf_threshold = covar_parameters$ctf_threshold
    )
    
    effect <- effect[[var_names[j]]] + yjitter[i]
    
    # plot temperature against actual average temps
    if (j == 6)
      covar_set <- seq(temp_range[1, i], temp_range[2, i], length = nplot)
    
    # add each curve
    lines(effect ~ covar_set, col = col_pal[i], lwd = 2)
    
    
  }
  
  if (j == 6) {
    par(mar = rep(0, 4))
    plot(seq(0, 1, length = 10) ~ seq(0, 200, length = 10),
         type = "n",
         xaxt = "n",
         yaxt = "n",
         xlab = "",
         ylab = "",
         bty = "n")
    legend_text <- c(
      "goulburn" = "Goulburn River at Trawool", 
      "king" = "King River below Lake William Hovell",
      "ovens" = "Ovens River", 
      "mitta" = "Mitta Mitta River and Lake Dartmouth", 
      "sevens" = "Seven Creeks below Polly McQuinn Weir",
      "yarra" = "Yarra River at Warrandyte"
    )
    legend_text <- legend_text[sites]
    legend(
      x = 100,
      y = 0.5, 
      xjust = 0.5,
      yjust = 0.5,
      ncol = 2,
      legend = legend_text,
      lty = 1,
      lwd = 2, 
      col = col_pal, 
      bty = "n", 
      cex = 1.1,
      y.intersp = 0.75,
      x.intersp = 1,
      xpd = TRUE
    )
  }
  
}

dev.off()

# create table of expected change in EMPS under different best vs no actions
emps_diff <- lapply(
  pop_persist,
  get_emps_change
)

# repeat for pr_persist
persist_diff <- lapply(
  pop_persist,
  get_persist_change
)

# filter species outputs to just the best/do-nothing scenarios
spp_filtered <- lapply(
  pop_persist,
  get_spp_filtered
)
spp_filtered <- mapply(
  function(x, y) lapply(x, function(.x) data.frame(.x, neff = y)),
  x = spp_filtered,
  y = effective_popsize,
  SIMPLIFY = FALSE
)
spp_best <- lapply(spp_filtered, function(x) list(t(sapply(x, function(y) c(y$emps[1], y$persist[1] * y$neff[1]))),
                                                  t(sapply(x, function(y) c(y$emps[2], y$persist[2] * y$neff[2])))))
spp_best$goulburn[[2]] <- ifelse(spp_best$goulburn[[2]], 0, 0)
spp_best$king[[2]] <- ifelse(spp_best$king[[2]], 0, 0)

idx <- expand.grid(1:2, 1:2, 1:2, 1:2, 1:2, 1:2)
spp_summary <- vector("list", length = nrow(idx))
for (i in seq_len(nrow(idx))) {
  tmp <- mapply(function(x, y) x[[y]], x = spp_best, y = idx[i, ], SIMPLIFY = FALSE)
  spp_summary[[i]] <- apply(abind::abind(tmp, along = 3), c(1, 2), sum)
  rownames(spp_summary[[i]]) <- climates
  colnames(spp_summary[[i]]) <- c("emps", "expected_neff")
}
idx_info <- -1 * idx + 2
colnames(idx_info) <- names(spp_best)

idx_info <- cbind(
  idx_info,
  t(sapply(spp_summary, function(x) x[1, ])),
  t(sapply(spp_summary, function(x) x[2, ])),
  t(sapply(spp_summary, function(x) x[3, ])),
  t(sapply(spp_summary, function(x) x[4, ])),
  t(sapply(spp_summary, function(x) x[5, ]))
)
colnames(idx_info)[7:ncol(idx_info)] <- paste0(colnames(idx_info)[7:ncol(idx_info)], c("_1975", "_1975", "_1997", "_1997", "_RCPlow", "_RCPlow", "_RCPmed", "_RCPmed", "_RCPhigh", "_RCPhigh"))
idx_info <- idx_info[order(idx_info$emps_1975, decreasing = TRUE), ]

# plot species summaries (main text)
png(file = "outputs/figs/Fig5.png", width = 10, height = 7, units = "in", res = 600, pointsize = 12)
par(mar = c(0.5, 3.1, 3.5, 8.1), mfrow = c(2, 1))
plot_spp_summary(idx_info, clim = "1975", label = "Historical")
plot_spp_summary(idx_info, clim = "RCPmed", label = "RCP8.5 (medium)")
dev.off()

# plot species summaries (supp mat, all climates)
png(file = "outputs/figs/FigS9.png", width = 40/3, height = 7, units = "in", res = 600, pointsize = 12)
par(mar = c(0.5, 4.1, 3.5, 9.1), mfrow = c(3, 2))
cc_label <- c("Historical", "Post-1997", "RCP8.5 (low)", "RCP8.5 (medium)", "RCP8.5 (high)")
cc_scen <- c("1975", "1997", "RCPlow", "RCPmed", "RCPhigh")
for (i in seq_along(cc_scen))
  plot_spp_summary(idx_info, clim = cc_scen[i], label = cc_label[i])
dev.off()

# barplot of frequency of action inclusion under each climate for each pop
pop_outcomes <- read.csv("outputs/tables/population_outcomes.csv")
pop_outcomes$sites <- sapply(strsplit(pop_outcomes$X, "\\."), function(x) x[1])
sites_to_plot <- c("mitta", "yarra", "sevens", "ovens", "king", "goulburn")
action_freq <- vector("list", length = length(sites_to_plot))
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
  action_freq[[i]] <- matrix(0, nrow = length(climates), ncol = 10)
  colnames(action_freq[[i]]) <- c("n", all_actions)
  for (j in seq_along(climates)) {
    idx <- pop_outcomes$sites == sites_to_plot[i] & pop_outcomes$climate == climates[j]
    x_sub <- pop_outcomes[idx, ]
    idy <- x_sub$emps >= (0.98 * max(x_sub$emps))
    freq_table <- table(unlist(x_sub[idy, grepl("actions", colnames(x_sub))])) / sum(idy)
    if (sites_to_plot[i] %in% c("king", "goulburn")) {
      freq_table <- c(1.0, freq_table)
      names(freq_table)[1] <- "pop_establish"
    }
    freq_table <- c(sum(idy), freq_table)
    names(freq_table)[1] <- "n"
    freq_table <- freq_table[names(freq_table) != "none"]
    action_freq[[i]][j, match(names(freq_table), colnames(action_freq[[i]]))] <- freq_table
  }
}

legend_text <- c(
  "mitta" = "Lake Dartmouth", 
  "yarra" = "Yarra River",
  "sevens" = "Seven Creeks",
  "ovens" = "Ovens River", 
  "king" = "King River",
  "goulburn" = "Goulburn River"
)

# remove fishing regs from Sevens if they're included
if ("fishing_regulations" %in% colnames(action_freq[[3]])) {
  action_freq[[3]][, colnames(action_freq[[3]]) == "fishing_regulations"] <- 0
}

# plot results
png(file = "outputs/figs/FigS8.png", width = 8, height = (8 / 3) * 2.35, units = "in", res = 600, pointsize = 12)

# set layout
laymat <- matrix(c(seq_along(sites_to_plot), rep(length(sites_to_plot) + 1, 3)), byrow = TRUE, ncol = 3)
layout(laymat, heights = c(rep(1, 2), 0.35))

# set plot margins
old_mar <- par()$mar
par(mar = c(6.1, 4.1, 2.1, 0.5))

col_pal <- RColorBrewer::brewer.pal(length(all_actions), "Set3")

for (i in seq_along(sites_to_plot)) {
  
  # sort action_freq
  act_tmp <- action_freq[[i]][, order(colnames(action_freq[[i]]))]
  
  # remove dud columns
  idx <- !colnames(act_tmp) %in% c("none", "n", "null")
  act_tmp <- act_tmp[, idx]
  act_tmp <- act_tmp[, match(all_actions, colnames(act_tmp))]
  
  # match actions and colours
  idy <- match(colnames(act_tmp), all_actions)  

  # plot it
  rownames(act_tmp) <- c("Historical", "Post-1997", "RCP8.5 low", "RCP8.5 med", "RCP8.5 high")
  barplot(
    t(as.matrix(act_tmp)),
    beside = TRUE,
    space = c(0, 2),
    col = col_pal[idy],
    cex.names = 0.95,
    las = 2,
    ylim = c(0, 1),
    xlab = "",
    ylab = ifelse(i %in% c(1, 4), "Proportional inclusion", "")
  )
  
  # add a label
  mtext(legend_text[sites_to_plot[i]], side = 3, adj = 0, line = 0.4, cex = 1)
  
}

# set plot margins for legend panel and add legend
action_label <- c(
  "pop_establish" = "Establish population",
  "gene_mixing" = "Assisted gene flow", 
  "stocking" = "Stocking", 
  "stocking5" = "Stocking (5 years)", 
  "stocking10" = "Stocking (10 years)",
  "stocking_aspirational" = "Stocking (10 years, high rate)",
  "fishing_regulations" = "Fishing regulations",
  "habitat_restoration" = "Instream habitat",
  "env_water" = "Environmental water",
  "exclude_exotic" = "Exclude exotic species"
)
par(mar = rep(0, 4))
plot(c(0, 1) ~ 1, bty = "n", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend(x = "center", fill = col_pal, legend = action_label[all_actions], cex = 1.2, bty = "n", horiz = FALSE, ncol = 3)

dev.off()

# plot results
png(file = "outputs/figs/Fig4.png", width = 6, height = (2.2 / 2) * 6, units = "in", res = 600, pointsize = 12)

# set layout
laymat <- matrix(c(seq_along(climates[c(1, 4)]), length(climates[c(1, 4)]) + 1), ncol = 1)
layout(laymat, heights = c(1, 1, 0.2))

# set plot margins
par(mar = c(3.2, 4.1, 2.3, 0.5))

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
col_pal <- RColorBrewer::brewer.pal(length(all_actions), "Set3")

climate_label <- c(
  "1975" = "Historical",
  "1997" = "Post-1997",
  "rcp85low" = "RCP8.5 (low)",
  "rcp85med" = "RCP8.5 (medium)",
  "rcp85high" = "RCP8.5 (high)"
)
legend_text_reduced <- c(
  "goulburn" = "Goulburn River", 
  "king" = "King River",
  "ovens" = "Ovens River", 
  "mitta" = "Lake Dartmouth", 
  "sevens" = "Seven Creeks",
  "yarra" = "Yarra River"
)
for (j in c(1, 4)) {
  
  act_combined <- matrix(NA, nrow = length(sites_to_plot), ncol = length(all_actions))
  colnames(act_combined) <- all_actions
  rownames(act_combined) <- sites_to_plot
  for (i in seq_along(sites_to_plot))
    act_combined[i, ] <- action_freq[[i]][j, match(all_actions, colnames(action_freq[[i]]))]
  
  rownames(act_combined) <- legend_text_reduced[sites_to_plot]
  barplot(
    t(as.matrix(act_combined)),
    beside = TRUE,
    space = c(0, 2),
    col = col_pal,
    las = 1,
    cex.names = 1.1,
    xlab = "",
    ylab = "",
    ylim = c(0, 1)
  )
  
  # add a label
  mtext("Proportional inclusion", side = 2, adj = 0.5, line = 2.7, cex = 1)
  mtext(climate_label[climates[j]], side = 3, adj = 0, line = 0.6, cex = 1.1)
  
}

# set plot margins for legend panel and add legend
par(mar = rep(0, 4))
plot(c(0, 1) ~ 1, bty = "n", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
legend(x = "center", fill = col_pal, legend = action_label[all_actions], cex = 1, bty = "n", horiz = FALSE, ncol = 3)

# reset plotting margins and layout
par(mar = old_mar)
layout(matrix(1))

dev.off()

# prepare plots of age distributions
traj_sites <- c(
  "Lake Dartmouth", 
  "Yarra River",
  "Seven Creeks",
  "Ovens River", 
  "King River",
  "Goulburn River" 
)
run_lengths <- rle(actions$site)
pinch_points <- run_lengths$lengths
site_id <- run_lengths$values
start_index <- c(1, cumsum(pinch_points)[-length(pinch_points)] + 1)
end_index <- cumsum(pinch_points)
sims_to_plot <- list(
  mitta = lapply(1:5, function(x) c(start_index[site_id == "mitta"][x], end_index[site_id == "mitta"][x])),
  yarra = lapply(1:5, function(x) c(start_index[site_id == "yarra"][x], end_index[site_id == "yarra"][x])),
  sevens = lapply(1:5, function(x) c(start_index[site_id == "sevens"][x], end_index[site_id == "sevens"][x])),
  ovens = lapply(1:5, function(x) c(start_index[site_id == "ovens"][x], end_index[site_id == "ovens"][x])),
  king = lapply(1:5, function(x) c(start_index[site_id == "king"][x], end_index[site_id == "king"][x])),
  goulburn = lapply(1:5, function(x) c(start_index[site_id == "goulburn"][x], end_index[site_id == "goulburn"][x]))
)
png(file = paste0("outputs/figs/FigS2.png"),
    height = (8 / 3) * 2.35,
    width = 8,
    res = 600,
    units = "in",
    pointsize = 12)
layout(
  mat = matrix(c(seq_along(sims_to_plot), rep(length(sims_to_plot) + 1, 3)), nrow = 3, byrow = TRUE),
  heights = c(rep(1, length(sims_to_plot) / 3), 0.35)
)
par(mar = c(4.1, 4.8, 1.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(5, "Set1")
ylim_traj <- c(1500, 1200, 50, 1500, 350, 1800)
for (i in seq_along(sims_to_plot)) {
  
  sim_plot <- list(
    best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][1], ".qs")),
    baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][2], ".qs"))
  )
  xplot <- 3:30
  yplot <- apply(sim_plot[[1]][, xplot, 51], 2, median)
  plot(yplot ~ xplot,
       ylim = c(0, ylim_traj[i]),
       type = "n",
       xlab = "Age class", ylab = "Abundance", las = 1, bty = "l")
  mtext(traj_sites[i], side = 3, line = -1, adj = 1)
  
  for (j in seq_along(sims_to_plot[[i]])) {
    
    sim_plot <- list(
      best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][1], ".qs")),
      baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][2], ".qs"))
    )
    traj_sub <- sample(seq_len(nsim), size = min(c(30, nsim)), replace = FALSE)
    
    for (w in seq_along(sim_plot)) {
      
      # plot trajectories
      for (k in traj_sub) {
        yplot <- sim_plot[[w]][k, xplot, 51]
        lines(yplot ~ xplot, lty = w, col = scales::alpha(col_pal[j], 0.05)) 
      }
      
      # add mean lines
      ymean <- apply(sim_plot[[w]][, xplot, 51], 2, median)
      lines(ymean ~ xplot, col = col_pal[j], lty = w)
      
    }
    
  }
  
  # add a legend
  if (i == length(sims_to_plot)) {
    x <- seq(0, 1, by = 0.1)
    plot(seq(0, 1, by = 0.1) ~ x,
         bty = "n", xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", type = "n")
    legend(x = 0.5, y = -0.2, 
           ncol = 2,
           xpd = TRUE,
           legend = c("Historical: do nothing", "Post-1997: do nothing", "RCP8.5 low: do nothing", "RCP8.5 medium: do nothing", "RCP8.5 high: do nothing",
                      "Historical: all interventions", "Post-1997: all interventions", "RCP8.5 low: all interventions", "RCP8.5 medium: all interventions", "RCP8.5 high: all interventions"),
           col = rep(col_pal, times = 2), 
           lty = rep(seq_along(sim_plot), each = 5), 
           lwd = 2, xjust = 0.5, yjust = 0.5)
  }
  
}

# close plotting device
dev.off()

png(file = paste0("outputs/figs/FigS3.png"),
    height = (8 / 3) * 2.35,
    width = 8,
    res = 600,
    units = "in",
    pointsize = 12)
layout(
  mat = matrix(c(seq_along(sims_to_plot), rep(length(sims_to_plot) + 1, 3)), nrow = 3, byrow = TRUE),
  heights = c(rep(1, length(sims_to_plot) / 3), 0.35)
)
par(mar = c(4.1, 4.8, 1.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(5, "Set1")
ylim_traj <- c(150, 120, 5, 150, 40, 150)
for (i in seq_along(sims_to_plot)) {
  
  sim_plot <- list(
    best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][1], ".qs")),
    baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[1]][2], ".qs"))
  )
  xplot <- 10:30
  yplot <- apply(sim_plot[[1]][, xplot, 51], 2, median)
  plot(yplot ~ xplot,
       ylim = c(0, ylim_traj[i]),
       type = "n",
       xlab = "Age class", ylab = "Abundance", las = 1, bty = "l")
  mtext(traj_sites[i], side = 3, line = -1, adj = 1)
  
  for (j in seq_along(sims_to_plot[[i]])) {
    
    sim_plot <- list(
      best = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][1], ".qs")),
      baseline = qread(paste0("outputs/simulations/sims_", sims_to_plot[[i]][[j]][2], ".qs"))
    )
    traj_sub <- sample(seq_len(nsim), size = min(c(30, nsim)), replace = FALSE)
    
    for (w in seq_along(sim_plot)) {
      
      # plot trajectories
      for (k in traj_sub) {
        yplot <- sim_plot[[w]][k, xplot, 51]
        lines(yplot ~ xplot, lty = w, col = scales::alpha(col_pal[j], 0.05)) 
      }
      
      # add mean lines
      ymean <- apply(sim_plot[[w]][, xplot, 51], 2, median)
      lines(ymean ~ xplot, col = col_pal[j], lty = w)
      
    }
    
  }
  
  # add a legend
  if (i == length(sims_to_plot)) {
    x <- seq(0, 1, by = 0.1)
    plot(seq(0, 1, by = 0.1) ~ x,
         bty = "n", xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", type = "n")
    legend(x = 0.5, y = -0.2, 
           ncol = 2,
           xpd = TRUE,
           legend = c("Historical: do nothing", "Post-1997: do nothing", "RCP8.5 low: do nothing", "RCP8.5 medium: do nothing", "RCP8.5 high: do nothing",
                      "Historical: all interventions", "Post-1997: all interventions", "RCP8.5 low: all interventions", "RCP8.5 medium: all interventions", "RCP8.5 high: all interventions"),
           col = rep(col_pal, times = 2), 
           lty = rep(seq_along(sim_plot), each = 5), 
           lwd = 2, xjust = 0.5, yjust = 0.5)
  }
  
}

# close plotting device
dev.off()
