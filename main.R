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
  goulburn = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration", "env_water"),
  king = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration"),
  ovens = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration"),
  mitta = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration", "exclude_exotic"),
  sevens = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration"),
  yarra = c("gene_mixing", "stocking5", "stocking10", "fishing_regulations", "habitat_restoration", "env_water")
)

# expand management actions
actions <- lapply(actions, expand_combn)

# don't want to keep multiple stocking actions in a single run
actions <- lapply(actions, remove_conflicts, conflict = c("stocking5", "stocking10"))

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
# Yarra River @ Millgrove 229212 (broken, using downloaded Warrandyte data in "data/")
#    Latitude: 37°45'09.3"S, Longitude: 145°39'19.7"E
#    Currently ignoring water temperature data (back-filled with Ovens data for convenience)
#    Loading discharge data from saved file due to errors in ratings table in Yarra catchment
discharge_gauges <- c("goulburn" = "405201", "king" = "403228", "ovens" = "403200", "mitta" = "401203", "sevens" = "405234", "yarra" = "405201")
water_temp_gauges <- c("goulburn" = "405201", "king" = "403241", "ovens" = "403241", "mitta" = "401203", "sevens" = "405269", "yarra" = "403241")
site_names <- c("Trawool", "LWHovell", "Peechelba", "Hinnomunjie", "Kialla", "Millgrove")
site_codes <- c(405201, 403241, 403241, 401203, 405269, 403241)
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
  yarra_discharge <- read.csv("data/Warrandyte-Daily-River-Flow.csv")
  yarra_discharge$date_formatted <- parse_date_time(
    yarra_discharge$Date, orders = c("ymd_HM")
  )
  colnames(yarra_discharge)[4] <- "mean_discharge" 
  discharge$yarra <- discharge$yarra %>%
    left_join(yarra_discharge %>% select("date_formatted", "mean_discharge"),
              by = "date_formatted")
  discharge$yarra <- discharge$yarra %>% mutate(
    value = mean_discharge
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
  
  # some NaN values in Yarra discharge towards end (large chunk missing)
  # Fill with median over preceding two weeks
  discharge$yarra$value[is.nan(discharge$yarra$value)] <- 
    median(
      discharge$yarra$value[year(discharge$yarra$date_formatted) == 2020 &
                              month(discharge$yarra$date_formatted) == 12],
      na.rm = TRUE
    )
  
  # fill temperature based on nearby air temperatures, with final gaps filled
  #   with rolling means
  water_temp <- mapply(
    impute_temperature_temp,
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
  "goulburn" = 330,
  "king" = 330,
  "ovens" = 330,
  "mitta" = 307,
  "sevens" = 33,
  "yarra" = 344
)


## TODO: update habitat effect and fix exotic predators

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

# set covariate effects parameters
covar_parameters <- list(
  shift = 150, 
  survival_param = c(0.1, -0.1),
  recruit_param = -0.01,
  spawning_param = c(-0.1, -0.05),
  variability_param = 0.1
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
    stocking = 20000 * 0.5 * 0.13,    # account for fingerling survival and sex ratio
    fishing = 0.1,
    broodfish = 15,
    exotic = 0.75
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
  start = 11
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

# prepare plots of trajectories
traj_sites <- c(
  "Goulburn River at Trawool", 
  "King River below Lake William Hovell",
  "Ovens River", 
  "Mitta Mitta River and Lake Dartmouth", 
  "Seven Creeks below Polly McQuinn Weir",
  "Yarra River at Warrandyte"
)
run_lengths <- rle(actions$site)
pinch_points <- run_lengths$lengths
site_id <- run_lengths$values
start_index <- c(1, cumsum(pinch_points)[-length(pinch_points)] + 1)
end_index <- cumsum(pinch_points)
sims_to_plot <- list(
  goulburn = lapply(1:5, function(x) c(start_index[site_id == "goulburn"][x], end_index[site_id == "goulburn"][x])),
  king = lapply(1:5, function(x) c(start_index[site_id == "king"][x], end_index[site_id == "king"][x])),
  ovens = lapply(1:5, function(x) c(start_index[site_id == "ovens"][x], end_index[site_id == "ovens"][x])),
  mitta = lapply(1:5, function(x) c(start_index[site_id == "mitta"][x], end_index[site_id == "mitta"][x])),
  sevens = lapply(1:5, function(x) c(start_index[site_id == "sevens"][x], end_index[site_id == "sevens"][x])),
  yarra = lapply(1:5, function(x) c(start_index[site_id == "yarra"][x], end_index[site_id == "yarra"][x]))
)
ylim_traj <- c(1000, 500, 1000, 1200, 400, 1000)
png(file = paste0("outputs/figs/trajectories.png"),
     height = 13.2,
     width = 8.2,
     res = 600,
     units = "in",
     pointsize = 12)
layout(mat = seq_len(length(sims_to_plot) + 1), heights = c(rep(1, length(sims_to_plot)), 0.4))
par(mar = c(4.1, 4.8, 2.1, 1.1))
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
  mtext(traj_sites[i], side = 3, line = 0, adj = 1)
  
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
    plot(seq(0, 10, by = 1) ~ x,
         bty = "n", xlab = "", ylab = "",
         xaxt = "n", yaxt = "n", type = "n")
    legend(x = 0.5, y = 7.5, 
           ncol = 2,
           xpd = TRUE,
           legend = c("Historical: no actions", "Post-1997: no actions", "RCP8.5 low: no actions", "RCP8.5 medium: no actions", "RCP8.5 high: no actions",
                      "Historical: all actions", "Post-1997: all actions", "RCP8.5 low: all actions", "RCP8.5 medium: all actions", "RCP8.5 high: all actions"),
           col = rep(col_pal, times = 2), 
           lty = rep(seq_along(sim_plot), each = 5), 
           lwd = 2, xjust = 0.5, yjust = 0.8)
  }
  
}

# close plotting device
dev.off()

# prepare plots of risk curves
sites_to_plot <- c("goulburn", "king", "ovens", "mitta", "sevens", "yarra")
action_label <- c(
  "none" = "None",
  "gene_mixing" = "Genetic mixing", 
  "stocking" = "Stocking", 
  "stocking5" = "Stocking (5 years)", 
  "stocking10" = "Stocking (10 years)",
  "fishing_regulations" = "Fishing regulations",
  "habitat_restoration" = "Instream habitat",
  "env_water" = "Environmental water",
  "exclude_exotic" = "Exclude exotic species"
)
climate_label <- c("Post-1975", "Post-1997", "RCP8.5 low", "RCP8.5 medium", "RCP8.5 high")
max_thresh <- c(700, 120, 600, 700, 80, 600)
threshold <- lapply(seq_along(unique(actions$site)), function(x) seq_len(max_thresh[x]))
col_pal <- RColorBrewer::brewer.pal(8, "RdYlGn")
lty_set <- rep(1, 8)

# and plot it
for (i in seq_along(sites_to_plot)) {
  
  # initialise plot
  png(
    file = paste0("outputs/figs/risk_curve_", sites_to_plot[i], ".png"),
    units = "in",
    width = 6,
    height = 9,
    pointsize = 15,
    res = 300
  )
  layout(mat = cbind(c(1, 3, 5), c(2, 4, 6)))
  par(mar = c(4.1, 4.8, 2.1, 1.1))
  
  emps_sub <- NULL
  for (j in seq_along(climates)) {
    
    # subset to relevant risk curves
    idx <- actions$site == sites_to_plot[i] &
      actions$climate == climates[j]
    risk_sub <- risk_curves[idx, threshold[[i]]]
    
    # and actions
    actions_sub <- actions$actions[idx, ]
    
    # remove rows with continuous stocking
    idy <- 1:6
    if (nrow(actions_sub) == 48)
      idy <- 1:7
    idy <- c(idy, nrow(actions_sub))
    risk_sub <- risk_sub[idy, ]
    actions_sub <- actions_sub[idy, ]
    
    # add EMPS
    emps_sub <- cbind(emps_sub, emps[idx][idy])

    # set up plot
    plot(
      risk_sub[1, ] ~ threshold[[i]], 
      ylim = c(0, 1), 
      las = 1,
      type = "n", 
      bty = "l",
      xlim = c(0, max_thresh[i]),
      xlab = "Threshold population size (TPS)",
      ylab = "Pr(N < TPS)"
    )
    
    # add each curve
    for (k in seq_len(nrow(risk_sub))) {
      lines(risk_sub[k, ] ~ threshold[[i]], col = col_pal[k], lwd = 2, lty = lty_set[k])
    }
    
    # add climate label
    mtext(climate_label[j], side = 3, line = 0.5, adj = 1)
    
    # add legend
    if (j == 5) {
      par(mar = rep(0, 4))
      plot(seq(0, 1, length = 10) ~ seq(0, 200, length = 10),
           type = "n",
           xaxt = "n",
           yaxt = "n",
           xlab = "",
           ylab = "",
           bty = "n")
      actions_text <- apply(actions_sub, c(1, 2), function(x) action_label[x])
      actions_text <- apply(actions_text, 1, paste0, collapse = "/")
      actions_text <- gsub("/None", "", actions_text)
      actions_text <- gsub("/", " / ", actions_text)
      actions_text[nchar(actions_text) > 50] <- "All"
      legend_text <- actions_text
      legend(
        x = 100,
        y = 0.5, 
        xjust = 0.5,
        yjust = 0.5,
        ncol = 1,
        legend = legend_text,
        lty = lty_set,
        lwd = 2, 
        col = col_pal, 
        bty = "n", 
        cex = 0.89,
        y.intersp = 0.75,
        x.intersp = 1,
        xpd = TRUE
      )
    }
    
  }
  
  # shut down plot
  dev.off()
  
}

# extract pr_persist from risk curves, based on 25 % of K by system
pr_persist <- mapply(
  function(x, y) 1 - x$risk[1, y],
  x = sim_sum, 
  y = 100
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

# calculate species-level outcomes by considering all possible
#   combinations of site-level actions
costs <- c(
  "none" = 0,
  "gene_mixing" = 2000, 
  "stocking10" = 1000,
  "stocking5" = 500,
  "fishing_regulations" = 100,
  "habitat_restoration" = 5000,
  "env_water" = 10000
)

# calculate benefit
effective_popsize <- c(
  "goulburn" = 330,
  "king" = 330,
  "ovens" = 330,
  "mitta" = 307,
  "sevens" = 33,
  "yarra" = 344
)
species_benefit <- list(
  post1975 = calculate_species_benefit(
    pop_persist, 
    action_cost = costs,
    climate = "1975",
    neff = effective_popsize,
    max_action = 5
  ),
  post1997 = calculate_species_benefit(
    pop_persist, 
    action_cost = costs,
    climate = "1997",
    neff = effective_popsize,
    max_action = 5
  ),
  rcp85low = calculate_species_benefit(
    pop_persist, 
    action_cost = costs,
    climate = "rcp85low",
    neff = effective_popsize,
    max_action = 5
  ),
  rcp85med = calculate_species_benefit(
    pop_persist, 
    action_cost = costs,
    climate = "rcp85med",
    neff = effective_popsize,
    max_action = 5
  ),
  rcp85high = calculate_species_benefit(
    pop_persist, 
    action_cost = costs,
    climate = "rcp85high",
    neff = effective_popsize,
    max_action = 5
  )
)

# calculate maximum effective pop size
max_neff <- ceiling(max(sapply(species_benefit, function(x) max(x$benefit$weighted_neff))))

# population level tables
reorder_benefit <- function(x) {
  out <- x[order(x$climate, x$persist, x$emps, decreasing = c(FALSE, TRUE, TRUE)), ]
  out$emps <- round(out$emps)
  out
}
write.csv(
  do.call(
    rbind, lapply(
      pop_persist,
      reorder_benefit
    )
  ),
  file = "outputs/tables/population_outcomes_ordered.csv"
)

# species level tables
write.csv(
  do.call(
    rbind,
    lapply(species_benefit, function(x) cbind(round(x$benefit, 3), x$actions))
  ),
  file = "outputs/tables/species_outcomes_ordered.csv"
)

# plot the ordered actions
col_pal <- RColorBrewer::brewer.pal(5, "Set1")

# initialise plot
png(
  file = paste0("outputs/figs/species_benefit.png"),
  units = "in",
  width = 7.2,
  height = 9.2,
  pointsize = 15,
  res = 300
)
layout(mat = c(1:6), heights = c(rep(1, 5), 0.2))
par(mar = c(4.1, 4.8, 2.1, 4.5))
for (i in seq_along(species_benefit)) {
  
  benefit_sub <- species_benefit[[i]]$benefit
  xseq <- seq_len(nrow(benefit_sub))
  
  # set up plot
  plot(
    benefit_sub[, 1] ~ xseq, 
    ylim = c(0, 1), 
    las = 1,
    type = "n", 
    bty = "u",
    xlab = "Rank",
    ylab = "Expected benefit"
  )
  
  # add each curve
  for (k in seq_len(ncol(benefit_sub))) {
    scaled_val <- benefit_sub[, k]
    if (k == 5)
      scaled_val <- scaled_val / max_neff
    lines(scaled_val ~ xseq, col = col_pal[k], lwd = 2)
  }
  
  # add climate label
  mtext(climate_label[i], side = 3, line = 0.5, adj = 1)
  
  # add right-hand y-axis
  ne_seq <- seq(0, max_neff, by = 300)
  axis(4, at = seq(0, 1, length = length(ne_seq)), labels = ne_seq, las = 1)
  
  # add legend
  if (i == 5) {
    par(mar = rep(0, 4))
    plot(seq(0, 1, length = 10) ~ seq(0, 16, length = 10),
         type = "n",
         xaxt = "n",
         yaxt = "n",
         xlab = "",
         ylab = "",
         bty = "n")
    legend(
      x = 8,
      y = 0.5, 
      xjust = 0.5,
      yjust = 0.5,
      ncol = 2,
      legend = c("Average (geometric)", "Average", "Pr(any)", "Pr(all)", "Weighted Ne"),
      lty = lty_set,
      lwd = 2, 
      col = col_pal, 
      bty = "n", 
      cex = 0.89,
      y.intersp = 0.75,
      x.intersp = 1,
      xpd = TRUE
    )
  }
  
}

# shut down plot
dev.off()

# plot covariate effects
# calculate covariates in each system and climate
covars <- vector("list", length = length(sites))
names(covars) <- sites
for (i in seq_along(sites)) {
  
  # initialise output vector
  covars[[i]] <- vector("list", length = length(climates))
  names(covars[[i]]) <- climates
  
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
      stocking = 20000 * 0.5 * 0.13,    # account for fingerling survival and sex ratio
      fishing = 0.1,
      broodfish = 15,
      exotic = 0.75
    )
    
    # save covariates
    covars[[i]][[j]] <- sim_settings$covars
    
  }
  
}

# calculate range of covariates over all climates
min_val <- t(sapply(covars, function(x) apply(sapply(x, function(y) apply(y, 2, min)), 1, min)))
max_val <- t(sapply(covars, function(x) apply(sapply(x, function(y) apply(y, 2, max)), 1, max)))

# and plot for each river
nplot <- 100
var_names <- colnames(min_val)

# initialise plotting device
png(
  file = paste0("outputs/figs/discharge-effects.png"),
  width = 5,
  height = 7.5,
  units = "in",
  res = 300,
  pointsize = 14
)
layout(mat = cbind(c(1, 3, 5), c(2, 4, 6)))
par(mar = c(4.1, 4.8, 2.1, 1.1))
col_pal <- RColorBrewer::brewer.pal(6, "RdYlGn")

# pull out formatted variable names
formatted_var_name <- gsub("_", " ", var_names)
formatted_var_name <- gsub(" effect", "", formatted_var_name)
formatted_var_name <- paste0(toupper(substr(formatted_var_name, 1, 1)), substr(formatted_var_name, 2, nchar(formatted_var_name)))

# work through each site and variable
ymin <- c(0.4, 0.6, 0.4, 0.6, 0.6)
ymax <- c(1.05, 1.05, 1.05, 1.05, 1.55)
yjitter <- seq(-0.02, 0.02, length = length(sites))
for (j in seq_len(ncol(min_val))) {
  
  # generic sequence of covariate j to set up plot window
  covar_range <-  seq(min(min_val[, j]), max(max_val[, j]), length = nplot)
  
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
  
  for (i in seq_along(sites)) {
    
    # fill empty data.frame with sequence for a single variable
    covar_set <-  seq(min_val[i, j], max_val[i, j], length = nplot)
    mat <- 1
    x <- data.frame(
      average_daily_flow = rep(0, nplot),
      spawning_flow = rep(0, nplot),
      spawning_variability = rep(0, nplot),
      river_height_change = rep(0, nplot),
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
      survival_param = covar_parameters$survival_param
    )
    effect <- effect[[var_names[j]]] + yjitter[i]
    
    # add each curve
    lines(effect ~ covar_set, col = col_pal[i], lwd = 2)
    
    
  }
  
  if (j == 5) {
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
      ncol = 1,
      legend = legend_text,
      lty = lty_set,
      lwd = 2, 
      col = col_pal, 
      bty = "n", 
      cex = 0.89,
      y.intersp = 0.75,
      x.intersp = 1,
      xpd = TRUE
    )
  }
  
}

dev.off()
