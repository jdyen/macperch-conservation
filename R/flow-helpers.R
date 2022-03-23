# functions to load and format flow scenarios

# internal function to add a season ID to data.frames
add_season <- function(x, cool_season = 4:10) {
  x %>% mutate(
    season = ifelse(month(date_formatted) %in% cool_season, "cool", "warm")
  )
}

# function to calculate water temp effect on spawning based on deviations from mean value
#   Parameters from Tonkin et al. 2017 (Ecohydro)
#   rescale values from Ovens River @ Peechelba (mean and SD of water temp)
calculate_simplified_temp_metrics <- function(
  x, 
  date,
  coef = 0.243,
  rescale = c(19.79, 3.29),
  bounds = c(0, 2)
) {
  
  # calculate average temperature over the spawning period (Oct-Dec)
  average_temperature <- calculate(
    value = x,
    date = date,
    standardise = NULL,
    resolution = survey(season = 10:12)
  )$metric
  
  # work out standardised deviations 
  temperature_std <- (average_temperature - rescale[1]) / rescale[2]
  
  # multiply by coef to give proportional changes in recruitment
  effect <- exp(coef * temperature_std)  
  
  # threshold values
  effect[effect < bounds[1]] <- bounds[1]
  effect[effect > bounds[2]] <- bounds[2]

  # and return
  effect
    
}
# function to calculate and return annual flow metrics
calculate_flow_metrics <- function(x, date, rescale) {
  
  # calculate change in river height
  discharge_median <- calculate(
    value = x,
    date = date,
    standardise = NULL,
    resolution = survey(season = 7:18)
  )
  river_height_change <- discharge_median$metric[-1] / discharge_median$metric[-nrow(discharge_median)]
  
  # make this a percentage change
  river_height_change <- (river_height_change * 100) - 100
  
  # and add a leading 0 because year 1 is the base level
  river_height_change <- c(0, river_height_change)
  
  # combine everything into a data frame
  flow_covariates <- data.frame(
    
    # Average daily flow, annual
    average_daily_flow = calculate(
      value = x,
      date = date,
      standardise = NULL,
      resolution = survey(season = 7:18)
    )$metric,
    
    # Nov-Dec: standardised average daily flow
    spawning_flow = calculate(
      value = x,
      date = date,
      standardise = NULL,
      resolution = survey(season = 11:12)
    )$metric,
    
    # Nov-Dec: number of days with > 100 % change from previous
    spawning_variability = calculate(
      value = x,
      date = date,
      standardise = NULL,
      fun = flow_variability,
      resolution = survey(season = 11:12)
    )$metric,
    
    # change in water level relative to previous year
    river_height_change = river_height_change,
    
    # calculate low-flow value
    min_daily_flow = calculate(
      value = x,
      date = date,
      fun = min,
      standardise = NULL,
      resolution = survey(season = 7:18)
    )$metric

  )
  
  # check for NAs in spawning variability (ratio that can have zero denominator)
  if (anyNA(flow_covariates$spawning_variability)) {
    idx <- is.na(flow_covariates$spawning_variability)
    flow_covariates$spawning_variability[idx] <- max(flow_covariates$spawning_variability, na.rm = TRUE)
  }
      
  # rescale average daily and spawning flows
  flow_covariates$average_daily_flow <- 
    flow_covariates$average_daily_flow / rescale
  flow_covariates$spawning_flow <-
    flow_covariates$spawning_flow / rescale
  
  # add row names
  min_year <- min(lubridate::year(date))
  max_year <- max(lubridate::year(date)) - 1L
  rownames(flow_covariates) <- min_year:max_year
  
  # return
  flow_covariates
  
}

# calculate number of days with change > 100 %
flow_variability <- function(x) {
  sum(x[-1] / x[-length(x)] > 2)
}

# function to add climate change impacts to discharge sequences
add_climate_change_scenarios <- function(x, catchment, scenario, reference, type = "discharge", variable = "value") {
  
  # check catchment and scenario
  if (!catchment %in% c(
    "upper_murray", 
    "goulburn",
    "yarra", 
    "lower_murray", 
    "campaspe",
    "snowy",
    "werribee",
    "ovens",
    "wimmera"
  )) {
    stop("catchment not implemented; using default Victorian values", call. = FALSE)
  }
  
  # catchments (values for 2040 and 2065 relative to 1995)
  if (type == "discharge") {
    
    effect <- list(
      upper_murray = list(
        rcp85 = list(low = c(17.2, 13.5), medium = c(-8.4, -16.6), high = c(-23.3, -39.4)),
        rcp45 = list(low = c(14.5, 16.3), medium = c(-5.2, -5.6), high = c(-26.4, -37.4))
      ),
      lower_murray = list(
        rcp85 = list(low = c(32.8, 27.1), medium = c(-4.6, -11.4), high = c(-37.5, -47.0)),
        rcp45 = list(low = c(21.6, 29.3), medium = c(1.2, -5.4), high = c(-27.6, -39.3))
      ),
      goulburn = list(
        rcp85 = list(low = c(9.9, 1.3), medium = c(-9.5, -13.7), high = c(-29.1, -41.9)),
        rcp45 = list(low = c(12.2, 8.1), medium = c(-3.8, -11.7), high = c(-28.7, -33.1))
      ),
      campaspe = list(
        rcp85 = list(low = c(10.5, 1.0), medium = c(-12.3, -20.7), high = c(-37.1, -57.0)),
        rcp45 = list(low = c(20.8, 9.1), medium = c(-6.4, -12.3), high = c(-43.9, -43.6))
      ),
      yarra = list(
        rcp85 = list(low = c(10, 0.8), medium = c(-11, -16.4), high = c(-29.2, -44.3)),
        rcp45 = list(low = c(10.6, 8.7), medium = c(-3.1, -11.4), high = c(-30.0, -34.0))
      ),
      snowy = list(
        rcp85 = list(low = c(22.5, 21.0), medium = c(-7.1, -17.9), high = c(-25.3, -36.1)),
        rcp45 = list(low = c(13.8, 17.3), medium = c(-3.3, -9.3), high = c(-29.9, -31.8))
      ),
      werribee = list(
        rcp85 = list(low = c(11.8, 7.5), medium = c(-7.7, -18.1), high = c(-28.9, -45.5)),
        rcp45 = list(low = c(16.5, 10.7), medium = c(-3.8, -7.6), high = c(-35.0, -36.5))
      ),
      ovens = list(
        rcp85 = list(low = c(11.7, 1.2), medium = c(-10.8, -6.0), high = c(-23.3, -43.9)),
        rcp45 = list(low = c(12.8, 9.6), medium = c(-6.0, -15.9), high = c(-31.0, -34.4))
      ),
      wimmera = list(
        rcp85 = list(low = c(12.1, 12.3), medium = c(-6.5, -14.4), high = c(-32.3, -53.1)),
        rcp45 = list(low = c(21.0, 11.5), medium = c(-4.4, -12.0), high = c(-33.8, -38.6))
      ),
      default = list(
        rcp85 = list(low = c(8.7, 1.5), medium = c(-8.5, -15.9), high = c(-24.7, -43.8)),
        rcp45 = list(low = c(14.0, 9.4), medium = c(-1.6, -11.1), high = c(-29.1, -33.3))
      )
    )
    
  } else {
    
    if (type == "water_temperature") {
      
      # use air temperature multiplied by 0.65 based on linear model of observed
      #   air temperature against water temperature in the lower Murray, Goulburn,
      #   and Campaspe systems
      
      effect <- list(
        upper_murray = list(
          rcp85 = list(low = c(1.1, 1.9), medium = c(1.4, 2.6), high = c(1.7, 3.0)),
          rcp45 = list(low = c(0.8, 1.2), medium = c(1.1, 1.6), high = c(1.5, 2.1))
        ),
        lower_murray = list(
          rcp85 = list(low = c(1.1, 2.0), medium = c(1.5, 2.5), high = c(1.7, 3.1)),
          rcp45 = list(low = c(0.7, 1.2), medium = c(1.1, 1.6), high = c(1.4, 2.0))
        ),
        goulburn = list(
          rcp85 = list(low = c(1.0, 2.0), medium = c(1.4, 2.4), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.4, 1.9))
        ),
        campaspe = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.4), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.4), high = c(1.4, 1.8))
        ),
        yarra = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.8))
        ),
        snowy = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.4, 2.5), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.8, 1.2), medium = c(1.1, 1.5), high = c(1.4, 2.0))
        ),
        
        werribee = list(
          rcp85 = list(low = c(1.0, 1.8), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.4), high = c(1.2, 1.8))
        ),
        ovens = list(
          rcp85 = list(low = c(1.0, 2.0), medium = c(1.4, 2.5), high = c(1.6, 3.0)),
          rcp45 = list(low = c(0.7, 1.2), medium = c(1.1, 1.6), high = c(1.5, 2.0))
        ),
        wimmera = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.6, 2.9)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.9))
        ),
        default = list(
          rcp85 = list(low = c(1.0, 1.9), medium = c(1.3, 2.3), high = c(1.5, 2.8)),
          rcp45 = list(low = c(0.7, 1.1), medium = c(1.0, 1.5), high = c(1.3, 1.8))
        )
      )
      
      # apply 0.65 multiplier to all
      effect <- lapply(effect, function(x) lapply(x, function(y) lapply(y, function(z) z * 0.65)))
      
    } else {
      stop("only discharge and water temperature are currently implemented", call. = FALSE)
    }
  }
  
  # pull out relevant effect
  scaling_factor <- effect[[catchment]][[scenario]]
  
  # calculate scaled discharge
  if (reference > 2075)
    stop("scaling will not extend beyond 2075", call. = FALSE)
  lowz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$low[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$low[1] +
        (scaling_factor$low[2] - scaling_factor$low[1]) * (reference - 2040) / (2065 - 2040)
    )
  )
  mediumz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$medium[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$medium[1] +
        (scaling_factor$medium[2] - scaling_factor$medium[1]) * (reference - 2040) / (2065 - 2040)
    )
  )
  highz <- ifelse(
    reference < 1995, 
    0,
    ifelse(
      reference <= 2040,
      scaling_factor$high[1] * (reference - 1995) / (2040 - 1995),
      scaling_factor$high[1] +
        (scaling_factor$high[2] - scaling_factor$high[1]) * (reference - 2040) / (2065 - 2040)
    )
  )

  # rescale appropriately for different data types
  if (type == "discharge") {
    x[[paste0(variable, "_", scenario, "low")]] <-
      x[[variable]] * (1 + lowz / 100)
    x[[paste0(variable, "_", scenario, "med")]] <-
      x[[variable]] * (1 + mediumz / 100)
    x[[paste0(variable, "_", scenario, "high")]] <-
      x[[variable]] * (1 + highz / 100)
  } else {
    x[[paste0(variable, "_", scenario, "low")]] <-
      x[[variable]] + lowz
    x[[paste0(variable, "_", scenario, "med")]] <-
      x[[variable]] + mediumz
    x[[paste0(variable, "_", scenario, "high")]] <-
      x[[variable]] + highz
  }

  # and return
  x
  
}

# internal function to rescale discharge
rescale_discharge <- function(data, variable = "value") {
  
  # add a season ID to data.frames
  data <- data %>% mutate(
    season = ifelse(month(date_formatted) %in% 4:10, "cool", "warm")
  )
  
  # decile rescaling for years prior to 1975
  data[[paste0(variable, "_1975")]] <- rescale_segment(
    x = data[[variable]],
    idx = data$date_formatted > ymd("1975-06-30"),
    season = data$season
  )
  
  # repeat decile rescaling for years prior to 1997
  data[[paste0(variable, "_1997")]] <- rescale_segment(
    x = data[[variable]],
    idx = data$date_formatted > ymd("1997-06-30"),
    season = data$season
  )
  
  # remove added columns
  data <- data %>% select(-season)

  # return
  data
  
}

# function to rescale flows based on seasonal quantiles
seasonal_quantile_rescale <- function(
  x,
  reference,
  x_season = NULL,
  reference_season = NULL,
  ...
) {
  
  # only one season if NULL
  if (is.null(x_season)) {
    x_season <- rep(1, length(x))
    reference_season <- rep(1, length(reference))
  }
  
  # need the reference season
  if (is.null(reference_season))
    stop("reference_season must be provided", call. = FALSE)
  
  # use quantile_rescale in each season
  seasons <- unique(x_season)
  for (i in seq_along(seasons)) {
    idx <- x_season == seasons[i]
    idy <- reference_season == seasons[i]
    x[idx] <- quantile_rescale(
      x = x[idx], reference = reference[idy], ...
    )
  }
  
  # and return
  x
  
}

# function to rescale by quantile
quantile_rescale <- function(
  x, 
  reference, 
  probs = seq(0, 1, by = 0.1)
) {
  
  # work out quantiles in reference (e.g. future climate)
  reference_quantiles <- get_quantile(
    x = reference, probs = probs
  )
  
  # work out quantiles in observed
  x_quantiles <- get_quantile(
    x = x, probs = probs
  )
  
  # simplified change ratio: proportional change in each quantile
  change_ratio <-
    reference_quantiles$quantiles / x_quantiles$quantiles

  # apply this change to observed values  
  x <- x * change_ratio[x_quantiles$bins]
  
  # return, ensure positive
  ifelse(x < 0, 0, x)
  
}

# internal function to calculate quantiles and bins with 
#   closed intervals
get_quantile <- function(x, probs) {
  
  # calculate quantiles of x  
  breaks <- quantile(x, probs = probs)
  
  # reduce first quantile to close interval at both ends
  breaks[1] <- breaks[1] - 1e-3
  
  # bin x into its quantiles
  bins <- cut(
    x, breaks = breaks, labels = FALSE
  )
  
  # and return mean of each bin and bins
  list(
    quantiles = tapply(
      x, bins, mean
    ),
    bins = bins
  )
  
}

# function to rescale part of a sequence based on the remaining
#   years in that sequence
rescale_segment <- function(x, idx, season = NULL, ...) {

  # all one season if not specified
  if (is.null(season))
    season <- rep(1, length(x))
  
  # calculate and return
  x_segment <- seasonal_quantile_rescale(
    x = x[!idx],
    reference = x[idx],
    x_season = season[!idx],
    reference_season = season[idx] 
  )

  # update rescaled segment
  x[!idx] <- x_segment
  
  # and return
  x
  
}

# check years of data by site
check_available <- function(x) {
  years_available <- tapply(
    x$value, 
    list(x$site_code, year(x$date_formatted), x$variable_name), 
    function(x) ifelse(sum(!is.na(x)) > 300, 1, 0)
  )
  years_available[is.na(years_available)] <- 0
  years_available
}

# fill gaps in discharge by resampling
resample_discharge <- function(x, available, variable = "value") {
  
  year_set <- as.numeric(colnames(available))
  for (i in seq_len(nrow(available))) {
    
    site_id <- rownames(available)[i]
    
    target <- year_set[available[site_id, , 1] == 0]
    source <- year_set[available[site_id, , 1] == 1]
    
    if (length(target) > 0) {
      
      x[[variable]] <- 
        resample_discharge_internal(
          x[[variable]],
          x$date_formatted,
          target = target,
          source = source
        )
      
    }
    
  }
  
  # return
  x
  
}

# function to resample years of discharge, accounting for leap years
resample_discharge_internal <- function(value, date, target, source) {
  
  # create a data.frame to resample
  data <- data.frame(
    date_formatted = date,
    value = value
  )
  
  # sample years from source to replace target years,
  #   accounting for leap years
  idx <- leap_year(source)
  idy <- leap_year(target)
  filled_years <- rep(NA, length(target))
  
  # check if we need to sample with replacement (if there are too many years missing)
  replace_leap <- replace_normal <- FALSE
  if (sum(idx) < sum(idy))
    replace_leap <- TRUE
  if (sum(!idx) < sum(!idy))
    replace_normal <- TRUE
  
  # and do the actual resampling
  filled_years[idy] <- sample(source[idx], size = sum(idy), replace = replace_leap)
  filled_years[!idy] <- sample(source[!idx], size = sum(!idy), replace = replace_normal)
  
  # resample observed discharge years and collapse into a single data.frame
  resampled <- do.call(
    rbind,
    lapply(
      filled_years, 
      function(x, data) data %>% filter(year(date_formatted) == x),
      data = data
    )
  )
  
  # fix the dates on these years to match target years
  resampled <- resampled %>% 
    mutate(
      date_formatted = ymd(
        paste(
          rep(target, times = ifelse(leap_year(target), 366, 365)),
          month(date_formatted),
          day(date_formatted),
          sep = "-"
        )
      )
    )

  # combine with target years and arrange chronologically
  resampled <- resampled %>% 
    bind_rows(data %>% filter(year(date_formatted) %in% source)) %>%
    arrange(date_formatted)
  
  # and return value only
  resampled[, 2]
  
}

# fill remaining gaps with a fill_rolling_na
fill_missing <- function(x) {
  
  idx <- apply(x, 2, function(.x) sum(is.na(.x))) > 0
  to_fill <- names(idx)[idx]
  for (i in seq_along(to_fill)) {
    x[[to_fill[i]]] <- fill_na_rolling(
      x, variable = to_fill[i], recursive = TRUE, max_iter = 20
    )
  }
  
}

# function to fill missing data iteratively based on rolling means of preceding days
fill_na_rolling <- function(
  flow, variable = "value", recursive = FALSE, max_iter = 20
) {
  
  x <- flow[[variable]]
  
  if (any(is.na(x)) & !all(is.na(x))) {
    
    if (recursive) {
      
      counter <- 1
      while(any(is.na(x)) & counter < (max_iter + 1)) {
        
        n <- length(x)
        idx <- sapply(rev(seq_len(5)) - 1, function(x) rep(x, n))
        idx <- sweep(idx, 1, seq_len(n), "+")
        idx <- ifelse(idx > n, NA, idx)
        df <- matrix(x[idx], nrow = n)
        
        mean_val <- apply(df, 1, mean, na.rm = TRUE)    
        
        idx <- which(is.na(x))
        
        flow[[variable]][idx] <- mean_val[idx]
        
        x <- flow[[variable]]
        
        counter <- counter + 1
        
      }
      
    } else {
      
      n <- length(x)
      idx <- sapply(rev(seq_len(5)) - 1, function(x) rep(x, n))
      idx <- sweep(idx, 1, seq_len(n), "+")
      idx <- ifelse(idx > n, NA, idx)
      df <- matrix(x[idx], nrow = n)
      
      mean_val <- apply(df, 1, mean, na.rm = TRUE)    
      
      idx <- which(is.na(x))
      
      flow[[variable]][idx] <- mean_val[idx]
      
    }
    
  }
  
  flow
  
}

# helper functions to impute temperature
add_bom_data <- function(flow, temp, variable = "max_temperature") {
  
  temp <- temp %>% mutate(
    date_formatted = parse_date_time(paste(temp$year, temp$month, temp$day, sep = "-"), orders = c("ymd"))  
  )
  
  temp <- temp %>% select(date_formatted, !!sym(variable)) %>% as.data.frame
  
  left_join(flow, temp, by = "date_formatted")
  
}

filter_qc <- function(x, variable = "value", quality = "quality_code", threshold = 150) {
  
  idx <- x[[quality]] > threshold
  x[[variable]][idx] <- NA
  
  x
  
}

filter_varcode <- function(x, priority = "141.00") {
  
  varcodes <- unique(x$variable_code)
  
  if (length(varcodes) > 1) {

    if (priority %in% varcodes) {
      x <- x[x$variable_code == priority, ]
    } else {
      x <- x[x$variable_code == varcodes[1], ]
      warning(
        "priority variable code not in data, using first unique",
        " variable code, which is ",
        varcodes[1],
        ". Alternative variable codes for this data set: ", 
        varcodes[-1],
        call. = FALSE
      )
    }

  }
  
  x
  
}

fill_na <- function(flow, mod, response = "value") {
  
}

# function to impute temperatures based on air temperature
impute_temperature <- function(data, target, site, latitude, longitude, response = "value") {
  
  # add a month factor for the regression
  data <- data %>% mutate(
    day_of_year = yday(date_formatted),
    month = month(date_formatted),
    month_fac = factor(month),
    discharge = target
  )
  
  # regression: water temp ~ air temp + factor(month)
  idx <- is.na(data[[response]])
  if (sum(idx) > 0) {
    temp_lm <- lm(
      as.formula(paste0(response, " ~ ", "discharge + day_of_year + day_of_year^2 + month_fac")),
      data = data
    ) 
    data[[response]][idx] <- predict(temp_lm, newdata = data[idx, ])
  }
  
  # if still missing, use previous 5 days recursively
  data <- fill_na_rolling(data, recursive = TRUE)
  
  # remove added columns
  data <- data %>% select(
    -day_of_year, -month, -month_fac, - discharge
  )
  
  # return
  data
  
}


# function to revise flows according to environmental watering plan
#   for the Yarra River
add_environmental_water <- function(x, date, system = "yarra") {
  
  # id each month and year
  month_id <- month(date)
  year_id <- year(date)
  
  # loop over each year (skip first because working with
  #   water years and will lag below)
  unique_years <- unique(year_id)[-1]
  for (i in seq_along(unique_years)) {
    
    # subset to correct year/month combos (shifted water year, starting in June)
    idx <- (year_id == (unique_years[i] - 1) & month_id %in% c(6:12)) |
      (year_id == unique_years[i] & month_id %in% c(1:5))
    
    # update with env watering plans
    x[idx] <- add_environmental_water_annual(
      x = x[idx],
      month_id = month_id[idx],
      system = system
    )
    
  }
  
  # return
  x
  
}

# internal function to add env watering plan for a single year
#   (shifted water years: June-May)
add_environmental_water_annual <- function(x, month_id, system) {
  
  # set parameters based on SWPs for each system
  if (!system %in% c("yarra", "goulburn"))
    stop("environmental water is implemented for the Yarra and Goulburn rivers only", call. = FALSE)
  if (system == "yarra") {
    pars <- list(
      baseflow_one = list(month = c(12, 1:5), range = c(80, 200)),
      baseflow_two = list(month = c(6:11), range = c(200, 350)),
      fresh_one = list(month = c(12, 1:3), range = c(350, 750), n = 3, duration = c(2, 4)),
      fresh_two = list(month = c(4:5), range = c(560, 1300), n = 1, duration = c(7, 14)),
      fresh_three = list(month = c(6:11), range = c(700, 2500), n = 2, duration = c(3, 7)),
      fresh_four = list(month = c(9), range = c(700, 2500), n = 1, duration = c(14, 14))
    )
  } else {
    pars <- list(
      baseflow_one = list(month = c(1:12), range = c(500, 830)),
      baseflow_two = list(month = c(4:9), range = c(350, 450)),
      fresh_one = list(month = c(7:10), range = c(6600, 7500), n = 2, duration = c(14, 14)),
      fresh_two = list(month = c(6:7), range = c(6600, 15000), n = 1, duration = c(14, 14)),
      fresh_three = list(month = c(9:12, 1:2), range = c(6600, 7500), n = 1, duration = c(1, 1)),
      fresh_four = NULL
    )
  }
  
  # summer/autumn baseflows
  #   maintain 80-200 ML/day from Dec-May
  x <- sample_baseflow(
    x = x, 
    month_id = month_id, 
    month_subset = pars$baseflow_one$month,
    range = pars$baseflow_one$range
  )
  
  #   maintain 200-350 ML/day from June-Nov
  x <- sample_baseflow(
    x = x, 
    month_id = month_id, 
    month_subset = pars$baseflow_two$month,
    range = pars$baseflow_two$range
  )
  
  #   1-3 freshes of 350-750 ML/day for 2-4 days during Dec-May
  x <- sample_fresh(
    x = x,
    month_id = month_id,
    month_subset = pars$fresh_one$month,
    n = pars$fresh_one$n,
    range = pars$fresh_one$range,
    duration = pars$fresh_one$duration
  )
  
  #   fresh of 560-1300 ML/day for 7-14 days during April-May
  x <- sample_fresh(
    x = x,
    month_id = month_id,
    month_subset = pars$fresh_two$month,
    n = pars$fresh_two$n,
    range = pars$fresh_two$range,
    duration = pars$fresh_two$duration
  )
  
  #   1-2 freshes of 700-2500 ML/day for 3-7 days during June-Nov
  x <- sample_fresh(
    x = x,
    month_id = month_id,
    month_subset = pars$fresh_three$month,
    n = pars$fresh_three$n,
    range = pars$fresh_three$range,
    duration = pars$fresh_three$duration
  )
  
  #   fresh of 700-2500 ML/day for 14 days in Sept
  if (!is.null(pars$fresh_four)) {
    x <- sample_fresh(
      x = x,
      month_id = month_id,
      month_subset = pars$fresh_four$month,
      n = pars$fresh_four$n,
      range = pars$fresh_four$range,
      duration = pars$fresh_four$duration
    )
  }
  
  # return
  x
  
}

# function to simulate baseflows with a set range in particular
#   months
sample_baseflow <- function(x, month_id, month_subset, range) {
  
  # subset to correct months
  idx <- month_id %in% month_subset
  
  # create a sequence of flow sizes to work with below
  range_seq <- range[1]:range[2]
  
  # fill with sampled baseflow within range if flow is below
  #   minimum baseflow level
  x[idx] <- ifelse(
    x[idx] < range[1], 
    sample(range_seq,
           size = sum(x[idx] < range[1]),
           replace = TRUE,
           prob = (1 / (range_seq - range[1] + 1))),
    x[idx]
  )
  
  # and return
  x
  
}

# function to simulate a specified number of freshes of set range
#   and duration, checking whether these events occur naturally
#   before adding
sample_fresh <- function(x, month_id, month_subset, n, range, duration) {
  
  # subset to correct months
  idx <- month_id %in% month_subset
  
  # check to see if there is already a fresh of the required
  #   duration and range (n_obs > 1)
  is_fresh <- x[idx] >= range[1]
  fresh_length <- rle(is_fresh)
  fresh_length <- fresh_length$lengths[fresh_length$values]
  n_obs <- sum(fresh_length >= duration[1])
  
  # if not, add between 1 and n freshes (minus observed freshes)
  if (n_obs == 0) {
    
    n_fresh <- sample(1:n, size = 1)
    
    # for each added fresh, sample the size
    fresh_size <- runif(n = n_fresh, min = range[1], max = range[2])
    
    # and sample the duration
    fresh_duration <- resample(
      duration[1]:duration[2], size = n_fresh, replace = TRUE
    )
    
    # sample time period randomly
    already_fresh <- NULL
    for (i in seq_len(n_fresh)) {
      available_days <- sum(idx) - fresh_duration[i] + 1
      start <- resample(seq_len(available_days), size = 1)
      end <- start + fresh_duration[i] - 1
      ntry <- 0
      while (any(start:end) %in% already_fresh & ntry < 100) {
        start <- resample(seq_along(x[idx]), size = 1)
        end <- start + fresh_duration[i] - 1
        ntry <- ntry + 1
      }
      xtmp <- x[idx][start:end]
      xtmp <- ifelse(
        xtmp < fresh_size[i],
        fresh_size[i],
        xtmp
      )
      
      x[idx][start:end] <- xtmp
      already_fresh <- c(already_fresh, start:end)
    }
    
  }
  
  # return
  x
  
}

# helper function to avoid issue in sample when vector is a single
#   integer value
resample <- function(x, ...) x[sample.int(length(x), ...)]

# calculate e-water allocations for a single site
calculate_ewater_contribution <- function(x, y) {
  out <- x %>% select(contains("value_")) - y %>% select(contains("value_"))
  apply(out, 2, sum)
}

# function to define effects of discharge on recruitment and survival (in rivers)
extract_covar_effects <- function(
  mat, 
  x, 
  spawning_param = c(-0.01, -0.05), 
  variability_param = -0.003,
  recruit_param = -0.01,
  shift = 200, 
  survival_param = c(0.2, -0.2),
  ctf_param = 1,
  ctf_threshold = 3,
  ...
) {
  
  log_flow <- log(x$spawning_flow + 0.01)
  spawning_flow_effect <- exp(spawning_param[1] * log_flow + spawning_param[2] * 
                                (log_flow^2))
  spawning_flow_effect[spawning_flow_effect > 1] <- 1
  spawning_flow_effect[spawning_flow_effect < 0] <- 0
  spawning_var_effect <- exp(-variability_param * x$spawning_variability)
  
  temp_effect <- x$temperature_effect
  
  height_effect <- 0.4 + (1/(1 + exp(recruit_param * (x$river_height_change + shift))))
  
  daily_log_flow <- log(x$average_daily_flow + 0.01)
  daily_flow_effect <- exp(survival_param[1] * daily_log_flow + 
                             survival_param[2] * (daily_log_flow^2))
  daily_flow_effect[daily_flow_effect > 1] <- 1
  daily_flow_effect[daily_flow_effect < 0] <- 0
  
  ctf_effect <- ifelse(
    x$min_daily_flow < ctf_threshold,
    exp(ctf_param * x$min_daily_flow) / exp(ctf_threshold * ctf_param),
    1
  )
  
  # return
  list(
    average_daily_flow = daily_flow_effect,
    spawning_flow = spawning_flow_effect,
    spawning_variability = spawning_var_effect,
    river_height_change = height_effect,
    min_daily_flow = ctf_effect,
    temperature_effect = temp_effect
  )
  
}
