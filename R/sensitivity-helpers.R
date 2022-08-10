# function to plot EMPS and risk curves under different parameter sensitivity tests
plot_sensitivity <- function(
    emps,
    risk,
    sens,
    par,
    emps_file = NULL, 
    risk_file = NULL, 
    group = "site", 
    panel_names = NULL,
    secondary_group = NULL
) {
  
  # work out a formatted parameter name
  par_formatted <- gsub("_", " ", par)
  par_formatted <- paste0(
    toupper(substr(par_formatted, 1, 1)),
    tolower(substr(par_formatted, 2, nchar(par_formatted)))
  )
  
  # and add default panel names if needed
  if (is.null(panel_names)) {
    panel_names <- c(
      "mitta" = "Lake Dartmouth",
      "yarra" = "Yarra River",
      "sevens" = "Seven Creeks",
      "ovens" = "Ovens River",
      "king" = "King River",
      "goulburn" = "Goulburn River"
    )
  }
  
  # set a discrete colour palette for emps plot
  sens$colour <- as.integer(as.factor(sens[[par]]))
  sens$group_col <- as.integer(as.factor(sens[[group]]))
  col_pal <- RColorBrewer::brewer.pal(length(unique(sens[[group]])), "Set2")

  # set plotting character based on secondary group if required
  pch_set <- 16
  if (!is.null(secondary_group))
    pch_set <- as.integer(as.factor(secondary_group)) + 14
  
  # plot emps
  if (!is.null(emps_file)) {
    png(file = emps_file,
        height = 6,
        width = 6 * 1.3,
        res = 600,
        units = "in",
        pointsize = 12)
  }
  
  # set up a plot layout
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.3))
  
  # plot it
  plot(
    emps ~ sens[[par]],
    col = col_pal[sens$group_col], 
    pch = pch_set,
    bty = "l",
    las = 1,
    xlab = par_formatted,
    ylab = "EMPS"
  )
  
  # add a legend
  plot(c(0, 1), c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  pch_legend <- 16
  nrep <- 1
  col_set <- col_pal
  legend_set <- panel_names[levels(factor(sens[[group]]))]
  title = "Location"
  if (!is.null(secondary_group)) {
    pch_legend <- rep(unique(pch_set), each = length(unique(sens[[group]])))
    legend_set <- paste0(
      legend_set,
      rep(c(" (no excl.)", " (excl.)"), each = length(unique(sens[[group]])))
    )
    col_set <- rep(col_set, times = 2)
    title <- "Effectiveness"
  }
  legend(
    x = -1.5,
    y = 0.5,
    title = title,
    legend = legend_set,
    pch = pch_legend,
    yjust = 0.5,
    col = col_set,
    xpd = TRUE
  )
  
  # close plotting device if open
  if (!is.null(emps_file))
    dev.off()
  
  # set a colour palette for risk curves
  col_pal <- viridisLite::viridis(length(unique(sens$colour)))
  
  # create a file for the risk curve plot
  if (!is.null(risk_file)) {
    png(file = risk_file,
        height = 8,
        width = 2 * (8 / 3),
        res = 600,
        units = "in",
        pointsize = 12)
  }
  
  # grab old plot margins
  old_mar <- par()$mar
  
  # set plot margins
  par(mar = c(4.5, 4.7, 2.5, 0.5))
  
  # set up a plot layout
  layout(matrix(c(1, 3, 5, 7, 2, 4, 6, 7), ncol = 2), heights = c(rep(1, 3), 0.5))
  
  # response variable is a sequence of probabilities from 0 to 1 with
  #   0.01 between each
  yset <- seq(0, 1, length = 101)
  
  # set line type based on secondary group if required
  lty_set <- 1
  if (!is.null(secondary_group))
    lty_set <- abs(as.integer(as.factor(secondary_group)) - 3)
  
  # plot each group separately
  for (i in seq_along(panel_names)) {
    
    # pull out group i
    idx <- sens[[group]] == names(panel_names)[i]
    risk_sub <- risk[idx, ]
    
    # set default values if all are zero (line at y = 1)
    yset_sub <- yset
    if (all(risk_sub == 0)) {
      yset_sub <- rep(1, length(yset_sub))
      risk_sub <- sweep(risk_sub, 2, seq(0, 100, length = ncol(risk_sub)), "+")
    }
    
    # initialise plot
    plot(
      yset_sub ~ risk_sub[1, ], 
      type = "n", 
      xlim = range(c(risk_sub)),
      ylim = c(0, 1),
      bty = "l", 
      las = 1,
      xlab = "TPS",
      ylab = "Pr(N < TPS)"
    )
    
    # add lines for each value
    for (j in seq_len(nrow(risk_sub))) {
      lines(yset_sub ~ risk_sub[j, ], col = col_pal[sens$colour[j]], lty = lty_set)
    }
    
    # and add a formatted label for each panel
    mtext(panel_names[i], side = 3, adj = 1, line = 0.5)
    
  }
  
  # add a legend
  par(mar = c(5.1, 1.1, 0.5, 1.1))
  add_colourbar(col_pal, min = min(sens[[par]]), max = max(sens[[par]]), xlab = par_formatted)
  
  # reset plot settings
  par(mar = old_mar)
  
  # close plotting device if opened
  if (!is.null(risk_file))
    dev.off()
  
  # return nothing
  out <- NULL
  
}

# function to add colorbar to plots
#  adapted from: https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
add_colourbar <- function(
  palette,
  min,
  max,
  xlab = "",
  ...
) {

  # how many colours are we plotting?  
  ncolour <- length(palette) - 1L
  
  # and where are the tick marks?
  ticks <- seq(min, max, length = 11)
  
  # what's the span of each colour?
  scale <- ncolour / (max - min)
  
  # set up an empty plotting area
  plot(c(min, max), c(0, 1), type = "n", bty = "n", xaxt = "n", xlab = xlab, yaxt = "n", ylab = "")
  
  # add a y axis
  axis(1, ticks, las = 1)
  
  # work out and plot the y values with correct colour
  xvals <- min + (seq_len(ncolour) - 1) / scale 
  for (i in seq_len(ncolour))
    rect(xvals[i], 0, xvals[i] + 1 / scale, 0.3, col = palette[i], border = NA)
  
  # and return nothing
  out <- NULL
  
}

# function to calculate the frequency of intervention inclusion in each site/climate
get_frequency <- function(x, var = "emps_mid") {
  
  # pull out relevant rows
  idy <- x[[var]] >= (0.98 * max(x[[var]]))
  
  # create a table of frequency of inclusions in top 2% of models
  freq_table <- table(unlist(x[idy, grepl("actions", colnames(x))])) / sum(idy)
  
  # add population establishment for new sites
  if (sites_to_plot[i] %in% c("king", "goulburn")) {
    freq_table <- c(1.0, freq_table)
    names(freq_table)[1] <- "pop_establish"
  }
  
  # add total count of top-performing models
  freq_table <- c(sum(idy), freq_table)
  names(freq_table)[1] <- "n"
  
  # return
  freq_table[names(freq_table) != "none"]
  
}

# function to plot the frequency of intervention inclusion under each climate/scenario
plot_frequency <- function(x, sites, actions) {
  
  # set layout
  laymat <- matrix(c(seq_along(sites), rep(length(sites) + 1, 3)), byrow = TRUE, ncol = 3)
  layout(laymat, heights = c(rep(1, 2), 0.35))
  
  # define legend text
  legend_text <- c(
    "mitta" = "Lake Dartmouth", 
    "yarra" = "Yarra River",
    "sevens" = "Seven Creeks",
    "ovens" = "Ovens River", 
    "king" = "King River",
    "goulburn" = "Goulburn River"
  )
  action_label <- c(
    "pop_establish" = "Establish population",
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
  
  # set plot margins
  old_mar <- par()$mar
  par(mar = c(6.1, 4.1, 2.1, 0.5))
  
  # define colour palette
  col_pal <- RColorBrewer::brewer.pal(length(actions), "Set3")
  
  # plot barplots for each site
  for (i in seq_along(sites)) {
    
    # sort action_freq
    act_tmp <- x[[i]][, order(colnames(x[[i]]))]
    
    # remove dud columns
    idx <- !colnames(act_tmp) %in% c("none", "n", "null")
    act_tmp <- act_tmp[, idx]
    act_tmp <- act_tmp[, match(actions, colnames(act_tmp))]
    
    # match actions and colours
    idy <- match(colnames(act_tmp), actions)  
    
    # plot it
    rownames(act_tmp) <- c("Historical", "Drying", "Variable")
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
    mtext(legend_text[sites[i]], side = 3, adj = 0, line = 0.4, cex = 1)
    
  }
  
  # set plot margins for legend panel and add legend
  par(mar = rep(0, 4))
  plot(c(0, 1) ~ 1, bty = "n", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  legend(x = "center", fill = col_pal, legend = action_label[actions], cex = 1.2, bty = "n", horiz = FALSE, ncol = 3)
  
  # return nothing
  out <- NULL
  
}


