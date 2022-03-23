## Using population models to inform conservation prioritisation across multiple populations

This repository contains code and data to support the following manuscript:
Yen JDL et al. (in prep.) Assisted gene flow and population reintroductions critical to long-term persistence of an endangered freshwater fish.

Last updated: 23/03/2022


### Abstract

Conservation prioritisation is increasingly focused on outcomes for species rather than individual populations. Decisions to prioritise some populations, potentially at the expense of others, require reliable assessments of species-level outcomes under different management strategies and future threat landscapes, both of which can differ markedly among populations. Here, we develop an approach to prioritise management interventions across multiple populations in highly variable and uncertain environments. We used a stochastic model of population dynamics to assess the effectiveness of seven management interventions targeting six populations of the endangered Australian freshwater fish, the Macquarie perch (Macquaria australasica). All management interventions were predicted to improve Macquarie perch population outcomes (increased abundances, reduced risk of quasi-extinction) to some extent, with the risk of declines lowest when management included assisted gene flow and fishing regulations. The lowest levels of quasi-extinction risk occurred in locations with large populations and stable environmental conditions, with outcomes for the species (i.e. aggregated over populations) similar when management interventions targeted these locations at the expense of those with smaller populations and highly variable environmental conditions. This result suggests that reintroducing Macquarie perch to large, regulated rivers with access to suitable spawning habitat (shallow, fast-flowing riffles with gravel substrate) may be more beneficial for the species than continued investment in extant but highly threatened populations in small, unregulated streams. Our study demonstrates that population models can inform spatial conservation planning without sacrificing important information on demography, genetics, or spatial and temporal environmental variation.


### Usage

The `main.R` script runs the entire analysis and generates all outputs. Note that this script may take several hours to run. Several helper functions are provided in scripts in the `R` directory; these are sourced directly within `main.R`.

Population models are simulated with the `aae.pop` R package, with the template for the Macquarie perch model included in the `aae.pop.templates` R package. Flow data are downloaded and compiled with the `aae.hydro` R package. These three packages are not available on CRAN but can be installed directly from GitHub with the remotes R package (e.g. `remotes::install_github("aae-stats/aae.pop")`). 

Discharge and water temperature data are provided to avoid downloading and re-calculating flow metrics for each simulation. Downloads from the Victorian WMIS are unreliable for some gauges (the Yarra River, in particular).

Several additional directories will need to be created to store model outputs: `outputs/figs`, `outputs/simulations`, and `outputs/tables`.


### Contact

For additional information, please contact Jian Yen (jian.yen [at] delwp.vic.gov.au).

