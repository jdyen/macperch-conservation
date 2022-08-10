## Using population models to inform conservation prioritisation across multiple populations

This repository contains code and data to support the following manuscript:
Yen JDL et al. (submitted) Establishing new populations in water-secure locations may benefit species persistence more than interventions in water-stressed locations.

Last updated: 10 August 2022

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6378145.svg)](https://zenodo.org/badge/DOI/10.5281/zenodo.6378145.svg)


### Abstract

The effects of climate and land-use change, coupled with limited management resources, mean that managers are increasingly faced with decisions to conserve some populations at the expense of others. Given their potentially controversial nature, such decisions require predictions of likely outcomes under different management strategies and potential future climates. We develop an approach to prioritise management interventions across multiple populations in highly variable and uncertain environments, using a stochastic model of population dynamics structured to include the (dynamic) effects of environmental conditions on key vital rates. We demonstrate this approach with a case study on the endangered Australian freshwater fish, the Macquarie perch (Macquaria australasica), assessing the effectiveness of seven management interventions under five scenarios of climate change. All management interventions were predicted to improve Macquarie perch population outcomes (increased abundances, reduced risk of quasi-extinction), with the risk of declines lowest when management included augmented gene flow and enforcement of fishing regulations. The lowest levels of risk occurred in locations with stable environmental conditions, which suggests that reintroducing Macquarie perch to large, regulated rivers with suitable access to spawning habitat may be more beneficial than investment in small, unregulated streams where populations were unlikely to persist under any scenario. Population reintroductions provide opportunities to shift populations from water-stressed to water-secure locations where required (e.g. to preserve unique genetic diversity). This finding likely applies to many aquatic species, and highlights a potential need to supplement threat mitigation in water-stressed locations with efforts to establish new populations in water-secure locations.


### Usage

The `main.R` script runs the primary analysis and generates associated outputs. The `sensitivity.R` script runs additional sensitivity analyses reported in the main text but presented primarily in supplementary material. Note that these scripts may take several hours to run. Several helper functions are provided in scripts in the `R` directory; these are sourced directly within `main.R` and `sensitivity.R`.

Population models are simulated with the `aae.pop` R package, with the template for the Macquarie perch model included in the `aae.pop.templates` R package. Flow data are downloaded and compiled with the `aae.hydro` R package. These three packages are not available on CRAN but can be installed directly from GitHub with the remotes R package (e.g. `remotes::install_github("aae-stats/aae.pop")`). 

Discharge and water temperature data are provided to avoid downloading and re-calculating flow metrics for each simulation. Downloads from the Victorian WMIS are unreliable for some gauges (the Yarra River, in particular).

Several additional directories will need to be created to store model outputs: `outputs/figs`, `outputs/simulations`, and `outputs/tables`.


### Contact

For additional information, please contact Jian Yen (jian.yen [at] delwp.vic.gov.au).

