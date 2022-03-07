
x <- TRACERx_NEJM_2017_REVOLVER$dataset[[2]]
x <- TRACERx_NEJM_2017

devtools::install_github("caravagn/evoverse.datasets")
library(evoverse.datasets)
data('TRACERx_NEJM_2017', package = 'evoverse.datasets')
available_data()