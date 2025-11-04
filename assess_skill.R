# - Write a script to compare the predictive skill of two seasonal precipitation forecasts: one is raw NMME data, the other is bias-adjusted.
# - Comparisons should be on a per-gridcell basis, meaning, the result should always be grids (rasters).
# - We're going to use the {stars}. You have access to its documentation.
# - The files of the forecasts have a third dimension "L", which represents the lead month
# - Comparisons should be only for the month of July.
# - We're going to compare the skill of different NMME models
# - We want to compare the mean absolute difference between forecasted value and observed. Observations come from ERA5
# - We also want to compare their anomalies in relation to climatology, which will be provided.
# - Use the functions in the general_tools.R file to download necessary files to the location "/mnt/pd-ssd-50/residuals/"
# - Observational data is found in this google bucket location: "gs://clim_data_reg_useast1/era5/monthly_means/total_precipitation/"
# - An example of a filename is: "era5_total-precipitation_mon_2021-01-01.nc"
# - Observational climatologies are in "gs://clim_data_reg_useast1/era5/monthly_means/climatologies/"
# - An example of a file name is "era5_total-precipitation_mon_gamma-params_1991-2020_02.nc" -- the "02" corresponds to the month.
# - Note that 6 months of observational data ahead of July need to be processed to be able to measure skill for all lead times
# - The raw NMME data is in "gs://clim_data_reg_useast1/nmme/monthly/<model name>/precipitation/"
# - An example of a filename is "nmme_<model name>_precipitation_mon_ic-2021-01-01_leads-6.nc"
# - The bias adjusted NMME data is in "gs://clim_data_reg_useast1/nmme/monthly/<model name>/precipitation_biasadj/"
# - An example of a filename is "nmme_<model name>_precipitation_mon_ic-2021-01-01_leads-6_biasadj.nc"
# - The NMME climatologies are in "gs://clim_data_reg_useast1/nmme/climatologies/<model name>/".
# - An example of a filename is "nmme_<model name>_precipitation_mon_gamma-params_1991-2020_01_leads-6.nc
# - The observational and NMME climatologies consist of stars objects with two variables, each representing a parameter of a gamma distribution.
# - We only want to run the analysis from 2021 to 2024 (4 years)
