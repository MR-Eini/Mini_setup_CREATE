# # If the package remotes is not installed run first:
# install.packages("remotes")
# # The installation of `SWATrunR`.
# remotes::install_github("chrisschuerz/SWATrunR@remove_legacy", dependencies = TRUE)
# # The installation of `SWATtunR`.
# remotes::install_github("biopsichas/SWATtunR")

#Check here for more details
#https://biopsichas.github.io/SWATtunR/  

## Required libraries to run workflow
library(SWATtunR)
library(SWATrunR)
library(tidyverse)
library(tibble)
library(purrr)


# initialize_softcal(project_name = 'softcal_CS6', 
#                    path = 'A:/Users/Eini/CS6/3_Crop/Calibration_crop/', 
#                    model_path = 'A:/Users/Eini/CS6/3_Crop/Calibration_crop/clean_setup')


# Parameter definition ----------------------------------------------------
# Path to the SWAT+ project folder.
model_path <- "./clean_setup"

# Set the number of cores for parallel model execution
n_cores <- 5 # Inf uses all cores. Set lower value if preferred.

# Set the number parameter combinations for the LHS sampling of crop parameters
n_combinations <- 10

# Path to the observed crop yields.
# This file must be updated with case study specific records!
yield_obs_path <- './observation/crop_yields.csv'

# Load and prepare data ---------------------------------------------------
# Load the yield observations
yield_obs  <- read.csv(yield_obs_path)

# Define the crops which should be used in the calibration.
# Default is all crops which are defined in yield_obs.
# Please define manually if only selected crops should be considered.
crop_names <- yield_obs$plant_name

# Optional reset of plants.plt --------------------------------------------
# In the case the crop calibration workflow should be redone after the last step
# of this script was already executed and the plants.plt was overwritten the
# plants.plt should be reset to its initial condition. To perform the reset set
# reset <- FALSE
reset <- FALSE
if(reset) {
  file.copy('./backup/plants.plt',
            paste0(model_path, '/plants.plt'),
            overwrite = TRUE)
} else if (!file.exists('./backup/plants.plt')){
  file.copy(paste0(model_path, '/plants.plt'),
            './backup/plants.plt',
            overwrite = TRUE)
}


par_dmat <- sample_days_mat(crop_names)

# Run the SWAT+ model with the generated days_mat parameter set
run_swatplus(project_path = model_path,
             output = list(yld = define_output(file = 'mgtout',
                                               variable = 'yld',
                                               label = crop_names),
                           bms = define_output(file = 'mgtout',
                                               variable = 'bioms',
                                               label = crop_names),
                           phu = define_output(file = 'mgtout',
                                               variable = 'phu',
                                               label = crop_names)
             ),
             parameter        = par_dmat,
             start_date       = NULL, # Change if necessary.
             end_date         = NULL, # Change if necessary.
             years_skip       = NULL, # Change if necessary.
             n_thread         = n_cores,
             save_path        = './simulation',
             save_file        = add_timestamp('sim_dmat'),
             return_output    = FALSE,
             time_out         = 3600 # seconds, change if run-time differs
)


# Load the most recent dmat simulation results
dmat_sims <- list.files('./simulation/', pattern = '[0-9]{12}_sim_dmat')
dmat_path <- paste0('./simulation/', dmat_sims[length(dmat_sims)])
ylds_phu_dmat <- load_swat_run(dmat_path, add_date = FALSE)

# Plot PHU, crop yields and biomass over adjusted days to maturity values.
plot_phu_yld_bms(ylds_phu_dmat, yield_obs)


# Set days to maturity values for all selected crops based on the figure above.
dmat_sel <- tibble(
  plant_name = c("rye", "oats", "canp", "fesc", "wwht"),
  'days_mat.pdb | change = absval' = c(140, 140, 140, 100, 100)
)

# Check if user defined days to maturity values for all crops.
stopifnot(all(crop_names %in% dmat_sel$plant_name))
# Update names of dmat_sel to be used as SWATrunR parameters
dmat_sel <- prepare_plant_parameter(dmat_sel)

# Additional parameters
par_bnd <- tibble('lai_pot.pdb | change = relchg'  = c(-0.3, 0.3),
                  'harv_idx.pdb | change = relchg' = c(-0.3, 0.3),
                  'tmp_base.pdb | change = abschg' = c(-5.5, 5.5),
                  'bm_e.pdb | change = relchg'     = c(-0.3, 0.1))

## The number of samples can be adjusted based on the available computational resources.
## Recommended number of samples is 50-100.
n_combinations <- 50
par_crop <- sample_lhs(par_bnd, n_combinations)
# Add updated days to maturity values to parameter set
par_crop <- bind_cols(par_crop, dmat_sel)

# Run the SWAT+ model with the additional parameter set
run_swatplus(project_path = model_path,
             output = list(yld = define_output(file = 'mgtout',
                                               variable = 'yld',
                                               label = crop_names)),
             parameter = par_crop,
             start_date       = NULL, # Change if necessary.
             end_date         = NULL, # Change if necessary.
             years_skip       = NULL, # Change if necessary.
             n_thread         = n_cores,
             save_path        = './simulation',
             save_file        = add_timestamp('sim_yld'),
             return_output    = FALSE,
             time_out         = 3600 # seconds, change if run-time differs
)

# Load the most recent yield simulation results
yld_sims <- list.files('./simulation/', pattern = '[0-9]{12}_sim_yld')
yld_path <- paste0('./simulation/', yld_sims[length(yld_sims)])
yld_sim  <- load_swat_run(yld_path, add_date = FALSE)
# Remove days to maturity parameter columns before plotting.
yld_sim$parameter$values <- yld_sim$parameter$values[, 1:4]

## Plot dotty figures for the selected crops
plot_dotty_yields(yld_sim, yield_obs)

######################

# Fix the parameter changes you want to apply to the crops
crop_par_sel <- tibble::tibble(
  plant_name = 						 c("rye", "oats", "canp", "fesc", "wwht"),   #Adjust
  'bm_e.pdb | change = relchg'     = c(-0.1,  0.1,      0,       0.09,   -0.21),
  'harv_idx.pdb | change = relchg' = c(0.1,     0,      0.24,    0.25,    0.26),
  'lai_pot.pdb | change = relchg'  = c(0.3,   0.3,      0.07,    0.28,    0.28),
  'tmp_base.pdb | change = abschg' = c( -1,    -1,      1,       -0.3,   -0.14)
)

# Check if user defined days to maturity values for all crops.
stopifnot(all(crop_names %in% crop_par_sel$plant_name))
# Restructure the set parameter changes to SWATrunR
crop_par_sel <- prepare_plant_parameter(crop_par_sel)

######################
# Run the simulations
run_swatplus(project_path = model_path,
             output = list(yld = define_output(file = 'mgtout',
                                               variable = 'yld',
                                               label = crop_names),
                           bms = define_output(file = 'mgtout',
                                               variable = 'bioms',
                                               label = crop_names),
                           phu = define_output(file = 'mgtout',
                                               variable = 'phu',
                                               label = crop_names)
             ),
             parameter        = crop_par_sel,
             start_date       = NULL, # Change if necessary.
             end_date         = NULL, # Change if necessary.
             years_skip       = NULL, # Change if necessary.
             n_thread         = n_cores,
             save_path        = './simulation',
             save_file        = add_timestamp('sim_check01'),
             return_output    = FALSE,
             time_out         = 3600, # seconds, change if run-time differs
             keep_folder      = TRUE
)

# Load the most recent check simulation results
check_sims <- list.files('./simulation/', pattern = '[0-9]{12}_sim_check01')
check_path <- paste0('./simulation/', check_sims[length(check_sims)])
check_sim  <- load_swat_run(check_path, add_date = FALSE)

# Plot PHU, crop yields and biomass for final simulation run.
plot_phu_yld_bms(check_sim, yield_obs, 0.3)

# Write ‘plants.plt’
file.copy(paste0(model_path, '/.model_run/thread_1/plants.plt'), model_path,
          overwrite = TRUE)
unlink(paste0(model_path, '/.model_run'), recursive = TRUE)




