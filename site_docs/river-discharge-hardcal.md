# 4_RiverDischarge / hardcal workflow

## Purpose

The `4_RiverDischarge/hardcal` folder is designed for hard calibration of river discharge and selected in-stream water-quality variables using SWATrunR and SWATtunR.

The workflow runs many SWAT+ simulations with different hydrological parameter combinations and evaluates the runs against observed discharge and/or water-quality data.

!!! warning "Template status"
    The scripts are templates. Several parts require case-study-specific editing. Some scripts are not directly runnable without additional input files or earlier objects loaded in the R session.

## Files and workflow scope

| File / folder | Role |
|---|---|
| `hardcal.Rproj` | RStudio project file. |
| `workflow/01_define_parameter.R` | Defines calibration parameters and samples 900 LHS combinations. |
| `workflow/02_define_output.R` | Defines SWATrunR outputs to extract from SWAT+ simulations. |
| `workflow/03_run_swat.R` | Runs SWAT+ over 2004–2023 and starts printed outputs from 2007. |
| `workflow/04_analyze_results.R` | Single-station analysis template for discharge and water-quality variables at channel 5. |
| `workflow/04_Multi_stations_analyze_results.R` | Weighted multi-gauge discharge-only analysis for channels 44, 72, and 76. |
| `workflow/04_Multi_stations_analyze_results_no_w.R` | Unweighted multi-gauge discharge-only analysis. |
| `workflow/04_Multi_stations_analyze_results_final_all.R` | Integrated multi-gauge, multi-variable workflow for Q, NO3, PO4-P, and sediment. |
| `workflow/functions.R` | Local helper functions. |
| `workflow/05_validate.R` | Validation template for selected parameter sets. |
| `observation/*.csv` | Observed water-quality concentration files for channel 5. |

## Expected folder structure

```text
4_RiverDischarge/hardcal/
├── hardcal.Rproj
├── clean_setup/                  # SWAT+ TxtInOut-style project folder
├── workflow/
│   ├── 01_define_parameter.R
│   ├── 02_define_output.R
│   ├── 03_run_swat.R
│   ├── 04_analyze_results.R
│   ├── 04_Multi_stations_analyze_results.R
│   ├── 04_Multi_stations_analyze_results_no_w.R
│   ├── 04_Multi_stations_analyze_results_final_all.R
│   ├── 05_validate.R
│   └── functions.R
├── observation/
│   ├── TN_cha5.csv
│   ├── TP_cha5.csv
│   ├── no3n_cha5.csv
│   └── suspended_sediment_cha5.csv
└── simulation/                   # created by SWATrunR
```

## Required packages

| Dependency | Purpose |
|---|---|
| `SWATrunR` | Runs SWAT+ with parameter tables and loads saved simulation folders. |
| `SWATtunR` | Provides sampling, GOF calculation, FDC analysis, ranking, and calibration-file writing. |
| `hydroGOF` | Provides NSE, KGE, PBIAS, and MAE functions. |
| `tidyverse`, `dplyr`, `readr`, `tibble`, `lubridate` | Data manipulation, date handling, joins, and plotting. |
| `zoo` | Optional observation gap filling. |
| SWAT+ executable and complete `clean_setup` | Required for runnable simulations. |

## Overall calibration logic

1. Define calibration parameter boundaries.
2. Generate sampled parameter combinations.
3. Define SWAT+ output files, variables, and spatial units to extract.
4. Run SWAT+ repeatedly using SWATrunR.
5. Load the most recent simulation folder.
6. Align simulated and observed time series over the calibration period.
7. Calculate objective functions such as NSE, KGE, PBIAS, MAE, and FDC-segment RSR.
8. Rank parameter sets and inspect diagnostic plots.
9. Select candidate runs.
10. Write selected parameter sets to `calibration.cal` or pass them to validation.

## Step 1 — Define calibration parameters

Script: `workflow/01_define_parameter.R`

| Process group | Parameters | Interpretation |
|---|---|---|
| Snow | `snomelt_tmp.hru`, `snofall_tmp.hru` | Relevant when snowfall contributes meaningfully. |
| Evapotranspiration / soil water | `esco.hru`, `awc.sol` | Controls soil evaporation compensation and available water capacity. |
| Surface runoff | `cn2.hru`, `cn3_swf.hru`, `surlag.bsn` | Controls curve-number behavior and runoff lag. |
| Lateral flow | `lat_len.hru`, `latq_co.hru`, `bd.sol`, `k.sol` | Controls lateral-flow length/coefficient and soil hydraulic properties. |
| Percolation / aquifer | `perco.hru`, `flo_min.aqu`, `revap_co.aqu`, `revap_min.aqu`, `alpha.aqu` | Controls percolation, aquifer recession/revap, and baseflow. |
| Nutrient/sediment parameters | Commented-out examples | Optional nutrient, channel, and phosphorus parameters. |

The script uses:

```r
parameter_set <- sample_lhs(parameter_boundaries, n_combinations)
```

with `n_combinations = 900`.

### HRU-group-specific treatment

| Parameter | Current treatment | Practical implication |
|---|---|---|
| `perco.hru` | Groups are created, but translation call is commented out. | Remains controlled by the global absolute range unless activated. |
| `cn3_swf.hru` | Translated into low/mod/high group-specific ranges. | Original normalized column is replaced by group-specific columns. |
| `latq_co.hru` | Translated into low/mod/high group-specific ranges. | Original normalized column is replaced by group-specific columns. |

## Step 2 — Define SWAT+ outputs

Script: `workflow/02_define_output.R`

The current setup extracts daily channel outputs for `cha_ids = 5`.

| Output object | SWAT+ output file | Variable | Meaning |
|---|---|---|---|
| `flo_day` | `channel_sd_day` | `flo_out` | Daily discharge. |
| `no3_day` | `channel_sd_day` | `no3_out` | Daily nitrate output. |
| `orgn_day` | `channel_sd_day` | `orgn_out` | Daily organic nitrogen output. |
| `nh3_day` | `channel_sd_day` | `nh3_out` | Daily ammonia output. |
| `no2_day` | `channel_sd_day` | `no2_out` | Daily nitrite output. |
| `solp_day` | `channel_sd_day` | `solp_out` | Daily soluble phosphorus output. |
| `sedp_day` | `channel_sd_day` | `sedp_out` | Daily sediment-attached phosphorus output. |
| `sed_day` | `channel_sd_day` | `sed_out` | Daily sediment output. |

## Step 3 — Run SWAT+ simulations

Script: `workflow/03_run_swat.R`

| Setting | Current value | Explanation |
|---|---|---|
| `model_path` | `./clean_setup` | SWAT+ project folder. |
| `start_date` | `2004-01-01` | Start of model simulation. |
| `end_date` | `2023-12-31` | End of model simulation. |
| `start_date_print` | `2007-01-01` | First day for printed outputs; 2004–2006 are warm-up. |
| `n_cores` | `10` | Parallel execution threads. |
| `save_path` | `./simulation` | Folder for timestamped simulation results. |
| `split_units` | `FALSE` | Script notes that `TRUE` may be better for many output units. |
| `time_out` | `3600` seconds | Per-run timeout. |

```r
run_swatplus(
  project_path     = model_path,
  output           = outputs,
  parameter        = parameter_set,
  start_date       = start_date,
  end_date         = end_date,
  start_date_print = start_date_print,
  n_thread         = n_cores,
  save_path        = save_path,
  save_file        = save_file_name,
  return_output    = FALSE,
  split_units      = FALSE,
  time_out         = 3600
)
```

## Observation data

The inspected observation folder contains four water-quality CSV files for channel 5. Each file has two columns: `date` and `value`.

| File | Variable | Use |
|---|---|---|
| `TN_cha5.csv` | Total nitrogen concentration | Used as `ntot_path`. |
| `TP_cha5.csv` | Total phosphorus concentration | Used as `ptot_path`. |
| `no3n_cha5.csv` | Nitrate-nitrogen concentration | Referenced in comments; not active in the single-station script. |
| `suspended_sediment_cha5.csv` | Suspended sediment concentration | Used as `sused_path`. |

!!! warning "Missing discharge observation files"
    The single-station and multi-station scripts reference discharge CSV files such as `q_cha5_cms.csv`, `q_cha44_cms.csv`, `q_cha72_cms.csv`, and `q_cha76_cms.csv`. These were not found in the inspected non-MP file list described in the manual.

## Local helper functions

Script: `workflow/functions.R`

| Function | Purpose |
|---|---|
| `as_date()` | Converts input to `Date`. |
| `filter_period()` | Filters a data frame by date. |
| `norm_run()` | Normalizes run labels to `run_0001` style. |
| `run_to_id()` | Extracts integer IDs from run labels. |
| `limit_runs()` | Keeps date plus the first N run columns. |
| `read_obs()` | Reads an observation CSV. |
| `fill_obs()` | Optionally fills observation gaps. |
| `align_sim_obs()` | Aligns simulation and observation tables. |
| `prep_for_fdc()` | Prepares tables for flow-duration-curve calculations. |
| `safe_fdc_rsr()` | Runs `calc_fdc_rsr()` inside `tryCatch`. |
| `rn()` | Rank-normalizes a vector to `[0,1]`. |
| `choose_gauge_key()` | Accepts a gauge name or index and returns a valid key. |

## Analysis scripts

### Single-station analysis

Script: `workflow/04_analyze_results.R`

This template extracts discharge, nitrogen, phosphorus, and suspended sediment at channel 5. It calculates NSE, KGE, PBIAS, MAE, and discharge FDC RSR, then prepares diagnostic plots and selected calibration files.

### Weighted multi-station discharge analysis

Script: `workflow/04_Multi_stations_analyze_results.R`

| Setting | Value |
|---|---|
| `N_RUNS` | 600 |
| `TOP_K` | 600 |
| `ALIGN_MODE` | `overlap_seq` |
| `FILL_METHOD` | `none` |
| Gauges | `cha44`, `cha72`, `cha76` |
| Period | 2004-01-01 to 2013-12-31 |
| Weights | `cha44 = 0.6`, `cha72 = 0.2`, `cha76 = 0.2` |

### Unweighted multi-station analysis

Script: `workflow/04_Multi_stations_analyze_results_no_w.R`

This averages normalized scores across gauges instead of applying explicit weights. It duplicates helper functions inside the script rather than sourcing `workflow/functions.R`.

### Integrated multi-variable analysis

Script: `workflow/04_Multi_stations_analyze_results_final_all.R`

| Variable | Simulation treatment | Observation expectation |
|---|---|---|
| Q | Simulated discharge directly from `flo_day_*`. | m³ s⁻¹ |
| NO3 | Simulated `no3_day_*` load converted to concentration using discharge. | mg L⁻¹ |
| PO4P | Simulated `solp_day_*` load converted to concentration using discharge. | mg L⁻¹ |
| SED | Simulated `sed_day_*` load converted to concentration using discharge. | mg L⁻¹ |

## Validation workflow

Script: `workflow/05_validate.R`

This is a template, not a ready-to-run script.

| Item | Current state | Required action |
|---|---|---|
| `period_valid` | Placeholder | Replace with validation dates. |
| `start_date` | Placeholder | Start at least two years before validation output begins. |
| `parameter_set_valid` | `parameter_set[run_ids, ]` | Requires `parameter_set` and `run_ids`. |
| `outputs` | `outputs` | Requires `workflow/02_define_output.R` or equivalent object. |
| `flow_path` | Placeholder | Replace with a real observation file. |

## Recommended calibration strategy

1. Decide whether the demonstration is single-station or multi-station.
2. Update `02_define_output.R` so `cha_ids` matches the intended analysis script.
3. Check that all observation CSV paths used by the chosen analysis script exist.
4. Run `01_define_parameter.R` and inspect `parameter_set`.
5. Test `03_run_swat.R` with 5–10 parameter combinations before launching 900 runs.
6. Load the newest simulation folder and verify expected objects in `sim$simulation`.
7. Run the analysis only after unit consistency is verified.
8. Select candidate runs.
9. Write `calibration.cal` or pass selected parameter sets to validation.

