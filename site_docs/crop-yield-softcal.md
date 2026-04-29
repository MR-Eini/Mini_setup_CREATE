# 3_CropYield / softcal workflow

## Purpose

The `3_CropYield/softcal` folder prepares a SWAT+ setup for later hard calibration by checking and adjusting process-relevant parameters before formal calibration.

The scripts focus on three tasks:

1. adjusting crop days to maturity (`days_mat`) and PHU behavior,
2. adjusting crop-growth parameters to improve simulated yields,
3. adjusting water-yield-related parameters, especially `esco` and optionally `epco`.

The workflow is deliberately semi-automatic. It runs parameter samples, produces diagnostic plots, and then requires manual selection of parameter values based on plausibility and comparison with observed crop-yield ranges.

## Key files and folders

| Path | Role |
|---|---|
| `softcal.Rproj` | RStudio project file. |
| `workflow/01_crop_phu.R` | Samples `days_mat`, runs SWAT+, and supports manual selection of days-to-maturity values. |
| `workflow/02_crop_yield.R` | Samples crop-growth parameters and writes accepted `plants.plt`. |
| `workflow/02_wateryield.R` | Samples `esco` or `esco + epco` and optionally writes `hydrology.hyd`. |
| `observation/crop_yields.csv` | Main observed-yield table used by the scripts. |
| `observation/crop_yield_all.csv` | Extended observed-yield table; not directly used by the current scripts. |
| `./clean_setup` | Expected SWAT+ model folder. |
| `./simulation` | Created/used by SWATrunR for simulation outputs. |
| `./backup` | Expected location for original `plants.plt` and `hydrology.hyd`. |

## Observation table

The main observation file is `observation/crop_yields.csv`.

| Crop | Yield minimum | Yield maximum | Yield mean |
|---|---:|---:|---:|
| `fesc` | 4 | 6 | 5.5 |
| `lupn` | 1.72 | 2.15 | 1.935 |
| `oats` | 2.58 | 3.612 | 2.99925 |
| `pota` | 4.4 | 8.8 | 6.215 |
| `canp` | 1.82 | 3.64 | 2.764125 |
| `wwht` | 3.01 | 4.73 | 3.88075 |

## Requirements

| Requirement | Details |
|---|---|
| R packages | `SWATtunR`, `SWATrunR`, `dplyr`, `tibble`; `02_wateryield.R` also requires `SWATreadR`. |
| SWAT+ model folder | `./clean_setup` must contain a runnable SWAT+ setup. |
| SWAT+ outputs | `mgtout`, `basin_wb_aa`, and `basin_aqu_aa` outputs must be available or requested correctly. |
| Observed crop-yield statistics | `./observation/crop_yields.csv` must contain crop names matching SWAT+ plant names. |
| Backup files | `./backup/plants.plt` and `./backup/hydrology.hyd` are used to restore original parameters. |
| Compute resources | Scripts use `n_cores = 5`, `10`, and `25`; adapt to the machine. |
| Execution mode | Interactive RStudio execution is preferred. |

## Recommended execution order

1. Confirm that `./clean_setup` is runnable.
2. Confirm that `./observation/crop_yields.csv` is correct.
3. Run `workflow/01_crop_phu.R`.
4. Inspect PHU, yield, and biomass plots.
5. Manually select `dmat_sel`.
6. Keep `dmat_sel` in memory or save it before restarting R.
7. Run `workflow/02_crop_yield.R`.
8. Inspect dotty yield plots.
9. Manually select `crop_par_sel`.
10. Run the final crop-parameter check.
11. Overwrite `plants.plt` only when results are acceptable.
12. Run `workflow/02_wateryield.R`.
13. Select alternative A or B.
14. Check water-yield ratio and crop-yield behavior.
15. Write `hydrology.hyd` only after both water-yield and crop-yield checks are acceptable.

## `workflow/01_crop_phu.R`

This script adjusts `days_mat` values so that simulated crop development is compatible with crop characteristics and management schedules.

| Part | What the script does |
|---|---|
| Load packages | Loads `SWATtunR`, `SWATrunR`, `dplyr`, and `tibble`. |
| Define paths | Sets `model_path` and `yield_obs_path`. |
| Select crops | Uses all `plant_name` entries in `crop_yields.csv`. |
| Reset `plants.plt` | Copies the backup file when `reset = TRUE`. |
| Sample `days_mat` | Generates parameter sets using `sample_days_mat()`. |
| Run simulations | Runs SWAT+ and extracts yield, biomass, and PHU from `mgtout`. |
| Load latest run | Loads the latest timestamped `sim_dmat` run. |
| Plot diagnostics | Plots PHU, yield, and biomass against adjusted `days_mat`. |
| Manual selection | Defines `dmat_sel` manually and prepares plant parameters. |

!!! note
    The target is not simply to maximize yield. The first crop-calibration step should make crop maturity plausible.

## `workflow/02_crop_yield.R`

This script assumes that `dmat_sel` already exists from `01_crop_phu.R`.

| Part | What the script does |
|---|---|
| Dependency | Requires `dmat_sel` in the active R environment. |
| Parameter bounds | Defines changes for `lai_pot`, `harv_idx`, `tmp_base`, and `bm_e`. |
| Sampling | Generates 40 Latin Hypercube samples. |
| Add PHU settings | Combines crop-growth samples with selected `days_mat` values. |
| Run simulations | Runs SWAT+ and extracts yield from `mgtout`. |
| Plot dotty figures | Plots yield responses against sampled parameter changes. |
| Manual selection | Defines `crop_par_sel`. |
| Final check | Runs final simulation and plots PHU, yield, and biomass. |
| Write `plants.plt` | Copies accepted `plants.plt` into `clean_setup`. |

## `workflow/02_wateryield.R`

This script calibrates the water-yield ratio using either:

- Alternative A: `esco` only,
- Alternative B: `esco` and `epco` together.

| Part | What the script does |
|---|---|
| Choose alternative | Default is alternative B. |
| Define target | Sets target water-yield ratio to `0.133`. |
| Load observations | Loads crop-yield observations for yield checks. |
| Optional reset | Restores `hydrology.hyd` from backup. |
| Generate parameter set | Samples `esco` or full `esco x epco` grid. |
| Run simulations | Extracts basin precipitation, runoff, lateral flow, tile flow, and aquifer flow. |
| Plot WYR response | Compares simulated water-yield ratio to target. |
| Check selected set | Runs follow-up simulation with crop outputs. |
| Write `hydrology.hyd` | Uses `SWATreadR` to write fixed `esco` and `epco`. |
| Optional calibration ranges | Creates parameter boundaries for later calibration. |

## Outputs

| Output | Produced by | Meaning |
|---|---|---|
| `./simulation/<timestamp>_sim_dmat` | `01_crop_phu.R` | Batch for sampled `days_mat` values. |
| `./simulation/<timestamp>_sim_yld` | `02_crop_yield.R` | Batch for sampled crop-growth parameters. |
| `./simulation/<timestamp>_sim_check01` | `02_crop_yield.R` | Final crop-parameter check run. |
| `./simulation/<timestamp>_sim_wbal` | `02_wateryield.R` | Batch for `esco/epco` combinations. |
| `./simulation/<timestamp>_sim_check02` | `02_wateryield.R` | Final water-yield parameter check run. |
| `./clean_setup/plants.plt` | `02_crop_yield.R` | Final crop-parameter table. |
| `./clean_setup/hydrology.hyd` | `02_wateryield.R` | Final hydrology parameter table, if written. |

## Interpretation rules

- PHU fractions should be plausible for the crop and harvest schedule.
- Simulated crop yields should fall within or near the observed minimum–maximum range unless deviation is justified.
- Biomass should be agronomically plausible and consistent with yield response.
- Water-yield ratio should approach the target value, but selected `esco/epco` values must also preserve crop-yield plausibility.
- Parameter changes should remain biologically and hydrologically defensible.
