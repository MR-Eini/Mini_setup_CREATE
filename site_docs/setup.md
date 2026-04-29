# 1_Setup workflow

## Purpose

The `1_Setup` folder is an end-to-end SWAT+ setup generation workflow. Its purpose is to regenerate a complete SWAT+ model setup from pre-processed GIS, soil, weather, and management inputs.

The workflow:

1. builds a base SWAT+ project,
2. adds climate and deposition data,
3. writes SWAT+ text input files,
4. modifies selected SWAT+ input tables,
5. prepares management operations with SWATfarmR,
6. runs SWAT+,
7. exports a clean `TxtInOut`-style setup folder.

!!! warning "Input preparation is outside this workflow"
    DEM, soil, land-use, river, basin, weather, crop, and management files must already be prepared and checked before running this folder.

## Key files

| Item | Description |
|---|---|
| Main driver | `setup_workflow.R` |
| Main configuration | `settings.R` |
| Utility functions | `functions.R` |
| Base setup builder | `Libraries/buildr_script/swatbuildr.R` |
| Management input preparation | `Libraries/farmR_input/write_SWATfarmR_input.R` |
| SWAT+ writer executable | `Libraries/write.exe` |
| SWAT+ model executable | Configured by `swat_exe` in `settings.R` |
| Default final output | `Temp/clean_setup/`, depending on `res_path` |

## Expected folder structure

```text
1_Setup/
├── README.md
├── setup_workflow.R
├── settings.R
├── functions.R
├── wf.Rproj
├── Data/
│   ├── for_buildr/
│   ├── for_prepr/
│   └── for_farmr_input/
├── Libraries/
│   ├── write.exe
│   ├── SWATp_jan_sept.exe
│   ├── buildr_script/
│   ├── farmR_input/
│   ├── files_to_overwrite_at_the_end/
│   └── calibration_cal/
└── Temp/                  # generated output; not source material
    ├── buildr_project/
    ├── farmR_input/
    └── clean_setup/
```

## Required software and packages

| Requirement | Role |
|---|---|
| R 4.x | Runs all R scripts. |
| RStudio | Recommended interactive environment. |
| SWAT+ executable | Runs the model in the setup workflow. |
| `write.exe` | Writes SWAT+ text input files from the SQLite database. |
| WhiteboxTools | Used by SWATbuildR for terrain and connectivity analyses. |
| SWAT+ Editor or DB Browser for SQLite | Optional inspection and troubleshooting. |

Packages used directly or indirectly include:

- `remotes`
- `RNetCDF`
- `tidyverse`
- `mapview`
- `sf`
- `gstat`
- `rstudioapi`
- `whitebox`
- `DBI`
- `RSQLite`
- `vroom`
- `lubridate`
- `reshape2`
- `data.table`
- `HighFreq`
- `SWATtunR`
- `SWATprepR`
- `SWATfarmR`
- `SWATdoctR`
- `SWATmeasR`
- `SWATreadR`

## Important settings in `settings.R`

| Variable | Example value / role |
|---|---|
| `swat_exe` | Name of SWAT+ executable in `Libraries`, e.g. `SWATp_jan_sept.exe`. |
| `res_path` | Output folder, e.g. `Temp`. |
| `data_path` | Root folder for pre-processed inputs, e.g. `Data`. |
| `lib_path` | Root folder for helper scripts and executables, e.g. `Libraries`. |
| `st_year` | Simulation start year, e.g. `2004`. |
| `end_year` | Simulation end year, e.g. `2023`. |
| `weather_path` | RDS weather object used by SWATprepR. |
| `pnt_path` | Point-source template. Set to `NULL` to skip point-source processing. |
| `lab_p` | Single labile phosphorus value written into `nutrients.sol`. |

## Main workflow steps

| Step | Name | Main action |
|---:|---|---|
| 1 | Initialize workflow | Load settings, functions, packages, WhiteboxTools, and create a clean result folder. |
| 2 | Run SWATbuildR | Build the base SWAT+ project and SQLite database. |
| 3 | Back up database | Create a zipped backup of the generated SQLite database. |
| 4 | Add weather data | Add weather and WGN information using SWATprepR. |
| 5 | Patch SQLite | Set `project_config$input_files_dir` to `"."` for `write.exe`. |
| 6 | Write SWAT+ text files | Copy and run `write.exe`. |
| 7 | Check land connectivity | Create a routing-unit connection shapefile for visual checking. |
| 8 | Add atmospheric deposition | Add annual atmospheric deposition inputs. |
| 9 | Link aquifers and channels | Create or update `aqu_cha.lin`. |
| 10 | Add point source data | Apply point-source template if `pnt_path` is not `NULL`. |
| 11 | Prepare SWATfarmR input | Generate `farmR_input.csv` and related check files. |
| 12 | Update `landuse.lum` | Modify land-use pointer columns using project-specific prefix rules. |
| 13 | Update nutrients and HRU data | Modify `nutrients.sol` and `hru-data.hru`. |
| 14 | Update `time.sim` | Write selected simulation years. |
| 15 | Run SWAT+ setup | Run the model before final management-file generation. |
| 16 | Generate management files | Use SWATfarmR to write management operations. |
| 17 | Fix unconnected reservoirs | Modify reservoir connectivity and hydrology defaults. |
| 18 | Optional edits | Placeholder for additional SWAT+ file edits. |
| 19 | Run final SWAT+ setup | Validate the modified setup by rerunning SWAT+. |
| 20 | Export clean setup | Copy filtered input files into `Temp/clean_setup`. |
| 21 | Optional `calibration.cal` | Disabled by default; enable only with a valid calibration file. |

## Outputs

| Output | Created by | Purpose |
|---|---|---|
| `Temp/buildr_project/` | SWATbuildR | Base generated SWAT+ project and database. |
| `Temp/db_backup.zip` | Step 3 | Archive of the generated SQLite database. |
| Text input folder under `buildr_project/Mini_CREATE` | `write.exe` | Working SWAT+ input folder. |
| `land_connections_as_lines.shp` | Step 7 | Visual inspection of land routing connectivity. |
| `Temp/farmR_input/` | Step 11 | Saved SWATfarmR input and check files. |
| `Temp/clean_setup/` | Step 20 | Final clean input-only setup for later stages. |

## Fragile points

- Do not keep SQLite files open in SWAT+ Editor, DB Browser, or Windows Explorer preview while `write.exe` runs.
- The setup stops if zero or more than one matching project SQLite database is found.
- `landuse.lum` prefix rules are project-specific and must be reviewed before reuse.
- The current `lab_p` implementation applies one value catchment-wide.
- Some helper functions rely on global variables and should be checked before teaching.
- Package installation during workshops is risky; pre-install and freeze package versions when possible.

