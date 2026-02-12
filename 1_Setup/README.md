# SWAT+ Uncalibrated Setup Preparation Workflow

Workflow to regenerate a complete **SWAT+** model setup from *pre-processed inputs* to a **pre-calibrated / pre-prepared** (and optionally calibrated) setup in (ideally) a single run.

This repository is centered around the script:

- `setup_workflow.R` — main end-to-end workflow driver  
- `settings.R` — project configuration (paths, years, SWATbuildR and SWATfarmR settings)

> The preparation of input data is **not** part of this workflow. You must provide the required pre-processed layers/tables externally.

---

## Contents

- [What this workflow does](#what-this-workflow-does)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Configuration](#configuration)
- [How to run](#how-to-run)
- [Outputs](#outputs)
- [Optional steps](#optional-steps)
- [Troubleshooting](#troubleshooting)
- [Reproducibility notes](#reproducibility-notes)
- [Credits](#credits)

---

## What this workflow does

The `setup_workflow.R` script automates the following steps (mirrors the numbered blocks in the script):

1. **Initialize**
   - Loads R packages (`SWATprepR`, `SWATfarmR`, `SWATtunR`, `SWATdoctR`)
   - Sources `settings.R` and `functions.R`
   - Recreates a clean results directory (`res_path`)

2. **Run SWATbuildR**
   - Executes `Libraries/buildr_script/swatbuildr.R` to build a base SWAT+ project from GIS inputs

3. **Add weather + weather generator stats**
   - Loads interpolated weather (`weather_path`, RDS)
   - Computes WGN statistics (`prepare_wgn`)
   - Writes weather (and atmospheric deposition, as supported by SWATprepR) to the model SQLite database (`add_weather`)

4. **Patch SQLite for `write.exe`**
   - Sets `project_config.input_files_dir` to `"."` to ensure `write.exe` can write text inputs

5. **Write SWAT+ text inputs (`TxtInOut`)**
   - Copies and runs `write.exe` via `exe_copy_run()`

6. **Add atmospheric deposition**
   - Downloads deposition data (`get_atmo_dep(...)`)
   - Writes it into the text-based setup (`add_atmo_dep(..., t_ext="annual")`)

7. **Link aquifer ↔ channels (geomorphic flow)**
   - Creates/updates `aqu_cha.lin` through `link_aquifer_channels()`

8. **(Optional) Point sources**
   - Template-based point sources support is scaffolded but disabled by default in the script

9. **Prepare SWATfarmR input**
   - Runs `Libraries/farmR_input/write_SWATfarmR_input.R`
   - Collects generated `.csv` files into `Temp_/farmR_input/`

10. **(Optional) Modify `farmR_input.csv`**
   - Placeholder section to duplicate schedules for drained areas (commented out)

11. **Update `landuse.lum`**
   - Backs up original (`landuse.lum.bak`)
   - Runs `Libraries/read_and_modify_landuse_lum.R`

12. **Update `nutrients.sol` + connect to HRUs**
   - Updates labile phosphorus (`lab_p`) in `nutrients.sol`
   - Sets `hru-data.hru` to use `soilplant1` initialization

13. **Update `time.sim`**
   - Enforces simulation start/end years (`st_year`, `end_year`)

14. **Create connectivity lines shapefile**
   - Generates `land_connections_as_lines.shp` for visual/manual checks of routing connectivity

15. **Run SWAT+ once**
   - Copies and runs SWAT executable (`swat_exe`) via `exe_copy_run()`

16. **Generate management operations with SWATfarmR**
   - Creates a `.farm` project
   - Adds an API (antecedent precipitation index) variable using `variable_decay()`
   - Reads management (`farmR_input.csv`), schedules ops, writes management files for the simulation period

17. **Fix unconnected reservoirs**
   - Backs up reservoir files and enforces consistent connectivity defaults

18. **(Optional) Other file edits**
   - Template section for patching additional SWAT+ input files

19. **Run SWAT+ again**
   - Executes a final run after all modifications

20. **Export a clean setup**
   - Copies only SWAT+ **input** files into `Temp_/clean_setup/` (filters out outputs, zips, sqlite, etc.)

21. **(Optional) Add `calibration.cal`**
   - Script contains a hard `stop()` to prevent accidental use; remove it if you want to copy a calibration file

---

## Requirements

### Software
- **R** (recommended: R 4.x)
- **SWAT+ executable** (the model binary you want to run; configured as `swat_exe` in `settings.R`)
- **write.exe** (SWAT+ input writer, used to generate `TxtInOut` from the SQLite database)

### R packages
Installed from GitHub (as indicated in `setup_workflow.R` comments):

- `SWATtunR`
- `SWATprepR`
- `SWATfarmR`
- `SWATrunR`
- `SWATdoctR`
- `swatmeasr` (from UFZ Git)

Example installation snippet:

```r
install.packages("remotes")

remotes::install_github("biopsichas/SWATtunR")
remotes::install_github("biopsichas/SWATprepR")
remotes::install_github("tkdweber/euptf2")
remotes::install_github("chrisschuerz/SWATfarmR")
remotes::install_github("chrisschuerz/SWATrunR")
remotes::install_github("biopsichas/SWATdoctR")
remotes::install_git("https://git.ufz.de/schuerz/swatmeasr")
```

> Note: `setup_workflow.R` expects **specific SWATfarmR versions** (>= 4.* or special 3.2.0 build referenced in the script). See [Troubleshooting](#troubleshooting).

---

## Repository structure

A typical folder layout expected by the scripts:

```
.
├── setup_workflow.R
├── settings.R
├── functions.R                      # required (sourced by setup_workflow.R)
├── Data/
│   ├── for_buildr/
│   │   ├── DEM1.tif
│   │   ├── soil1.tif
│   │   ├── Soil_SWAT_cod.csv
│   │   ├── usersoil_lrew.csv
│   │   ├── land1.shp (+ sidecars)
│   │   ├── river1.shp (+ sidecars)
│   │   └── basin1.shp (+ sidecars)
│   ├── for_prepr/
│   │   └── met_int.rds              # weather (interpolated), RDS
│   └── for_farmr_input/
│       ├── crops1.shp (+ sidecars)
│       ├── mgt_crops.csv
│       └── mgt_generic.csv
└── Libraries/
    ├── buildr_script/
    │   └── swatbuildr.R
    ├── farmR_input/
    │   └── write_SWATfarmR_input.R
    ├── files_to_overwrite_at_the_end/
    │   ├── plants.plt               # overwritten in step 9
    │   └── (optional other files)
    ├── read_and_modify_landuse_lum.R
    ├── create_connectivity_line_shape.R
    ├── calibration_cal/
    │   └── calibration1.cal         # optional
    ├── write.exe
    └── SWATp_jan_sept.exe           # example; name must match swat_exe
```

Adjust names/paths in `settings.R` if your structure differs.

---

## Configuration

All user-editable settings are in `settings.R`. Key parameters:

### Global
- `res_path` — results folder (default: `Temp_`)
- `data_path` — input data folder (default: `Data`)
- `lib_path` — scripts/executables folder (default: `Libraries`)
- `st_year`, `end_year` — simulation period written into `time.sim`
- `weather_path` — path to `met_int.rds` (interpolated weather RDS)
- `lab_p` — labile phosphorus value written into `nutrients.sol` for the catchment

### SWATbuildR inputs
- `project_name`, `project_path`, `txt_path`
- `dem_path`, `soil_layer_path`, `soil_lookup_path`, `soil_data_path`
- `land_path`, `channel_path`, `basin_path`
- outlet definition: `id_cha_out` (or `id_res_out`)
- routing simplification: `frc_thres`
- `wetland_landuse`, `max_point_dist`

### SWATfarmR input preparation
- `lu_shp`, `mgt_csv`, `lu_generic_csv`
- `hru_crops` — HRU land-use prefix that indicates cropland HRUs
- multi-year grass settings: `m_yr_sch_existing`, `crop_myr`, `max_yr`, `crop_s`, `additional_h_yr_sch_existing`

---

## How to run

From the repository root:

### Option A: Interactive (RStudio)
1. Open the project folder.
2. Edit `settings.R` to match your data paths and years.
3. Run:

```r
source("setup_workflow.R")
```

### Option B: Command line
```bash
Rscript setup_workflow.R
```

Important runtime notes:
- `setup_workflow.R` **deletes** the existing `res_path` directory if it exists.
- The script searches for a single SQLite DB named `<project_name>.sqlite` inside `res_path`. If multiple exist, it stops.

---

## Outputs

The workflow writes results under `res_path` (default `Temp_`), including:

- `Temp_/buildr_project/...` — SWATbuildR project artifacts
- `Temp_/db_backup.zip` — backup of the generated SQLite DB (copied into `res_path` then zipped)
- `Temp_/farmR_input/` — exported SWATfarmR input CSV files
- `Temp_/clean_setup/` — final cleaned SWAT+ input setup, containing only SWAT+ input files (no outputs/sqlite/zips)

The script prints the final clean setup location at the end.

---

## Optional steps

### Point sources
Point source support exists but is commented out in `setup_workflow.R`. To enable it:
1. Provide a point source template file (e.g., `Data/for_prepr/pnt_data.xlsx`)
2. Set `pnt_path` in `settings.R`
3. Uncomment the point source block in step 8

### Add `calibration.cal`
Step 21 is guarded by:

```r
stop("Remove this if you have calibration.cal file")
```

Remove that line to allow copying a calibration file into `clean_setup/` and updating `file.cio`.

---

## Troubleshooting

### “You have more than one database named `<project>.sqlite`…”
The script expects exactly **one** SQLite database that matches `project_name`. Delete/rename duplicates in `res_path` or set `db_path` manually.

### “No database found!”
SWATbuildR did not generate the database, or `res_path`/`project_name` is inconsistent. Confirm:
- `settings.R` values
- SWATbuildR ran successfully (`Libraries/buildr_script/swatbuildr.R`)
- output is written under your expected `res_path`

### `write.exe` fails
The workflow applies a known fix (step 4) by setting `input_files_dir="."`. If it still fails:
- try writing inputs via **SWAT+ Editor** as a workaround
- check that `write.exe` is present in `Libraries/` and is compatible with your SWAT+ build

### SWATfarmR version error
The script supports:
- SWATfarmR **4.x**
- SWATfarmR **3.2.0** (special package referenced in the script)

If you get the stop message about version, install the required version accordingly.

### Reservoir connectivity issues
Step 17 forces defaults for reservoirs with missing/invalid connections. If you still see reservoir routing problems:
- inspect `reservoir.con`, `reservoir.res`, `hydrology.res`
- use `land_connections_as_lines.shp` to inspect the connectivity network and update `rout_unit.con` as needed

---

## Reproducibility notes

- A zip backup of the SQLite DB is created at `Temp_/db_backup.zip`.
- Several files are backed up with `.bak` / `.bkp0` suffixes before being modified (e.g., `landuse.lum`, reservoir files, fertilizer/tillage tables).

---

