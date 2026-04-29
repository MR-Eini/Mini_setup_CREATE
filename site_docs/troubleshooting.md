# Troubleshooting

## Setup workflow

| Problem | Likely cause | Action |
|---|---|---|
| More than one project SQLite database is detected | Old generated files are still present in the result folder. | Clean the result folder before rerunning. |
| `write.exe` fails | SQLite is open, database path is wrong, or `project_config` is incompatible. | Close SWAT+ Editor, DB Browser, and Explorer preview; check `input_files_dir`. |
| SWAT+ executable does not run | Wrong executable name or missing executable in `Libraries`. | Check `swat_exe` in `settings.R`. |
| Connectivity looks wrong | Input topology, CRS, routing, or channel/land geometry problem. | Inspect the connectivity shapefile generated in Step 7. |
| `landuse.lum` pointers are wrong | Prefix rules do not match current land-use names. | Edit `read_and_modify_landuse_lum.R`. |
| Reservoir files behave unexpectedly | Template fixes are project-specific. | Review `reservoir.con`, `reservoir.res`, and `hydrology.res`. |

## SWATdoctR workflow

| Problem | Likely cause | Action |
|---|---|---|
| `plots.pdf` is incomplete | Error occurred before `dev.off()`. | Run `dev.off()` manually and rerun the failed plot. |
| Climate plots show missing or unrealistic values | Weather input formatting or unit problem. | Check weather files and SWAT+ weather assignment. |
| Snow plots look implausible | Temperature threshold, elevation, or snow parameter issue. | Check climate inputs and snow settings. |
| Yield is nearly zero | Management, crop parameters, stress, or harvest operation issue. | Compare stress/no-stress runs and PHU plots. |

## Crop-yield soft calibration

| Problem | Likely cause | Action |
|---|---|---|
| `dmat_sel` not found | R session was restarted after `01_crop_phu.R`. | Save and reload `dmat_sel`, or rerun the first script. |
| Yield improves but PHU is implausible | Parameters are fitting yield only. | Recheck crop maturity and harvest timing. |
| Water-yield ratio improves but crop yield worsens | `esco/epco` values affect plant water use. | Use the final crop-yield check before writing `hydrology.hyd`. |

## River-discharge hard calibration

| Problem | Likely cause | Action |
|---|---|---|
| Observation file not found | Script references case-study-specific CSV names. | Update paths or add the required observation CSV files. |
| Expected output object is missing | `cha_ids` in `02_define_output.R` does not match the analysis script. | Align output definitions with gauge IDs. |
| 900-run batch is too slow | Too many parameter combinations for workshop/demo. | Test first with 5–10 runs. |
| Metrics are not comparable across stations | Unit or alignment inconsistency. | Verify units, dates, and alignment mode before ranking. |

## Website build

| Problem | Cause | Action |
|---|---|---|
| `mkdocs build --strict` fails | Broken internal link or missing file. | Read the error line and fix the link or file path. |
| GitHub Pages does not deploy | Pages source is not set to GitHub Actions, or workflow failed. | Go to **Settings → Pages → Build and deployment → Source → GitHub Actions**. |
| The site builds locally but not on GitHub | Dependency or case-sensitive path issue. | Check the Actions log and file capitalization. |

