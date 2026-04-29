# Getting started

## Clone the repository

```bash
git clone https://github.com/MR-Eini/Mini_setup_CREATE.git
cd Mini_setup_CREATE
```

## Workflow sequence

| Step | Folder | Purpose |
|---:|---|---|
| 1 | `1_Setup/` | Generate a complete SWAT+ setup from prepared inputs. |
| 2 | `2_SWATdoctR/` | Verify the generated clean setup with diagnostic runs and plots. |
| 3 | `3_CropYield/softcal/` | Soft-calibrate crop development, yield, and water-yield parameters. |
| 4 | `4_RiverDischarge/hardcal/` | Hard-calibrate river discharge and optionally water-quality variables. |
| 5 | `5_NBS/` | Optional nature-based solution scenario workflows. |

## Core requirements

| Requirement | Notes |
|---|---|
| R 4.x | Use one fixed R version for all steps. |
| RStudio | Recommended for interactive script execution and manual parameter selection. |
| SWAT+ executable | Must match the model files and workflow configuration. |
| `write.exe` | Used in the setup workflow to generate SWAT+ text input files from the SQLite database. |
| WhiteboxTools | Used by SWATbuildR-related terrain and connectivity processing. |
| SWAT-related R packages | Typical packages include `SWATprepR`, `SWATfarmR`, `SWATtunR`, `SWATrunR`, and `SWATdoctR`. |

## Important execution principle

The repository is designed as a staged workflow, not as one fully automatic command.

- `1_Setup` assumes that input GIS, soil, land-use, weather, crop, and management data are already prepared.
- `2_SWATdoctR` assumes that a runnable `clean_setup` already exists.
- `3_CropYield/softcal` depends on manual interpretation of plots.
- `4_RiverDischarge/hardcal` requires observation files and station/output consistency checks.

