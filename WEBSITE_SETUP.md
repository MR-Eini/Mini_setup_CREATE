# Website setup instructions for `Mini_setup_CREATE`

This package adds a MkDocs Material documentation website to:

```text
https://github.com/MR-Eini/Mini_setup_CREATE
```

The expected public website URL after deployment is:

```text
https://mr-eini.github.io/Mini_setup_CREATE/
```

## 1. Files included in this package

Copy these files/folders into the **root** of the repository:

```text
mkdocs.yml
requirements-docs.txt
site_docs/
overrides/
.github/workflows/deploy-docs.yml
WEBSITE_SETUP.md
```

The package intentionally uses `site_docs/` instead of `docs/` because the repository already has a `Docs/` folder. This avoids case-sensitivity problems on Windows.

## 2. Optional `.gitignore` additions

Append these lines to the existing repository `.gitignore`:

```gitignore
# MkDocs local build output
site/
.cache/
.venv/
```

## 3. Preview locally

From the repository root:

### Windows PowerShell

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements-docs.txt
mkdocs serve
```

### macOS / Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements-docs.txt
mkdocs serve
```

Then open the URL printed by MkDocs. It is commonly:

```text
http://127.0.0.1:8000/Mini_setup_CREATE/
```

## 4. Validate the build

Run:

```bash
mkdocs build --strict
```

This should complete without warnings. Fix any broken links before pushing.

## 5. Commit and push

```bash
git add mkdocs.yml requirements-docs.txt site_docs .github/workflows/deploy-docs.yml WEBSITE_SETUP.md
git commit -m "Add MkDocs documentation website"
git push origin main
```

## 6. Enable GitHub Pages

On GitHub:

1. Open the repository.
2. Go to **Settings**.
3. Open **Pages**.
4. Under **Build and deployment**, set **Source** to **GitHub Actions**.
5. Open the **Actions** tab and run or inspect **Deploy documentation website**.
6. After deployment succeeds, open:

```text
https://mr-eini.github.io/Mini_setup_CREATE/
```

## 7. Update repository About section

On the GitHub repository main page, click the gear icon next to **About** and set:

```text
Website: https://mr-eini.github.io/Mini_setup_CREATE/
Description: SWAT+ setup and scenario preparation toolkit for OPTAIN-style workflows.
Topics: swatplus, swat, hydrology, calibration, water-quality, mkdocs, r
```

## 8. Future maintenance

When you edit any Markdown page in `site_docs/` or edit `mkdocs.yml`, push to `main`. The GitHub Actions workflow will rebuild and redeploy the website automatically.

Do not manually edit the generated `site/` folder.


## Visual assets added in the updated package

The updated website package copies the following assets into `site_docs/assets/`:

```text
swativerse_update.png
Water4all_0.png
EN_co_fundedvertical_RGB_POS.png
logo_CREATE_horizontal.png
swatplus_logo.svg
```

The first four files are taken from the repository `assets/` folder. The SWAT+ SVG is a local display badge because the repository did not contain an official SWAT+ logo file. Replace it with an official file if you have permission to reuse one.

The file `overrides/main.html` adds a logo footer to every generated webpage. Do not delete the `overrides/` folder unless you also remove `custom_dir: overrides` from `mkdocs.yml`.
