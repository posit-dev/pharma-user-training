# User-training example content

This repository contains a collection of content that facilitates learning how to use Posit Connect, Posit Package Manager, and Posit Workbench with RStudio, JupyterLab, and VS Code.

Each example produces the same recognizable pharma outputs: 

- a TEAE frequency **table**
- a **listing** of severe or serious TEAES
- a Kaplan-Meier **figure** of time to first TEAE


| # | Content | File | Run it |
|---|---------|------|--------|
| 1 | R script | `r_script.R` |  |
| 2 | Quarto report | `quarto_report.qmd` |  |
| 3 | Dashboard (Quarto) | `quarto_dashboard.qmd` |  |
| 4 | Shiny app | `app.R` |  |
| 5 | Python script | `py_script.py` |  |
| 6 | Jupyter notebook | `py_notebook_.ipynb` |  |
| 7 | requirements.txt | `requirements.txt` | `pip install -r requirements.txt` |

R packages: dplyr, tidyr, ggsurvfit (plus gt for the report and dashboard, shiny for the app).

Expected result in every example: 254 subjects, 217 with a TEAE; PRURITUS is the most frequent term (8, 21, and 25 subjects on placebo, low dose, and high dose); the listing has 42 rows.