# teae_outputs.R
# Safety summary for study CDISCPILOT01 (public pharmaverse pilot data):
#   1. TABLE   - subjects with treatment-emergent adverse events (TEAEs) by preferred term
#   2. LIST    - listing of severe or serious TEAEs
#   3. FIGURE  - Kaplan-Meier curve of time to first TEAE
#
# Run from this folder:  Rscript teae_outputs.R   (or open in Positron and Run Source)
# Input : data/ADSL.csv (one row per subject), data/ADAE.csv (one row per TEAE)
# Output: output/teae_table.csv, output/severe_serious_listing.csv, output/km_time_to_first_teae.png

library(readr)
library(dplyr)
library(tidyr)
library(ggsurvfit)

# Protocol order
arms <- c("Placebo", "Xanomeline Low Dose", "Xanomeline High Dose")

# ---- Read the data -----------------------------------------------------------
# Both files are already limited to the safety population (SAFFL = "Y");
# ADAE contains treatment-emergent events only (TRTEMFL = "Y").
adsl <- read_csv("data/ADSL.csv", show_col_types = FALSE) |>
  mutate(
    TRT01A = factor(TRT01A, levels = arms),
    TRTSDT = as.Date(TRTSDT),
    TRTEDT = as.Date(TRTEDT),
    LSTALVDT = as.Date(LSTALVDT)
  )

adae <- read_csv("data/ADAE.csv", show_col_types = FALSE) |>
  mutate(
    TRT01A = factor(TRT01A, levels = arms),
    ASTDT = as.Date(ASTDT)
  )

dir.create("output", showWarnings = FALSE)

# ---- 1. TABLE: TEAEs by preferred term, n (%) of subjects --------------------
# Denominator for each arm
arm_n <- count(adsl, TRT01A, name = "N")

top10 <- adae |>
  # A subject counts once per term
  distinct(USUBJID, AEDECOD) |>
  count(AEDECOD, sort = TRUE) |>
  slice_head(n = 10) |>
  pull(AEDECOD)

teae_table <- adae |>
  filter(AEDECOD %in% top10) |>
  distinct(USUBJID, TRT01A, AEDECOD) |>
  count(AEDECOD, TRT01A, name = "n") |>
  # Show 0, never blank
  complete(AEDECOD, TRT01A, fill = list(n = 0L)) |>
  left_join(arm_n, by = "TRT01A") |>
  mutate(cell = sprintf("%d (%.1f%%)", n, 100 * n / N)) |>
  select(AEDECOD, TRT01A, cell) |>
  pivot_wider(names_from = TRT01A, values_from = cell) |>
  mutate(AEDECOD = factor(AEDECOD, levels = top10)) |>
  arrange(AEDECOD) |>
  rename(`Preferred term` = AEDECOD)

teae_table
write_csv(teae_table, "output/teae_table.csv")

# ---- 2. LIST: severe or serious TEAEs ----------------------------------------
severe_listing <- adae |>
  filter(AESEV == "SEVERE" | AESER == "Y") |>
  transmute(
    `Subject ID` = USUBJID,
    Arm = TRT01A,
    `Preferred term` = AEDECOD,
    Onset = ASTDT,
    Severity = AESEV,
    Serious = AESER,
    Relationship = AEREL
  ) |>
  arrange(`Subject ID`, Onset)

head(severe_listing, 10)
write_csv(severe_listing, "output/severe_serious_listing.csv")

# ---- 3. FIGURE: time to first TEAE (Kaplan-Meier) ----------------------------
first_teae <- adae |>
  group_by(USUBJID) |>
  summarise(teae_date = min(ASTDT), .groups = "drop")

tte <- adsl |>
  left_join(first_teae, by = "USUBJID") |>
  mutate(
    # 1 = had a TEAE
    event = as.integer(!is.na(teae_date)),
    # Censor 30 days after last dose
    end_date = if_else(!is.na(TRTEDT), TRTEDT + 30, LSTALVDT),
    time = pmax(
      1,
      as.numeric(if_else(event == 1, teae_date, end_date) - TRTSDT) + 1
    )
  )

km_plot <- survfit2(Surv(time, event) ~ TRT01A, data = tte) |>
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  labs(
    x = "Days since first dose",
    y = "Probability of no TEAE",
    title = "Time to First Treatment-Emergent Adverse Event by Treatment Arm"
  )

ggsave(
  "output/km_time_to_first_teae.png",
  km_plot,
  width = 8,
  height = 7,
  dpi = 150
)
message(
  "Saved: output/teae_table.csv, output/severe_serious_listing.csv, ",
  "output/km_time_to_first_teae.png"
)
