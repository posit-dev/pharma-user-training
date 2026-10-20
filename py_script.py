"""teae_outputs.py

Safety summary for study CDISCPILOT01 (public pharmaverse pilot data):
  1. TABLE  - subjects with treatment-emergent adverse events (TEAEs) by preferred term
  2. LIST   - listing of severe or serious TEAEs
  3. FIGURE - Kaplan-Meier curve of time to first TEAE

Run from this folder:  python teae_outputs.py
Input : data/ADSL.csv (one row per subject), data/ADAE.csv (one row per TEAE)
Output: output/teae_table.csv, output/severe_serious_listing.csv,
        output/km_time_to_first_teae.png
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter

ARMS = ["Placebo", "Xanomeline Low Dose", "Xanomeline High Dose"]  # protocol order

# ---- Read the data -----------------------------------------------------------
# Both files are already limited to the safety population (SAFFL = "Y");
# ADAE contains treatment-emergent events only (TRTEMFL = "Y").
adsl = pd.read_csv("data/ADSL.csv", parse_dates=["TRTSDT", "TRTEDT", "LSTALVDT"])
adae = pd.read_csv("data/ADAE.csv", parse_dates=["ASTDT"])

Path("output").mkdir(exist_ok=True)

# ---- 1. TABLE: TEAEs by preferred term, n (%) of subjects --------------------
arm_n = adsl["TRT01A"].value_counts().reindex(ARMS)          # denominator for each arm

subject_term = adae[["USUBJID", "TRT01A", "AEDECOD"]].drop_duplicates()  # a subject counts once per term

# Ten terms reported by the most subjects (ties stay in alphabetical order)
term_counts = subject_term.groupby("AEDECOD")["USUBJID"].nunique()
top10 = term_counts.sort_values(ascending=False, kind="stable").head(10).index

counts = (
    subject_term[subject_term["AEDECOD"].isin(top10)]
    .groupby(["AEDECOD", "TRT01A"])
    .size()
    .unstack(fill_value=0)          # fill_value=0: show 0, never blank
    .reindex(index=top10, columns=ARMS, fill_value=0)
)
percent = 100 * counts / arm_n      # divide each column by its arm's N
teae_table = counts.astype(str) + " (" + percent.map("{:.1f}%".format) + ")"
teae_table.index.name = "Preferred term"

print(teae_table.to_string())
teae_table.to_csv("output/teae_table.csv")

# ---- 2. LIST: severe or serious TEAEs ----------------------------------------
severe_listing = (
    adae[(adae["AESEV"] == "SEVERE") | (adae["AESER"] == "Y")]
    .rename(columns={"USUBJID": "Subject ID", "TRT01A": "Arm", "AEDECOD": "Preferred term",
                     "ASTDT": "Onset", "AESEV": "Severity", "AESER": "Serious",
                     "AEREL": "Relationship"})
    [["Subject ID", "Arm", "Preferred term", "Onset", "Severity", "Serious", "Relationship"]]
    .sort_values(["Subject ID", "Onset"])
)

print(severe_listing.head(10).to_string(index=False))
severe_listing.to_csv("output/severe_serious_listing.csv", index=False)

# ---- 3. FIGURE: time to first TEAE (Kaplan-Meier) ----------------------------
first_teae = adae.groupby("USUBJID", as_index=False)["ASTDT"].min().rename(columns={"ASTDT": "teae_date"})

tte = adsl.merge(first_teae, on="USUBJID", how="left")
tte["event"] = tte["teae_date"].notna().astype(int)                  # 1 = had a TEAE
end_date = (tte["TRTEDT"] + pd.Timedelta(days=30)).fillna(tte["LSTALVDT"])  # censor 30 days after last dose
stop_date = tte["teae_date"].where(tte["event"] == 1, end_date)
tte["time"] = ((stop_date - tte["TRTSDT"]).dt.days + 1).clip(lower=1)

fig, ax = plt.subplots(figsize=(8, 5))
for arm in ARMS:
    arm_data = tte[tte["TRT01A"] == arm]
    KaplanMeierFitter().fit(arm_data["time"], arm_data["event"], label=arm).plot_survival_function(ax=ax)
ax.set_xlabel("Days since first dose")
ax.set_ylabel("Probability of no TEAE")
ax.set_title("Time to First Treatment-Emergent Adverse Event by Treatment Arm")
fig.tight_layout()
fig.savefig("output/km_time_to_first_teae.png", dpi=150)

print("\nSaved: output/teae_table.csv, output/severe_serious_listing.csv, output/km_time_to_first_teae.png")
