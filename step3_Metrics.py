import pandas as pd
import numpy as np

# ============================================================
# Load Data
# ============================================================
df = pd.read_csv("step2_Fhigh_from_pHigh_small1.csv")
df.columns = df.columns.str.strip()

# ============================================================
# Use correct columns
# ============================================================
df["F_high"] = df["F_high_final"]
df["p_high"] = df["p_high_smoothed"]

# ============================================================
# Residue-level averaging
# ============================================================
avg = (
    df.groupby(["condition", "residue"], as_index=False)
      .agg(Fmean=("F_high", "mean"))
)

# ============================================================
# Pivot
# ============================================================
Fmean_wide = avg.pivot(index="residue", columns="condition", values="Fmean")

Fmean_wide = Fmean_wide.rename(columns={
    "Empty": "Fmean_Empty",
    "GroEL": "Fmean_GroEL",
    "Seq576": "Fmean_Seq576"
})

res = Fmean_wide.copy()

# ============================================================
# Ratios
# ============================================================
eps = 1e-9

res["ratioF_GroEL"]  = res["Fmean_GroEL"]  / (res["Fmean_Empty"] + eps)
res["ratioF_Seq576"] = res["Fmean_Seq576"] / (res["Fmean_Empty"] + eps)

res["log2ratioF_GroEL"]  = np.log2(res["ratioF_GroEL"] + eps)
res["log2ratioF_Seq576"] = np.log2(res["ratioF_Seq576"] + eps)

# ============================================================
# Save
# ============================================================
res = res.reset_index()
res.to_csv("Step3_Residue_Folding_Metrics_small1.csv", index=False)

print("saved")
