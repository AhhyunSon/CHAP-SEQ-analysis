import pandas as pd
import numpy as np

# =========================
# INPUT
# =========================
PVT_CSV = "per_sample_strict_raw.csv"
FACS_XLSX = "FACS_Statistics_Round1.xlsx"
OUT_CSV = "step1_p_high_frequency_based_with_lowerF_small1.csv"

# =========================
# LOAD
# =========================
pvt = pd.read_csv(PVT_CSV)
fac = pd.read_excel(FACS_XLSX)

# =========================
# CLEAN TYPES (PVT)
# =========================
pvt["bin"] = pvt["bin"].astype(str).str.strip()
pvt["condition"] = pvt["condition"].astype(str).str.strip()
pvt["replicate"] = pd.to_numeric(pvt["replicate"], errors="coerce")
pvt["residue"] = pd.to_numeric(pvt["residue"], errors="coerce")

# =========================
# CLEAN TYPES (FACS)
# =========================
# NOTE:
# Name      -> condition
# Replicate -> replicate
# Gate == B-Q4, Minimum -> lowerF
fac["condition"] = fac["Name"].astype(str).str.strip()
fac["replicate"] = pd.to_numeric(fac["Replicate"], errors="coerce")
fac["Gate"] = fac["Gate"].astype(str).str.strip()

# =========================
# LONG -> WIDE (AD / DP)
# =========================
wide = (
    pvt.pivot_table(
        index=["condition", "replicate", "residue"],
        columns="bin",
        values=["AD", "DP"],
        aggfunc="sum",
        fill_value=0
    )
)

wide.columns = [f"{metric}_{b[0].upper()}" for metric, b in wide.columns]
wide = wide.reset_index()

# ensure expected columns exist
for col in ["AD_H", "AD_L", "DP_H", "DP_L"]:
    if col not in wide.columns:
        wide[col] = 0

# =========================
# RAW FREQUENCY
# =========================
wide["Freq_H"] = np.where(
    wide["DP_H"] > 0,
    wide["AD_H"] / wide["DP_H"],
    np.nan
)

wide["Freq_L"] = np.where(
    wide["DP_L"] > 0,
    wide["AD_L"] / wide["DP_L"],
    np.nan
)
# explicit zero fill for downstream comparison
wide["Freq_H"] = wide["Freq_H"].fillna(0.0)
wide["Freq_L"] = wide["Freq_L"].fillna(0.0)

# no bin signal at all (after processing)
wide["flag_no_bin_signal"] = (wide["Freq_H"] == 0) & (wide["Freq_L"] == 0)

# likely filtered by strict codon / QC
wide["flag_likely_filtered"] = wide["flag_no_bin_signal"] & (
    (wide["DP_H"] == 0) | (wide["DP_L"] == 0)
)

# =========================
# RAW p_high
# =========================
den_freq = wide["Freq_H"] + wide["Freq_L"]
wide["p_high_raw"] = np.where(
    den_freq > 0,
    wide["Freq_H"] / den_freq,
    np.nan
)

# =========================
# FLAGS (RAW)
# =========================
wide["flag_raw_extreme"] = wide["p_high_raw"].isin([0.0, 1.0])

wide["flag_zero_AD_or_DP"] = (
    (wide["AD_H"] == 0) | (wide["AD_L"] == 0) |
    (wide["DP_H"] == 0) | (wide["DP_L"] == 0)
)

# =========================
# p_high smoothing (round-level scale-adaptive stabilization)
# =========================

fH = wide["Freq_H"].fillna(0.0)
fL = wide["Freq_L"].fillna(0.0)
den = fH + fL

# Round-level small constant (5% quantile / 2)
nonzero_den = den[den > 0]

if len(nonzero_den) > 0:
    small = np.quantile(nonzero_den, 0.05) / 2.0
else:
    small = 1e-6  # fallback safety

wide["p_high_smoothed"] = np.where(
    den == 0,
    0.5,
    (fH + small) / (den + 2.0 * small)
)





# =========================
# FACS lowerF (replicate-level, B-Q4 Minimum ONLY)
# =========================
fac_q4 = fac[fac["Gate"] == "B-Q4"].copy()
fac_q4["lowerF"] = fac_q4["Minimum"]

fac_wide = (
    fac_q4[["condition", "replicate", "lowerF"]]
    .drop_duplicates()
)

# =========================
# MERGE FACS lowerF
# =========================
wide = wide.merge(
    fac_wide,
    on=["condition", "replicate"],
    how="left"
)

# =========================
# REMOVE RESIDUE 1
# =========================
wide = wide[wide["residue"] != 1].copy()

# =========================
# COLUMN ORDER
# =========================
cols_order = [
    "condition", "replicate", "residue",
    "AD_H", "DP_H", "Freq_H",
    "AD_L", "DP_L", "Freq_L",
    "lowerF",
    "p_high_raw", "p_high_smoothed",
    "flag_raw_extreme", "flag_zero_AD_or_DP",
]

wide = wide[cols_order]

# =========================
# WRITE
# =========================
wide.to_csv(OUT_CSV, index=False)

print("Final CSV written:")
print(" ", OUT_CSV)
wide.groupby(["condition","replicate"])["lowerF"].nunique()
