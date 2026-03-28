import pandas as pd
import numpy as np
from scipy.special import erfinv

input_file = "step1_p_high_frequency_based_with_lowerF_small1.csv"
output_file = "step2_Fhigh_from_pHigh_small1.csv"

sigma_overall_med = 0.230437356
eps = 1e-6

def fmean_from_p(p, fg, sigma):
    try:
        return fg * np.exp(
            0.5 * sigma**2
            - sigma * np.sqrt(2) * erfinv(1 - 2 * p)
        )
    except Exception:
        return np.nan

df = pd.read_csv(input_file)

if "p_high_smoothed" not in df.columns:
    raise ValueError("p_high_smoothed not found")

df["p_high_raw_for_integral"] = pd.to_numeric(df["p_high_smoothed"], errors="coerce")
df["lowerF"] = pd.to_numeric(df["lowerF"], errors="coerce")

# unclipped
df["F_high_unclipped"] = df.apply(
    lambda r: fmean_from_p(r["p_high_raw_for_integral"], r["lowerF"], sigma_overall_med),
    axis=1
)

# clipped
df["p_high_clipped"] = df["p_high_raw_for_integral"].clip(lower=eps, upper=1 - eps)
df["F_high_clipped"] = df.apply(
    lambda r: fmean_from_p(r["p_high_clipped"], r["lowerF"], sigma_overall_med),
    axis=1
)

# flags
df["flag_unclipped_nan_or_inf"] = (
    df["F_high_unclipped"].isna()
    | (~np.isfinite(df["F_high_unclipped"]))
)

df["flag_clipping_changed_p"] = df["p_high_clipped"] != df["p_high_raw_for_integral"]

df["flag_clipping_changed_F"] = (
    np.isfinite(df["F_high_unclipped"])
    & np.isfinite(df["F_high_clipped"])
    & (df["F_high_unclipped"] != df["F_high_clipped"])
)
df["p_high_raw_for_integral"] = df["p_high_raw_for_integral"].fillna(0.5)

# choose final F_high: prefer unclipped unless it fails
df["F_high_final"] = np.where(
    df["flag_unclipped_nan_or_inf"],
    df["F_high_clipped"],
    df["F_high_unclipped"]
)

# optional: remove residue 1
df = df[df["residue"] != 1].copy()

df.to_csv(output_file, index=False)

print("Wrote:", output_file)
print("Unclipped failed rows:", int(df["flag_unclipped_nan_or_inf"].sum()))
print("Clipping changed p rows:", int(df["flag_clipping_changed_p"].sum()))
