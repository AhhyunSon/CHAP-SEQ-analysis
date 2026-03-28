#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, json, re
from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# sample parsing
# ============================================================
COND_PAT = re.compile(r"(Empty|GroEL|Seq576)", re.IGNORECASE)
BIN_PAT  = re.compile(r"(High|Low)", re.IGNORECASE)
REP_PAT  = re.compile(r"(\d+)")

def parse_sample(sample: str):
    s = str(sample)
    rep = int(REP_PAT.search(s).group(1)) if REP_PAT.search(s) else 1

    mcond = COND_PAT.search(s)
    cond = mcond.group(1).title() if mcond else "Empty"
    if cond.lower() == "groel": cond = "GroEL"
    if cond.lower() == "seq576": cond = "Seq576"
    if cond.lower() == "empty": cond = "Empty"

    mbin = BIN_PAT.search(s)
    binlab = mbin.group(1).title() if mbin else "Low"

    return cond, binlab, rep

# ============================================================
# WT geometry
# ============================================================
def build_residue_geom(wt_seq, cds_start, target="GCG"):
    rows = []
    n_res = len(wt_seq) // 3
    for r in range(1, n_res + 1):
        i = (r - 1) * 3
        wt = wt_seq[i:i+3]
        pos = [cds_start + i + j for j in range(3)]
        need_idx = [j for j in range(3) if wt[j] != target[j]]
        rows.append({"residue": r, "pos": pos, "need_idx": need_idx})
    return rows

# ============================================================
# STRICT raw engine
# ============================================================
def compute_strict_raw(df_sample, geom, target="GCG", alt_noise_floor=0.0):
    df_sample = df_sample.copy()
    df_sample["POS"] = pd.to_numeric(df_sample["POS"], errors="coerce")

    pos_groups = {}
    for pos_val, sub in df_sample.groupby("POS"):
        if pd.isna(pos_val):
            continue
        pos_groups[int(pos_val)] = sub

    records = []

    for g in geom:
        res = g["residue"]
        pos_list = g["pos"]
        need_idx = g["need_idx"]

        if len(need_idx) == 0:
            continue

        AD_raw_all = 0.0
        AD_raw_gcg = 0.0
        ad_vals = []
        dp_vals = []

        missing_required = False
        gcg_ok = True

        # ---- iterate ONLY required positions ----
        for j in need_idx:
            pos = pos_list[j]

            # position missing entirely
            if pos not in pos_groups:
                missing_required = True
                break

            sub = pos_groups[pos]
            DP_j = int(sub["DP"].iloc[0])

            if DP_j <= 0:
                missing_required = True
                break

            # raw AD regardless of GCG
            AD_all_j = int(sub["AD"].sum())
            AD_raw_all += AD_all_j

            # check GCG evidence
            hit = sub[sub["ALT"].astype(str).str.upper() == target[j]]
            if len(hit) == 0:
                gcg_ok = False
                break

            AD_gcg_j = int(hit["AD"].sum())
            AD_raw_gcg += AD_gcg_j
            ad_vals.append(AD_gcg_j)
            dp_vals.append(DP_j)

        if missing_required:
            records.append({
                "residue": res,
                "AD_raw_all": AD_raw_all,
                "AD_raw_gcg": 0.0,
                "AD": 0.0,
                "DP": 0.0,
                "freq_strict": np.nan,
                "flag": "missing_position"
            })
            continue

        if not gcg_ok:
            records.append({
                "residue": res,
                "AD_raw_all": AD_raw_all,
                "AD_raw_gcg": 0.0,
                "AD": 0.0,
                "DP": 0.0,
                "freq_strict": np.nan,
                "flag": "no_gcg_evidence"
            })
            continue

        # ---- valid alanine codon (n = len(need_idx)) ----
        AD = float(np.mean(ad_vals))
        DP = float(np.mean(dp_vals))
        freq_strict = AD / DP if DP > 0 else np.nan

        records.append({
            "residue": res,
            "AD_raw_all": AD_raw_all,
            "AD_raw_gcg": AD_raw_gcg,
            "AD": AD,
            "DP": DP,
            "freq_strict": freq_strict,
            "flag": "ok"
        })

    return pd.DataFrame(records)



# ============================================================
# MAIN
# ============================================================
def main():
    ap = argparse.ArgumentParser("Strict alanine raw-only pipeline")
    ap.add_argument("--merged_csv", default="merged_variants_raw.csv")
    ap.add_argument("--config", default="config.json")
    ap.add_argument("--target", default="GCG")
    ap.add_argument("--alt_noise_floor", type=float, default=0.0)
    args = ap.parse_args()


    cfg = json.load(open(args.config))
    cds_start = int(cfg["cds_start"])
    wt_seq = cfg["wt_cds_seq"].upper()

    df = pd.read_csv("merged_variants_raw.csv")
    geom = build_residue_geom(wt_seq, cds_start, args.target)

    all_rows = []

    for sample, sub in df.groupby("sample"):
        cond, binlab, rep = parse_sample(sample)
        strict = compute_strict_raw(sub, geom, args.target, args.alt_noise_floor)
        strict["sample"] = sample
        strict["condition"] = cond
        strict["bin"] = binlab
        strict["replicate"] = rep
        all_rows.append(strict)

    out = pd.concat(all_rows, ignore_index=True)
    out.to_csv("per_sample_strict_raw.csv", index=False)

    print("✅ DONE")
    print("Output:", "per_sample_strict_raw.csv")

if __name__ == "__main__":
    main()

