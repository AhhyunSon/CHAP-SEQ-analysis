# CHAP-seq analysis

This repository contains code for analyzing CHAP-seq data.

## Pipeline

1. step0: codon-level filtering  
2. step1: pHigh calculation and stabilization  
3. step2: Fmean inference from pHigh  
4. step3: residue-level metrics (Fmean, ratioF, log2ratioF)  

## Input files

- merged_variants_raw.csv
- FACS_Statistics_Round1.xlsx
- config.json

## Run

```bash
python step0_alanine_strict_full_pipeline.py
python step1_pHigh.py
python step2_FHigh.py
python step3_Metrics.py
