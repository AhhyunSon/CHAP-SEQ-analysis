#!/bin/bash

python3 step0_alanine_strict_full_pipeline.py
python3 step1_pHigh.py
python3 step2_Fhigh.py
python3 step3_Metrics.py
