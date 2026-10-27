#!/usr/bin/env bash
set -euo pipefail
cd /mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase7_review_states_20260917
while kill -0 2835691 2>/dev/null; do sleep 900; done
/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python download_expansion_20260917/code/mark_pause.py
