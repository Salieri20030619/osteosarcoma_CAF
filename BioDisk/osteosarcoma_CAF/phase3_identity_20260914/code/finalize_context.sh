#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n finalize_context.sh
# RUN_COMMAND: bash finalize_context.sh
# EXPECTED_SUCCESS_MARKER: PHASE3_CONTEXT_FINALIZED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
SOURCE=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/OS_CAF_phase2
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
echo $$ > "$ROOT/logs/finalize_context.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$ROOT/logs/finalize_context.exit_code"' EXIT
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 MPLBACKEND=Agg
"$PY" -m py_compile "$ROOT/code/integrate_identity.py" "$ROOT/code/validate_identity.py" "$ROOT/code/plot_cnv_summary.py" "$ROOT/code/build_delivery.py"
cp "$ROOT/identity/sample_identity_summary.tsv" "$ROOT/identity/sample_identity_summary_before_CTA_context.tsv"
if [[ -f "$ROOT/review/CTA_cluster_context_audit.tsv" ]];then cp "$ROOT/review/CTA_cluster_context_audit.tsv" "$ROOT/review/CTA_cluster_context_before_refinement.tsv";fi
if [[ -f "$ROOT/identity/COMPLETE" ]];then mv "$ROOT/identity/COMPLETE" "$ROOT/identity/COMPLETE_before_CTA_context";fi
"$PY" "$ROOT/code/integrate_identity.py" --root "$SOURCE" --phase3 "$ROOT" > "$ROOT/logs/identity_context_final.log" 2>&1
"$PY" "$ROOT/code/validate_identity.py" --root "$SOURCE" --phase3 "$ROOT" > "$ROOT/logs/identity_validation_context.log" 2>&1
"$PY" "$ROOT/code/audit_cta_context.py" "$SOURCE" "$ROOT" > "$ROOT/logs/CTA_context_final.log" 2>&1
# Reuse the unchanged actual CNV matrices and saved plot-data exports; only layout changes here.
"$PY" "$ROOT/code/plot_cnv_summary.py" "$ROOT" > "$ROOT/logs/cnv_layout_final.log" 2>&1
"$PY" "$ROOT/code/plot_caf_sets.py" "$SOURCE" "$ROOT" > "$ROOT/logs/caf_sets_context.log" 2>&1
"$PY" "$ROOT/code/redraw_marker_plots.py" "$ROOT" > "$ROOT/logs/marker_context_layout.log" 2>&1
"$PY" "$ROOT/code/write_figure_captions.py" "$ROOT" > "$ROOT/logs/figure_captions_final.log" 2>&1
"$PY" "$ROOT/code/build_delivery.py" --root "$SOURCE" --phase3 "$ROOT"
echo PHASE3_CONTEXT_FINALIZED
