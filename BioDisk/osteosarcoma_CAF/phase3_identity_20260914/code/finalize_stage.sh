#!/usr/bin/env bash
# SYNTAX_CHECK: bash -n finalize_stage.sh
# RUN_COMMAND: bash finalize_stage.sh
# EXPECTED_SUCCESS_MARKER: PHASE3_STAGE_FINALIZED
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/phase3_identity_20260914
SOURCE=/mnt/BioDisk/wangmingyuan/osteosarcoma_CAF/OS_CAF_phase2
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RBIN=/mnt/BioDisk/wangmingyuan/project26210_envs/phase11_infercnv_conda_v2/bin/Rscript
echo $$ > "$ROOT/logs/finalize_stage.pid"
trap 'rc=$?; printf "%s\n" "$rc" > "$ROOT/logs/finalize_stage.exit_code"' EXIT
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 MPLBACKEND=Agg
for method in infercnv copykat; do
 for s in S0058 S0059 S0217 S0218 SC069 SC072 SC073; do test -f "$ROOT/$method/$s/COMPLETE";done
done
for method in copykat_balanced scevan;do for s in S0058 SC072;do test -f "$ROOT/$method/$s/COMPLETE";done;done
"$PY" -m py_compile "$ROOT/code/integrate_identity.py" "$ROOT/code/validate_identity.py" "$ROOT/code/build_delivery.py" "$ROOT/code/plot_cnv_summary.py" "$ROOT/code/plot_caf_sets.py"
"$RBIN" --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$ROOT/code/export_cnv_plot_data.R"
if [[ -f "$ROOT/identity/COMPLETE" ]];then mv "$ROOT/identity/COMPLETE" "$ROOT/identity/COMPLETE_pre_metadata_update";fi
cp "$ROOT/identity/sample_identity_summary.tsv" "$ROOT/identity/sample_identity_summary_before_metadata.tsv"
/usr/bin/time -v -o "$ROOT/logs/identity_final_resources.txt" "$PY" "$ROOT/code/integrate_identity.py" --root "$SOURCE" --phase3 "$ROOT" > "$ROOT/logs/identity_final.log" 2>&1
# Keep the prior summary for audit; biological context refinements may deliberately change identities.
"$PY" "$ROOT/code/validate_identity.py" --root "$SOURCE" --phase3 "$ROOT" > "$ROOT/logs/identity_validation.log" 2>&1
"$RBIN" --vanilla "$ROOT/code/export_cnv_plot_data.R" > "$ROOT/logs/cnv_plot_export.log" 2>&1
"$PY" "$ROOT/code/plot_cnv_summary.py" "$ROOT" > "$ROOT/logs/cnv_plot_summary.log" 2>&1
"$PY" "$ROOT/code/plot_caf_sets.py" "$SOURCE" "$ROOT" > "$ROOT/logs/caf_set_plots.log" 2>&1
"$PY" -m py_compile "$ROOT/code/redraw_marker_plots.py"
"$PY" "$ROOT/code/redraw_marker_plots.py" "$ROOT" > "$ROOT/logs/marker_plot_layout.log" 2>&1
"$PY" -m py_compile "$ROOT/code/write_figure_captions.py"
"$PY" "$ROOT/code/write_figure_captions.py" "$ROOT" > "$ROOT/logs/figure_captions.log" 2>&1
"$PY" "$ROOT/code/build_delivery.py" --root "$SOURCE" --phase3 "$ROOT"
echo PHASE3_STAGE_FINALIZED
