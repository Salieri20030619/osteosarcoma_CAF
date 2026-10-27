#!/usr/bin/env bash
set -eo pipefail
CODE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$CODE_DIR/.." && pwd)"
MODE="${1:-start}"
mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/runs" "$PROJECT_DIR/raw/GSE270231" "$PROJECT_DIR/exports"
prepare() {
  bash -n "$CODE_DIR/run_server.sh"
  if ! command -v conda >/dev/null 2>&1; then
    for candidate in "$HOME/miniconda3/etc/profile.d/conda.sh" "$HOME/anaconda3/etc/profile.d/conda.sh"; do
      if [[ -f "$candidate" ]]; then source "$candidate"; break; fi
    done
  fi
  command -v conda >/dev/null 2>&1 || { echo 'ERROR: conda not found; open your normal conda shell.'; return 1; }
  source "$(conda info --base)/etc/profile.d/conda.sh"
  if conda run -n batch_bench python -c 'import sys' >/dev/null 2>&1; then
    conda activate batch_bench
  else
    conda create -y -n batch_bench python=3.12
    conda activate batch_bench
  fi
  # Keep pinned analysis dependencies inside the project, preserving batch_bench packages.
  if [[ ! -x "$PROJECT_DIR/runtime/bin/python" ]]; then
    if python -c 'import sys; assert sys.version_info[:2] == (3,12)' >/dev/null 2>&1; then
      python -m venv "$PROJECT_DIR/runtime"
    else
      conda create -y -p "$PROJECT_DIR/runtime" python=3.12 pip
    fi
  fi
  "$PROJECT_DIR/runtime/bin/python" -m py_compile "$CODE_DIR/phase2.py" "$CODE_DIR/pack_results.py" "$CODE_DIR/finalize_review.py"
  echo 'SYNTAX_CHECK=PASS'
  "$PROJECT_DIR/runtime/bin/python" -m pip install -r "$CODE_DIR/requirements.txt"
  "$PROJECT_DIR/runtime/bin/python" -c 'import scanpy,skimage,igraph,leidenalg,anndata; print("DEPENDENCIES_PASS")'
}
case "$MODE" in
  start)
    if [[ -f "$PROJECT_DIR/logs/worker.pid" ]]; then
      read -r OLD_PID < "$PROJECT_DIR/logs/worker.pid"
      if [[ "$OLD_PID" =~ ^[0-9]+$ ]] && kill -0 "$OLD_PID" 2>/dev/null; then
        echo "A run is active: PID=$OLD_PID"; exit 1
      fi
    fi
    prepare
    RUN_DIR="$(mktemp -d "$PROJECT_DIR/runs/run_$(date +%Y%m%d_%H%M%S)_XXXXXX")"
    printf '%s\n' "$RUN_DIR" > "$PROJECT_DIR/logs/latest_run.txt"
    nohup bash "$CODE_DIR/run_server.sh" worker "$RUN_DIR" > "$PROJECT_DIR/logs/phase2_latest.log" 2>&1 < /dev/null &
    PID=$!
    printf '%s\n' "$PID" > "$PROJECT_DIR/logs/worker.pid"
    echo "RUN_DIR=$RUN_DIR"
    echo "PID=$PID"
    echo 'RUN_COMMAND=bash code/run_server.sh start'
    echo 'EXPECTED_SUCCESS_MARKER=PHASE2_SERVER_COMPLETE'
    ;;
  worker)
    set -u
    RUN_DIR="${2:?Missing run directory}"
    trap 'rc=$?; if (( rc != 0 )); then echo "PHASE2_FAILED exit=$rc"; printf "%s\n" "$rc" > "$RUN_DIR/FAILED.txt"; fi' EXIT
    export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 NUMBA_NUM_THREADS=4 MPLBACKEND=Agg
    "$PROJECT_DIR/runtime/bin/python" "$CODE_DIR/phase2.py" --manifest "$CODE_DIR/manifest.tsv" --raw "$PROJECT_DIR/raw/GSE270231" --out "$RUN_DIR"
    "$PROJECT_DIR/runtime/bin/python" "$CODE_DIR/finalize_review.py" --results "$RUN_DIR"
    "$PROJECT_DIR/runtime/bin/python" "$CODE_DIR/pack_results.py" --results "$RUN_DIR" --output "$PROJECT_DIR/exports/phase2_review_latest.zip"
    echo 'PHASE2_SERVER_COMPLETE'
    ;;
  monitor)
    touch "$PROJECT_DIR/logs/phase2_latest.log"
    tail -n 30 -F "$PROJECT_DIR/logs/phase2_latest.log"
    ;;
  status)
    [[ ! -f "$PROJECT_DIR/logs/latest_run.txt" ]] || cat "$PROJECT_DIR/logs/latest_run.txt"
    [[ ! -f "$PROJECT_DIR/logs/worker.pid" ]] || ps -p "$(cat "$PROJECT_DIR/logs/worker.pid")" -o pid,etime,%cpu,%mem,args || true
    [[ ! -f "$PROJECT_DIR/logs/phase2_latest.log" ]] || tail -n 15 "$PROJECT_DIR/logs/phase2_latest.log"
    ;;
  *) echo 'Usage: bash code/run_server.sh {start|monitor|status}'; exit 2 ;;
esac
