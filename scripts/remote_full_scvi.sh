#!/usr/bin/env bash
set -euo pipefail

REMOTE="${REMOTE:-sabharwaly2@mia.ninds.nih.gov}"
REMOTE_DIR="${REMOTE_DIR:-~/olfactory}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
SESSION="${SESSION:-ora_full_scvi}"
SSH_OPTS="${SSH_OPTS:-}"
RSYNC_PROGRESS="${RSYNC_PROGRESS:---progress --partial}"
ACTION="${1:-help}"

usage() {
  cat <<'EOF'
Usage: scripts/remote_full_scvi.sh <check|push|launch|status|fetch>

Environment overrides:
  REMOTE=sabharwaly2@mia.ninds.nih.gov
  REMOTE_DIR=~/olfactory
  PYTHON_BIN=python3
  SESSION=ora_full_scvi
  SSH_OPTS="-o ControlPath=/Users/sabharwaly2/.ssh/cm/%r@%h:%p"

Actions:
  check   Test SSH and print host/python/GPU/memory summary.
  push    Sync repo source plus required raw Gateway H5AD to REMOTE_DIR.
  launch  Install environment if needed and start reduced full 4M scVI in tmux/nohup.
  status  Show running process and tail remote log.
  fetch   Sync reduced full 4M outputs and logs back to this workstation.
EOF
}

ssh_remote() {
  ssh ${SSH_OPTS} "${REMOTE}" "$@"
}

rsync_to_remote() {
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" \
    --exclude ".git/" \
    --exclude ".venv/" \
    --exclude ".mamba/" \
    --exclude "__pycache__/" \
    --exclude ".ruff_cache/" \
    --exclude "data/external/" \
    --exclude "data/processed/" \
    --exclude "results/" \
    ./ "${REMOTE}:${REMOTE_DIR}/"
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" data/raw/gateway.h5ad "${REMOTE}:${REMOTE_DIR}/data/raw/gateway.h5ad"
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" data/processed/gateway_scvi_stratified_250k.h5ad "${REMOTE}:${REMOTE_DIR}/data/processed/gateway_scvi_stratified_250k.h5ad"
}

remote_launch_command() {
  cat <<'EOF'
set -euo pipefail
cd "$REMOTE_DIR"
mkdir -p results/logs data/processed results/models results/tables
if [ ! -x .venv/bin/python ]; then
  "$PYTHON_BIN" -m venv .venv
fi
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e ".[dev,latent]"
PYTHON=.venv/bin/python make scvi-reduced-4m
PYTHON=.venv/bin/python make scvi-full-4m-reduced
PYTHON=.venv/bin/python make scvi-full-validation
PYTHON=.venv/bin/python make output-provenance
EOF
}

check_remote() {
  ssh_remote "hostname; uname -a; command -v ${PYTHON_BIN}; ${PYTHON_BIN} --version; free -h 2>/dev/null || true; nvidia-smi 2>/dev/null || true; df -h ."
}

push_remote() {
  ssh_remote "mkdir -p ${REMOTE_DIR}/data/raw"
  rsync_to_remote
}

launch_remote() {
  local tmp_script
  tmp_script="$(mktemp)"
  remote_launch_command > "${tmp_script}"
  ssh_remote "mkdir -p ${REMOTE_DIR}/results/logs"
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${tmp_script}" "${REMOTE}:${REMOTE_DIR}/results/logs/run_full_scvi_job.sh"
  rm -f "${tmp_script}"
  ssh_remote "chmod +x ${REMOTE_DIR}/results/logs/run_full_scvi_job.sh"
  if ssh_remote "command -v tmux >/dev/null 2>&1"; then
    ssh_remote "cd ${REMOTE_DIR} && REMOTE_DIR=${REMOTE_DIR} PYTHON_BIN=${PYTHON_BIN} tmux new-session -d -s ${SESSION} 'bash results/logs/run_full_scvi_job.sh > results/logs/scvi_full_4m.log 2>&1'"
    echo "Launched tmux session ${SESSION} on ${REMOTE}."
  else
    ssh_remote "cd ${REMOTE_DIR} && REMOTE_DIR=${REMOTE_DIR} PYTHON_BIN=${PYTHON_BIN} nohup bash results/logs/run_full_scvi_job.sh > results/logs/scvi_full_4m.log 2>&1 &"
    echo "Launched nohup job on ${REMOTE}."
  fi
}

status_remote() {
  ssh_remote "cd ${REMOTE_DIR} && { tmux ls 2>/dev/null || true; } && { pgrep -af 'build_reduced_h5ad|run_scvi_latent|make scvi-(reduced|full)' || true; } && tail -n 80 results/logs/scvi_full_4m.log 2>/dev/null || true"
}

fetch_remote() {
  mkdir -p data/processed results/models results/tables results/reports results/logs
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/data/processed/gateway_hvg3003_4m.h5ad" data/processed/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/data/processed/gateway_scvi_full_4m_reduced.h5ad" data/processed/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/data/processed/gateway_scvi_stratified_250k_seed23.h5ad" data/processed/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/models/gateway_scvi_full_4m_reduced/" results/models/gateway_scvi_full_4m_reduced/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/models/gateway_scvi_stratified_250k_seed23/" results/models/gateway_scvi_stratified_250k_seed23/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/tables/scvi_full_4m_reduced_validation.tsv" results/tables/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/tables/scvi_scaled_250k_seed23_validation.tsv" results/tables/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/reports/output_provenance.tsv" results/reports/ || true
  rsync -az ${RSYNC_PROGRESS} -e "ssh ${SSH_OPTS}" "${REMOTE}:${REMOTE_DIR}/results/logs/scvi_full_4m.log" results/logs/ || true
}

case "${ACTION}" in
  check) check_remote ;;
  push) push_remote ;;
  launch) launch_remote ;;
  status) status_remote ;;
  fetch) fetch_remote ;;
  help|-h|--help) usage ;;
  *) usage; exit 2 ;;
esac
