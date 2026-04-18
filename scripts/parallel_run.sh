#!/usr/bin/env bash
# Run a commands file in parallel using GNU parallel.
#
# Usage:
#   ./scripts/parallel_run.sh -c data/processed/group_fastas/tn93_commands.txt
#   ./scripts/parallel_run.sh -c <commands.txt> -j 16 --progress
#   ./scripts/parallel_run.sh -c <commands.txt> --retries 2 --resume-failed

set -euo pipefail

COMMANDS_FILE=""
JOBS=""
JOBLOG=""
RETRIES=0
TIMEOUT=""
PROGRESS=false
DRY_RUN=false
RESUME_FAILED=false
TMPDIR_ARG=""
NO_COMPRESS=false

cpu_count() {
  command -v nproc >/dev/null 2>&1 && nproc && return
  [[ "$OSTYPE" == darwin* ]] && sysctl -n hw.ncpu && return
  echo 1
}

print_usage() {
  cat <<'EOF'
Usage: scripts/parallel_run.sh -c COMMANDS_FILE [options]

Required:
  -c, --commands FILE    Path to commands file (one command per line)

Options:
  -j, --jobs N           Parallel workers (default: CPU count)
      --joblog FILE      GNU parallel joblog path (default: <commands>.joblog.tsv)
      --retries N        Retry failed jobs N times (default: 0)
      --timeout SECS     Per-job timeout in seconds
      --progress         Show live progress bar
      --resume-failed    Re-run only previously failed jobs from joblog
      --tmpdir DIR       Parallel per-job buffer dir (default: next to joblog,
                         or \$TMPDIR if set). Avoids filling /tmp.
      --no-compress      Disable parallel --compress (on by default)
      --dry-run          Print first 5 commands and exit
  -h, --help             Show this message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--commands) COMMANDS_FILE="${2:-}"; shift 2 ;;
    -j|--jobs)     JOBS="${2:-}"; shift 2 ;;
    --joblog)      JOBLOG="${2:-}"; shift 2 ;;
    --retries)     RETRIES="${2:-0}"; shift 2 ;;
    --timeout)     TIMEOUT="${2:-}"; shift 2 ;;
    --progress)    PROGRESS=true; shift ;;
    --resume-failed) RESUME_FAILED=true; shift ;;
    --tmpdir)      TMPDIR_ARG="${2:-}"; shift 2 ;;
    --no-compress) NO_COMPRESS=true; shift ;;
    --dry-run)     DRY_RUN=true; shift ;;
    -h|--help)     print_usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; print_usage; exit 2 ;;
  esac
done

[[ -z "$COMMANDS_FILE" ]] && { echo "Error: --commands is required." >&2; print_usage; exit 2; }
[[ -f "$COMMANDS_FILE" ]] || { echo "Error: not found: $COMMANDS_FILE" >&2; exit 2; }
command -v parallel >/dev/null 2>&1 || { echo "Error: GNU parallel not on PATH." >&2; exit 2; }

[[ -n "$JOBS" ]]   || JOBS="$(cpu_count)"
[[ -n "$JOBLOG" ]] || JOBLOG="${COMMANDS_FILE%.*}.joblog.tsv"

# Resolve tmpdir: explicit flag > $TMPDIR > a scratch dir next to the joblog.
# This keeps parallel's per-job stdout/stderr buffers off /tmp, which routinely
# fills up on long batches (e.g. samtools faidx over thousands of groups).
if [[ -z "$TMPDIR_ARG" ]]; then
  if [[ -n "${TMPDIR:-}" ]]; then
    TMPDIR_ARG="$TMPDIR"
  else
    TMPDIR_ARG="$(dirname "$JOBLOG")/.parallel-tmp"
  fi
fi
mkdir -p "$TMPDIR_ARG"

grep -q '[^[:space:]]' "$COMMANDS_FILE" || { echo "Commands file is empty."; exit 0; }

PAR_OPTS=(--jobs "$JOBS" --joblog "$JOBLOG" --will-cite --tmpdir "$TMPDIR_ARG")
[[ "$NO_COMPRESS"   != true ]]   && PAR_OPTS+=(--compress)
[[ "$PROGRESS"      == true ]]   && PAR_OPTS+=(--bar)
[[ "$RETRIES"       -gt 0 ]]     && PAR_OPTS+=(--retries "$RETRIES")
[[ -n "$TIMEOUT" ]]              && PAR_OPTS+=(--timeout "$TIMEOUT")
[[ "$RESUME_FAILED" == true ]]   && PAR_OPTS+=(--resume-failed)

echo "parallel ${PAR_OPTS[*]}"
echo "Commands: $COMMANDS_FILE  |  Jobs: $JOBS  |  Joblog: $JOBLOG"
echo "Tmpdir:   $TMPDIR_ARG  |  Compress: $([[ "$NO_COMPRESS" == true ]] && echo off || echo on)"

if $DRY_RUN; then
  echo "--- dry run: first 5 commands ---"
  head -n 5 "$COMMANDS_FILE"
  exit 0
fi

parallel "${PAR_OPTS[@]}" < "$COMMANDS_FILE"
echo "Done. Joblog: $JOBLOG"