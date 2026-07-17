#!/usr/bin/env bash
# Run a commands file in parallel using GNU parallel.
#
# Usage:
#   ./method/parallel_run.sh -c <commands.txt>
#   ./method/parallel_run.sh -c <commands.txt> -j 16 --progress
#   ./method/parallel_run.sh -c <commands.txt> --retries 2 --resume-failed

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

# Reusable integer validator.
require_uint() {
  # $1 = value, $2 = flag name
  [[ "$1" =~ ^[0-9]+$ ]] || { echo "Error: $2 must be a non-negative integer (got: '$1')" >&2; exit 2; }
}

print_usage() {
  cat <<'EOF'
Usage: scripts/parallel_run.sh -c COMMANDS_FILE [options]

Required:
  -c, --commands FILE    Path to commands file (one command per line)

Options:
  -j, --jobs N            Parallel workers (default: CPU count)
      --joblog FILE       GNU parallel joblog path (default: <commands>.joblog.tsv)
      --retries N         Total attempts per job (GNU parallel semantics: N=total
                            runs, not extra retries; N<=1 means no retry) (default: 0)
      --timeout SECS      Per-job timeout in seconds
      --progress          Show live progress bar
      --resume-failed     Re-run only previously failed jobs from an existing joblog
      --tmpdir DIR        Parallel per-job buffer dir (default: \$TMPDIR,
                            then /var/tmp, then /tmp, then next to joblog)
      --no-compress       Disable parallel --compress (on by default)
      --dry-run           Print first 5 commands and exit
  -h, --help              Show this message
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

# Validate numeric inputs before they reach [[ -gt ]] etc.
require_uint "$JOBS" "--jobs"
require_uint "$RETRIES" "--retries"
[[ -n "$TIMEOUT" ]] && require_uint "$TIMEOUT" "--timeout"

[[ -n "$JOBLOG" ]] || JOBLOG="${COMMANDS_FILE}.joblog.tsv"

if [[ "$RESUME_FAILED" == true && ! -f "$JOBLOG" ]]; then
  echo "Error: --resume-failed requires an existing joblog: $JOBLOG" >&2
  exit 2
fi

if [[ "$NO_COMPRESS" != true ]]; then
  if ! command -v pigz >/dev/null 2>&1 \
     && ! command -v gzip >/dev/null 2>&1 \
     && ! command -v zstd >/dev/null 2>&1 \
     && ! command -v lz4  >/dev/null 2>&1; then
    echo "Warning: --compress is on but no compressor (pigz/gzip/zstd/lz4) found;" >&2
    echo "         disabling compression. Use --no-compress to silence this." >&2
    NO_COMPRESS=true
  fi
fi

resolve_tmpdir() {
  local candidate

  if [[ -n "$TMPDIR_ARG" ]]; then
    if ! mkdir -p "$TMPDIR_ARG" 2>/dev/null; then
      echo "Error: cannot create --tmpdir: $TMPDIR_ARG" >&2
      exit 2
    fi
    [[ -w "$TMPDIR_ARG" ]] || { echo "Error: --tmpdir is not writable: $TMPDIR_ARG" >&2; exit 2; }
    return 0
  fi

  for candidate in "${TMPDIR:-}" /var/tmp /tmp "$(dirname "$JOBLOG")/.parallel-tmp"; do
    [[ -n "$candidate" ]] || continue
    mkdir -p "$candidate" 2>/dev/null || continue
    if [[ -w "$candidate" ]]; then
      TMPDIR_ARG="$candidate"
      if [[ "$candidate" == "$(dirname "$JOBLOG")/.parallel-tmp" ]]; then
        echo "Warning: using tmpdir on the project filesystem ($candidate);" >&2
        echo "         set \$TMPDIR or --tmpdir to keep buffers off it." >&2
      fi
      return 0
    fi
  done

  echo "Error: no writable tmpdir found (tried \$TMPDIR, /var/tmp, /tmp, and $(dirname "$JOBLOG")/.parallel-tmp)." >&2
  exit 2
}

# Keep GNU parallel's per-job stdout/stderr buffers off the project filesystem
# whenever possible; full temp volumes are the usual cause of this failure.
resolve_tmpdir

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

set +e
parallel "${PAR_OPTS[@]}" < "$COMMANDS_FILE"
status=$?
set -e

if [[ "$status" -eq 0 ]]; then
  echo "Done. Joblog: $JOBLOG (all jobs succeeded)"
elif [[ "$status" -ge 1 && "$status" -le 100 ]]; then
  echo "Done. Joblog: $JOBLOG ($status job(s) failed)" >&2
elif [[ "$status" -eq 101 ]]; then
  echo "Done. Joblog: $JOBLOG (more than 100 jobs failed)" >&2
else
  echo "Done. Joblog: $JOBLOG (parallel exited with status $status)" >&2
fi

exit "$status"
