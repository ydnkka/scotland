#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run all Chapter 4 observation/network build scripts.

Usage:
  bash observation_networks/run_all.sh [options]

Options:
  --workers N                  Worker processes for build_mixing. Default: build_mixing default.
  --log-level LEVEL            Python log level. Default: INFO.
  --conda-env NAME             Run via `conda run -n NAME python`.
  --python PATH                Python executable. Ignored if --conda-env is used.
  --max-windows N              Development cap for build_tables main windows.
  --max-transition-windows N   Development cap for build_tables transition windows.
  --max-mixing-windows N       Development cap for build_mixing windows.
  --mixing-permutations N      Permutations for compatibility assortativity p-values.
  --mixing-permutation-seed N  Base seed for compatibility permutation p-values.
  --mixing-missing-label LABEL Missing node-attribute label passed to build_mixing.
  --mixing-progress-every N    Log build_mixing progress every N pairwise files.
  --skip-simd                  Skip SIMD population-weighting validation.
  --skip-tables                Skip core observation/transition tables.
  --skip-mixing                Skip compatibility-network mixing.
  --skip-figures               Skip figure generation.
  --skip-transition            Pass through to build_tables.
  --dry-run                    Print commands without executing them.
  -h, --help                   Show this help.

Examples:
  bash observation_networks/run_all.sh --workers 5
  bash observation_networks/run_all.sh --conda-env PhD --workers 5
  bash observation_networks/run_all.sh --max-windows 2 --max-transition-windows 3 --max-mixing-windows 1
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
cd "${repo_root}"

workers=""
log_level="INFO"
conda_env=""
python_bin="${PYTHON:-python}"
max_windows=""
max_transition_windows=""
max_mixing_windows=""
mixing_permutations=""
mixing_permutation_seed=""
mixing_missing_label=""
mixing_progress_every=""
skip_simd=0
skip_tables=0
skip_mixing=0
skip_figures=0
skip_transition=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workers)
      workers="${2:?--workers requires a value}"
      shift 2
      ;;
    --log-level)
      log_level="${2:?--log-level requires a value}"
      shift 2
      ;;
    --conda-env)
      conda_env="${2:?--conda-env requires a value}"
      shift 2
      ;;
    --python)
      python_bin="${2:?--python requires a value}"
      shift 2
      ;;
    --max-windows)
      max_windows="${2:?--max-windows requires a value}"
      shift 2
      ;;
    --max-transition-windows)
      max_transition_windows="${2:?--max-transition-windows requires a value}"
      shift 2
      ;;
    --max-mixing-windows)
      max_mixing_windows="${2:?--max-mixing-windows requires a value}"
      shift 2
      ;;
    --mixing-permutations)
      mixing_permutations="${2:?--mixing-permutations requires a value}"
      shift 2
      ;;
    --mixing-permutation-seed)
      mixing_permutation_seed="${2:?--mixing-permutation-seed requires a value}"
      shift 2
      ;;
    --mixing-missing-label)
      mixing_missing_label="${2:?--mixing-missing-label requires a value}"
      shift 2
      ;;
    --mixing-progress-every)
      mixing_progress_every="${2:?--mixing-progress-every requires a value}"
      shift 2
      ;;
    --skip-simd)
      skip_simd=1
      shift
      ;;
    --skip-tables)
      skip_tables=1
      shift
      ;;
    --skip-mixing)
      skip_mixing=1
      shift
      ;;
    --skip-figures)
      skip_figures=1
      shift
      ;;
    --skip-transition)
      skip_transition=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -n "${conda_env}" ]]; then
  python_cmd=(conda run -n "${conda_env}" python)
else
  python_cmd=("${python_bin}")
fi

if [[ "${dry_run}" -eq 0 ]]; then
  export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
  export MPLCONFIGDIR="${MPLCONFIGDIR:-${repo_root}/observation_networks/results/intermediate/matplotlib}"
  mkdir -p "${MPLCONFIGDIR}"
fi

run_cmd() {
  printf '\n==> '
  printf '%q ' "$@"
  printf '\n'
  if [[ "${dry_run}" -eq 0 ]]; then
    "$@"
  fi
}

if [[ "${skip_simd}" -eq 0 ]]; then
  run_cmd "${python_cmd[@]}" -m observation_networks.build_simd_validation \
    --log-level "${log_level}"
fi

if [[ "${skip_tables}" -eq 0 ]]; then
  table_args=(--log-level "${log_level}")
  if [[ -n "${max_windows}" ]]; then
    table_args+=(--max-windows "${max_windows}")
  fi
  if [[ -n "${max_transition_windows}" ]]; then
    table_args+=(--max-transition-windows "${max_transition_windows}")
  fi
  if [[ "${skip_transition}" -eq 1 ]]; then
    table_args+=(--skip-transition)
  fi
  run_cmd "${python_cmd[@]}" -m observation_networks.build_tables "${table_args[@]}"
fi

if [[ "${skip_mixing}" -eq 0 ]]; then
  mixing_args=(--all-windows --log-level "${log_level}")
  if [[ -n "${workers}" ]]; then
    mixing_args+=(--workers "${workers}")
  fi
  if [[ -n "${max_mixing_windows}" ]]; then
    mixing_args+=(--max-windows "${max_mixing_windows}")
  fi
  if [[ -n "${mixing_permutations}" ]]; then
    mixing_args+=(--n-permutations "${mixing_permutations}")
  fi
  if [[ -n "${mixing_permutation_seed}" ]]; then
    mixing_args+=(--permutation-seed "${mixing_permutation_seed}")
  fi
  if [[ -n "${mixing_missing_label}" ]]; then
    mixing_args+=(--missing-label "${mixing_missing_label}")
  fi
  if [[ -n "${mixing_progress_every}" ]]; then
    mixing_args+=(--progress-every "${mixing_progress_every}")
  fi
  run_cmd "${python_cmd[@]}" -m observation_networks.build_mixing "${mixing_args[@]}"
fi

if [[ "${skip_figures}" -eq 0 ]]; then
  run_cmd "${python_cmd[@]}" -m observation_networks.make_figures \
    --skip-missing \
    --log-level "${log_level}"
fi

printf '\nChapter 4 observation/network build complete.\n'
