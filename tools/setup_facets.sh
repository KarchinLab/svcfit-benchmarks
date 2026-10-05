#!/usr/bin/env bash

# Create the included FACETS conda environment and compile the official
# mskcc/facets snp-pileup source shipped inside its installed R package.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
software_root="$(cd "$repo_root/.." && pwd)"
env_name="${ENV_FACET:-facet}"
env_file="${repo_root}/prostate_mixture/prostate_replicates/conda_envs/facet_portable.yml"
dest="${FACET_HOME:-${software_root}/facets}"
rebuild=0

usage() {
  cat <<EOF
Usage: tools/setup_facets.sh [--rebuild] [--dest DIRECTORY]

Creates conda environment '${env_name}', compiles the official FACETS
snp-pileup program, and prints config.local.sh settings.

Default destination: ${dest}
EOF
}

while (($#)); do
  case "$1" in
    --rebuild) rebuild=1; shift ;;
    --dest) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; dest="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

conda_exe="${CONDA:-$(command -v conda || true)}"
[[ -n "$conda_exe" && -x "$conda_exe" ]] || {
  echo "ERROR: conda is not available; set CONDA to its executable" >&2
  exit 1
}

if ! "$conda_exe" env list | awk '{print $1}' | grep -qx "$env_name"; then
  echo "Creating conda environment '$env_name' from $env_file"
  "$conda_exe" env create -n "$env_name" -f "$env_file"
else
  echo "Conda environment '$env_name' already exists"
fi

env_prefix="$("$conda_exe" env list | awk -v name="$env_name" '$1 == name {print $NF; exit}')"
[[ -n "$env_prefix" && -d "$env_prefix" ]] || {
  echo "ERROR: cannot resolve prefix for conda environment '$env_name'" >&2
  exit 1
}

extcode="$("$conda_exe" run -n "$env_name" Rscript -e \
  'cat(system.file("extcode", package="facets"))')"
[[ -f "$extcode/snp-pileup.cpp" && -f "$extcode/snp-pileup.h" ]] || {
  echo "ERROR: official snp-pileup source is absent from $extcode" >&2
  exit 1
}

mkdir -p "$dest"
binary="$dest/snp-pileup"
driver="$dest/facet.R"
if [[ -x "$binary" && $rebuild -eq 0 ]]; then
  echo "Existing executable retained: $binary"
else
  echo "Compiling official mskcc/facets snp-pileup source"
  "$conda_exe" run -n "$env_name" g++ -std=c++11 \
    -I"$extcode" -I"$env_prefix/include" \
    "$extcode/snp-pileup.cpp" \
    -L"$env_prefix/lib" -lhts -Wl,-rpath,"$env_prefix/lib" \
    -o "$binary"
  chmod 0755 "$binary"
fi

# Install a stable runnable copy while retaining tools/facet.R as the tracked
# source of truth. Re-running this helper updates the installed driver.
install -m 0755 "$repo_root/tools/facet.R" "$driver"

"$conda_exe" run -n "$env_name" Rscript -e \
  'stopifnot(requireNamespace("facets", quietly=TRUE)); cat("facets ", as.character(packageVersion("facets")), "\n", sep="")'
[[ -x "$binary" ]] || { echo "ERROR: build did not create $binary" >&2; exit 1; }
"$conda_exe" run -n "$env_name" Rscript "$driver" --help >/dev/null

cat <<EOF

FACETS setup complete. Put these lines in config.local.sh:

export FACET_HOME="$dest"
export FACET_SNP_PILEUP="\${FACET_HOME}/snp-pileup"
export FACET_R_SCRIPT="$driver"

Then run:
  ./check_deps.sh
EOF
