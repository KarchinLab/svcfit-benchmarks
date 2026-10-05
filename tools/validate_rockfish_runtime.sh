#!/usr/bin/env bash
# Validate the core ROCKFISH runtime using only absolute project-owned tools.

set -euo pipefail

readonly script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly workflow_repo="$(git -C "$script_dir" rev-parse --show-toplevel)"
readonly project_root="$(cd "$workflow_repo/../.." && pwd -P)"
readonly loader="$project_root/01_software/runtime/rockfish-runtime.sh"

[[ -r "$loader" ]] || { printf 'ERROR: missing loader: %s\n' "$loader" >&2; exit 1; }
if grep -En '/home/|\$HOME|~/' "$loader"; then
  printf 'ERROR: personal path found in runtime loader\n' >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$loader"

required_executables=(
  "$RSCRIPT" "$RETICULATE_PYTHON" "$SAMTOOLS" "$BCFTOOLS" "$BEDTOOLS"
  "$ROCKFISH_SVCFIT_ENV/bin/bgzip" "$ROCKFISH_SVCFIT_ENV/bin/tabix"
  "$ENV_VISOR/bin/VISOR" "$ENV_VISOR/bin/python" "$ENV_VISOR/bin/samtools" "$ENV_VISOR/bin/minimap2"
  "$ENV_MANTA/bin/configManta.py" "$ENV_MANTA/bin/convertInversion.py"
  "$ENV_SVTYPER/bin/svtyper" "$ENV_GATK/bin/gatk"
  "$ENV_FACET/bin/Rscript" "$ENV_SVCLONE/bin/svclone"
)
for executable in "${required_executables[@]}"; do
  [[ -x "$executable" ]] || { printf 'ERROR: missing executable: %s\n' "$executable" >&2; exit 1; }
  case "$executable" in
    /home/*) printf 'ERROR: personal executable: %s\n' "$executable" >&2; exit 1 ;;
  esac
done

tool_environments=("$ENV_VISOR" "$ENV_MANTA" "$ENV_SVTYPER" "$ENV_GATK" "$ENV_FACET" "$ENV_SVCLONE")
for environment in "${tool_environments[@]}"; do
  [[ "$environment" == "$ROCKFISH_ROOT"/01_software/conda-envs/* ]] || {
    printf 'ERROR: tool environment escapes project software root: %s\n' "$environment" >&2
    exit 1
  }
  [[ -s "$environment/conda-explicit-linux-64.txt" && -s "$environment/environment.yml" && -s "$environment/creation-record.txt" ]] || {
    printf 'ERROR: environment provenance is incomplete: %s\n' "$environment" >&2
    exit 1
  }
done
"$ENV_VISOR/bin/python" -c 'import pysam, VISOR'
"$ENV_FACET/bin/Rscript" -e 'stopifnot(requireNamespace("facets", quietly=TRUE))'
"$ENV_SVCLONE/bin/python" -c 'import importlib.metadata as m, numpy, pandas, sys; assert sys.version_info[:2] == (3, 10); assert int(numpy.__version__.split(".")[0]) < 2; assert int(pandas.__version__.split(".")[0]) < 3; assert int(m.version("setuptools").split(".")[0]) < 81'
"$ENV_SVCLONE/bin/Rscript" -e 'stopifnot(all(vapply(c("ccube", "dplyr", "optparse"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))))'

required_references=("$REF_AUTO" "$REF_CHRX" "$PROSTATE_REF")
for reference in "${required_references[@]}"; do
  [[ -s "$reference" ]] || { printf 'ERROR: missing reference: %s\n' "$reference" >&2; exit 1; }
done
[[ -s "$ROCKFISH_SVCFIT_ENV/conda-explicit-linux-64.txt" ]] || {
  printf 'ERROR: explicit runtime lock is missing\n' >&2
  exit 1
}

"$RETICULATE_PYTHON" -c 'import sklearn; print("sklearn", sklearn.__version__)'
RETICULATE_PYTHON="$RETICULATE_PYTHON" "$RSCRIPT" -e '
if (!identical(as.character(getRversion()), "4.4.3"))
  stop("unexpected R version: ", getRversion())
required <- c("SVCFit", "optparse", "dplyr", "tidyr", "stringr", "purrr",
              "GenomicRanges", "IRanges", "S4Vectors", "readr", "ggplot2",
              "RColorBrewer", "reticulate", "igraph", "testthat", "lme4",
              "patchwork")
missing <- required[!vapply(required, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(missing)) stop("missing R packages: ", paste(missing, collapse=", "))
reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required=TRUE)
reticulate::import("sklearn")
cat("ROCKFISH runtime validation: OK\n")
'
printf 'Validated on host %s at %s\n' "$(hostname -f)" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
