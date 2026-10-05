#!/usr/bin/env bash
# Install the project-owned ROCKFISH core SVCFit runtime.

set -euo pipefail
umask 0027

readonly script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly workflow_repo="$(git -C "$script_dir" rev-parse --show-toplevel)"
readonly project_root="$(cd "$workflow_repo/../.." && pwd -P)"
readonly software_root="$project_root/01_software"
readonly runtime_root="$software_root/runtime"
readonly environment_root="$software_root/conda-envs/rockfish-svcfit-20260920"
readonly package_cache="$software_root/conda-pkgs"
readonly miniforge_version=26.7.2-0
readonly installer_name="Miniforge3-${miniforge_version}-Linux-x86_64.sh"
readonly installer_sha256=281b0ac7d550802efc81af633225a5e6116d29ae72f3ab4eae7168c3931a4c05
readonly installer_url="https://github.com/conda-forge/miniforge/releases/download/${miniforge_version}/${installer_name}"
readonly installer_dir="$software_root/installers"
readonly installer="$installer_dir/$installer_name"
readonly miniforge_root="$software_root/miniforge3-${miniforge_version}"
readonly loader="$runtime_root/rockfish-runtime.sh"
readonly environment_yaml="$script_dir/rockfish-runtime.yml"
readonly svcfit_source="$software_root/SVCFit"
readonly lock_dir="$runtime_root/.install.lock"

[[ "$(pwd -P)" == "$project_root" ]] || {
  printf 'ERROR: run from %s\n' "$project_root" >&2
  exit 1
}
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || {
  printf 'ERROR: pinned runtime supports Linux x86_64 only\n' >&2
  exit 1
}
[[ -r "$environment_yaml" ]] || {
  printf 'ERROR: missing environment YAML: %s\n' "$environment_yaml" >&2
  exit 1
}
[[ "$(git -C "$svcfit_source" rev-parse HEAD)" == e0d7e0b8d704caa3cfb227bc3dbf1ee99d66fac7 ]] || {
  printf 'ERROR: SVCFit source is not at the required commit\n' >&2
  exit 1
}
[[ -z "$(git -C "$svcfit_source" status --porcelain)" ]] || {
  printf 'ERROR: SVCFit source is dirty\n' >&2
  exit 1
}

mkdir -p "$runtime_root" "$installer_dir" "$software_root/conda-envs" "$package_cache"
mkdir "$lock_dir" 2>/dev/null || {
  printf 'ERROR: runtime installation lock exists: %s\n' "$lock_dir" >&2
  exit 1
}
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT

# Keep all mutable package-manager state in the project. In particular, do not
# register environments or cache data below a user's home directory.
export CONDARC="$runtime_root/condarc"
export CONDA_PKGS_DIRS="$package_cache"
export CONDA_ENVS_PATH="$software_root/conda-envs"
export CONDA_REGISTER_ENVS=false
export MAMBA_ROOT_PREFIX="$miniforge_root"
export XDG_CACHE_HOME="$runtime_root/cache"
export PIP_CACHE_DIR="$runtime_root/cache/pip"
export R_USER="$runtime_root/r-user"
mkdir -p "$XDG_CACHE_HOME" "$PIP_CACHE_DIR" "$R_USER"
printf 'channels:\n  - conda-forge\n  - bioconda\nchannel_priority: strict\nregister_envs: false\n' > "$CONDARC"

if [[ ! -x "$miniforge_root/bin/conda" ]]; then
  if [[ ! -f "$installer" ]]; then
    curl --fail --location --output "$installer.part" "$installer_url"
    mv "$installer.part" "$installer"
  fi
  printf '%s  %s\n' "$installer_sha256" "$installer" | sha256sum --check --status || {
    printf 'ERROR: installer checksum mismatch: %s\n' "$installer" >&2
    exit 1
  }
  bash "$installer" -b -p "$miniforge_root"
fi

if [[ ! -x "$environment_root/bin/Rscript" ]]; then
  # Use Conda's classic solver/transaction implementation on ROCKFISH NFS. The
  # libmamba transaction path produced an incomplete Perl cache extraction in
  # the first preserved attempt from this session.
  "$miniforge_root/bin/conda" env create --yes --solver=classic \
    --prefix "$environment_root" --file "$environment_yaml"
fi

export PATH="$environment_root/bin:$PATH"

# The tree-evaluation notebook and publication figure need these packages.
# Install them explicitly when resuming an environment created from an earlier
# revision of the YAML; ccube is supplied by Bioconda rather than CRAN.
if ! "$environment_root/bin/Rscript" -e 'ok <- identical(as.character(getRversion()), "4.4.3") && all(vapply(c("lme4", "patchwork", "devtools", "rmarkdown", "ggpubr", "rstatix", "ccube", "data.table", "gridExtra", "gtools", "plyr"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))); quit(status=ifelse(ok, 0, 1))'; then
  "$miniforge_root/bin/mamba" install --yes \
    --prefix "$environment_root" -c conda-forge -c bioconda \
    'r-base=4.4.3' 'r-lme4' 'r-patchwork' 'r-devtools' 'r-rmarkdown' \
    'r-ggpubr' 'r-rstatix' 'r-ccube' 'r-data.table' 'r-gridextra' \
    'r-gtools' 'r-plyr'
fi

# Source-installed Bioconductor packages compile against zlib. Conda's
# libzlib runtime does not include zlib.h, so ensure the matching development
# payload is present even when resuming a partially completed installation.
if [[ ! -r "$environment_root/include/zlib.h" ]]; then
  "$miniforge_root/bin/mamba" install --yes \
    --prefix "$environment_root" -c conda-forge 'zlib=1.3.2'
fi

# Bioconda's GenomeInfoDbData post-link script creates a quoted temporary
# filename rejected by ROCKFISH NFS. Install the matching Bioconductor 3.20
# packages through R into this prefix instead.
if ! "$environment_root/bin/Rscript" -e 'quit(status=ifelse(requireNamespace("GenomicRanges", quietly=TRUE), 0, 1))'; then
  "$environment_root/bin/Rscript" -e '
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors"),
                     version="3.20", ask=FALSE, update=FALSE)
'
fi

# bcftools metadata pulls Perl for optional helper scripts. ROCKFISH NFS rejects
# Perl manpage filenames containing "::". HTSlib is already supplied by
# Samtools, so install the pinned bcftools payload without expanding that
# optional dependency chain.
if [[ ! -x "$environment_root/bin/bcftools" ]]; then
  "$miniforge_root/bin/mamba" install --yes --no-deps \
    --prefix "$environment_root" -c bioconda 'bcftools=1.23.1'
fi

required_executables=(R Rscript python bcftools samtools bedtools bgzip tabix)
for executable in "${required_executables[@]}"; do
  [[ -x "$environment_root/bin/$executable" ]] || {
    printf 'ERROR: missing executable: %s/bin/%s\n' "$environment_root" "$executable" >&2
    exit 1
  }
done

if ! "$environment_root/bin/Rscript" -e 'quit(status=ifelse(requireNamespace("SVCFit", quietly=TRUE), 0, 1))'; then
  "$environment_root/bin/R" CMD INSTALL --no-multiarch "$svcfit_source"
fi

RETICULATE_PYTHON="$environment_root/bin/python" "$environment_root/bin/Rscript" -e '
required <- c("SVCFit", "optparse", "dplyr", "tidyr", "stringr", "purrr",
              "GenomicRanges", "IRanges", "S4Vectors", "readr", "ggplot2",
              "RColorBrewer", "reticulate", "igraph", "testthat", "lme4",
              "patchwork", "devtools", "rmarkdown", "ggpubr", "rstatix",
              "ccube", "data.table", "gridExtra", "gtools", "plyr")
missing <- required[!vapply(required, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(missing)) stop("missing R packages: ", paste(missing, collapse=", "))
reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required=TRUE)
reticulate::import("sklearn")
cat("ROCKFISH R packages and reticulate -> sklearn: OK\n")
'

"$miniforge_root/bin/conda" list --prefix "$environment_root" --explicit \
  > "$environment_root/conda-explicit-linux-64.txt"
ROCKFISH_R_LOCK="$environment_root/r-packages.tsv" "$environment_root/bin/Rscript" -e '
x <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
write.table(x, file=Sys.getenv("ROCKFISH_R_LOCK"), sep="\t", quote=FALSE, row.names=FALSE)
'
cp -p "$environment_yaml" "$environment_root/environment.yml"
sha256sum "$environment_root/environment.yml" > "$environment_root/environment-source.sha256"
{
  printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname -f)"
  printf 'installer=%s\n' "$installer"
  printf 'installer_sha256=%s\n' "$installer_sha256"
  printf 'svcfit_commit=%s\n' "$(git -C "$svcfit_source" rev-parse HEAD)"
} > "$environment_root/creation-record.txt"

loader_tmp="$loader.tmp.$$"
cat > "$loader_tmp" <<EOF
# Generated project-owned ROCKFISH runtime. Source this file; do not activate by name.
export ROCKFISH_ROOT="$project_root"
export ROCKFISH_RUNTIME_ROOT="$runtime_root"
export ROCKFISH_CONDA_ROOT="$miniforge_root"
export ROCKFISH_SVCFIT_ENV="$environment_root"
export CONDA_SH="$miniforge_root/etc/profile.d/conda.sh"
export ENV_SVCFIT="$environment_root"
export ENV_VISOR="$software_root/conda-envs/visor"
export ENV_MANTA="$software_root/conda-envs/manta"
export ENV_SVTYPER="$software_root/conda-envs/svtyper"
export ENV_GATK="$software_root/conda-envs/gatk4"
export ENV_FACET="$software_root/conda-envs/facet"
export ENV_SVCLONE="$software_root/conda-envs/svclone"
export RETICULATE_PYTHON="$environment_root/bin/python"
export RSCRIPT="$environment_root/bin/Rscript"
export SAMTOOLS="$environment_root/bin/samtools"
export BCFTOOLS="$environment_root/bin/bcftools"
export BEDTOOLS="$environment_root/bin/bedtools"
export SVCFIT_PKG_DIR="$svcfit_source"
export REF_AUTO="$project_root/02_data/tree_eval/input_data/reference/chr1-2.fa"
export REF_CHRX="$project_root/02_data/references/genomes/chr22-X.fa"
export PROSTATE_REF="$project_root/02_data/prostate_mixture/ref/hs37d5.fa"
export SLURM_BIN="${SLURM_BIN:-}"
export PATH="$environment_root/bin:\$PATH"
EOF
chmod 0640 "$loader_tmp"
mv "$loader_tmp" "$loader"

chmod 2770 "$runtime_root" "$installer_dir" "$software_root/conda-envs" "$package_cache" "$software_root/r-libs"
chmod -R g+rX,o-rwx "$miniforge_root" "$runtime_root" "$software_root/conda-envs" "$package_cache" "$software_root/r-libs"
printf 'ROCKFISH core runtime ready: %s\n' "$loader"
