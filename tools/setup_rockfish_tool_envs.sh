#!/usr/bin/env bash
# Install the project-owned tool environments at absolute prefixes.

set -Eeuo pipefail
umask 0027

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"
project_root="$(cd "$repo_root/../.." && pwd -P)"
software_root="$project_root/01_software"
runtime_loader="$software_root/runtime/rockfish-runtime.sh"
environment_specs="$repo_root/prostate_mixture/prostate_replicates/conda_envs"
curated_specs="$script_dir/environments"
visor_source_url="https://github.com/davidebolo1993/VISOR.git"
visor_source_commit="11afd676a2bbdcc8eddd586e2f57480c39189044"
visor_source="$software_root/vendor/VISOR-$visor_source_commit"

[[ -r "$runtime_loader" ]] || {
  printf 'ERROR: core runtime loader is absent: %s\n' "$runtime_loader" >&2
  exit 1
}
# shellcheck disable=SC1090
source "$runtime_loader"

# The caller may already have an unrelated base environment activated. This
# installer always invokes the project Conda by absolute path, so discard only
# activation bookkeeping inherited from the caller before resolving packages.
unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL
unset CONDA_EXE CONDA_PYTHON_EXE

conda_bin="$ROCKFISH_CONDA_ROOT/bin/conda"
environment_parent="$software_root/conda-envs"
package_cache="$software_root/conda-pkgs"
lock_dir="$ROCKFISH_RUNTIME_ROOT/.tool-env-install.lock"
solver="${ROCKFISH_CONDA_SOLVER:-libmamba}"

[[ -x "$conda_bin" ]] || { printf 'ERROR: missing Conda: %s\n' "$conda_bin" >&2; exit 1; }
[[ "$project_root" == "$ROCKFISH_ROOT" ]] || {
  printf 'ERROR: derived project root %s differs from configured root %s\n' "$project_root" "$ROCKFISH_ROOT" >&2
  exit 1
}

mkdir -p "$environment_parent" "$package_cache"
mkdir "$lock_dir" 2>/dev/null || {
  printf 'ERROR: tool-environment installation lock exists: %s\n' "$lock_dir" >&2
  exit 1
}
trap 'rmdir "$lock_dir" 2>/dev/null || true' EXIT

export CONDARC="$ROCKFISH_RUNTIME_ROOT/condarc"
export CONDA_PKGS_DIRS="$package_cache"
export CONDA_ENVS_PATH="$environment_parent"
export CONDA_REGISTER_ENVS=false
export CONDA_CHANNEL_PRIORITY="${ROCKFISH_CONDA_CHANNEL_PRIORITY:-flexible}"
export MAMBA_ROOT_PREFIX="$ROCKFISH_CONDA_ROOT"
export XDG_CACHE_HOME="$ROCKFISH_RUNTIME_ROOT/cache"
export PIP_CACHE_DIR="$ROCKFISH_RUNTIME_ROOT/cache/pip"

requested=("$@")
((${#requested[@]})) || requested=(visor manta svtyper gatk4 facet svclone dnacopy)

valid_name() {
  case "$1" in
    visor|manta|svtyper|gatk4|facet|svclone|dnacopy) return 0 ;;
    *) return 1 ;;
  esac
}

environment_ready() {
  local name="$1" prefix="$2"
  case "$name" in
    visor) [[ -x "$prefix/bin/VISOR" && -x "$prefix/bin/bedtools" && -x "$prefix/bin/minimap2" && -x "$prefix/bin/samtools" ]] ;;
    manta) [[ -x "$prefix/bin/configManta.py" && -x "$prefix/bin/convertInversion.py" ]] ;;
    svtyper) [[ -x "$prefix/bin/svtyper" ]] ;;
    gatk4) [[ -x "$prefix/bin/gatk" ]] ;;
    facet) [[ -x "$prefix/bin/Rscript" ]] ;;
    dnacopy)
      [[ -x "$prefix/bin/Rscript" ]] &&
        "$prefix/bin/Rscript" -e 'quit(status=ifelse(requireNamespace("DNAcopy", quietly=TRUE), 0, 1))' >/dev/null 2>&1
      ;;
    svclone)
      [[ -x "$prefix/bin/svclone" && -x "$prefix/bin/Rscript" ]] &&
        "$prefix/bin/python" -c 'import importlib.metadata as m, numpy, pandas, sys; assert sys.version_info[:2] == (3, 10); assert int(numpy.__version__.split(".")[0]) < 2; assert int(pandas.__version__.split(".")[0]) < 3; assert int(m.version("setuptools").split(".")[0]) < 81' >/dev/null 2>&1 &&
        "$prefix/bin/Rscript" -e 'quit(status=ifelse(all(vapply(c("ccube", "dplyr", "optparse"), requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))), 0, 1))' >/dev/null 2>&1
      ;;
  esac
}

install_visor_source() {
  local prefix="$1" actual_commit
  if [[ ! -d "$visor_source/.git" ]]; then
    [[ ! -e "$visor_source" ]] || {
      printf 'ERROR: non-repository path blocks pinned VISOR source: %s\n' "$visor_source" >&2
      return 1
    }
    mkdir -p "$software_root/vendor"
    git clone "$visor_source_url" "$visor_source"
  fi
  git -C "$visor_source" fetch --tags origin
  git -C "$visor_source" checkout --detach "$visor_source_commit"
  actual_commit="$(git -C "$visor_source" rev-parse HEAD)"
  [[ "$actual_commit" == "$visor_source_commit" ]] || {
    printf 'ERROR: VISOR source is at %s, expected %s\n' "$actual_commit" "$visor_source_commit" >&2
    return 1
  }
  "$prefix/bin/python" -m pip install --no-deps "$visor_source"
  printf '%s  %s\n' "$visor_source_commit" "$visor_source_url" > "$prefix/visor-source.txt"
}

for name in "${requested[@]}"; do
  valid_name "$name" || { printf 'ERROR: unknown tool environment: %s\n' "$name" >&2; exit 2; }
  prefix="$environment_parent/$name"
  spec="$environment_specs/${name}_portable.yml"
  # A full GATK environment export bundles unrelated HTS utilities and Perl.
  # Keep GATK isolated; shared HTS commands come from the core runtime.
  [[ "$name" == gatk4 || "$name" == visor || "$name" == svclone || "$name" == dnacopy ]] && spec="$curated_specs/${name}.yml"
  [[ -r "$spec" ]] || { printf 'ERROR: missing environment definition: %s\n' "$spec" >&2; exit 1; }

  if ! environment_ready "$name" "$prefix"; then
    if [[ -d "$prefix/conda-meta" ]]; then
      printf 'Updating incomplete %s environment at %s\n' "$name" "$prefix"
      "$conda_bin" env update --solver="$solver" --prefix "$prefix" --file "$spec"
    else
      [[ ! -e "$prefix" ]] || { printf 'ERROR: incomplete non-Conda prefix exists: %s\n' "$prefix" >&2; exit 1; }
      printf 'Creating %s at %s\n' "$name" "$prefix"
      "$conda_bin" env create --yes --solver="$solver" --prefix "$prefix" --file "$spec"
    fi
    [[ "$name" != visor ]] || install_visor_source "$prefix"
  else
    printf 'Using existing environment: %s\n' "$prefix"
  fi

  environment_ready "$name" "$prefix" || {
    printf 'ERROR: required executables are missing after installing %s\n' "$name" >&2
    exit 1
  }

  "$conda_bin" list --prefix "$prefix" --explicit > "$prefix/conda-explicit-linux-64.txt"
  cp -p "$spec" "$prefix/environment.yml"
  sha256sum "$prefix/environment.yml" > "$prefix/environment-source.sha256"
  {
    printf 'recorded_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'host=%s\n' "$(hostname -f)"
    printf 'prefix=%s\n' "$prefix"
    printf 'source=%s\n' "$spec"
    printf 'workflow_commit=%s\n' "$(git -C "$repo_root" rev-parse HEAD)"
  } > "$prefix/creation-record.txt"
done

for name in "${requested[@]}"; do
  chmod -R g+rX,o-rwx "$environment_parent/$name"
done
printf 'Project tool environments are installed under %s\n' "$environment_parent"
