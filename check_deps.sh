#!/usr/bin/env bash
# Verify this machine has everything the pipelines need. Run ONCE after cloning,
# and again after changing conda envs.
#
#   ./check_deps.sh          # report
#   ./check_deps.sh -q       # exit status only, for scripting
#
# WHY THIS IS NOT IN config.local.sh. That file is sourced by every script and by
# every one of the 45 SLURM array tasks, so it must stay cheap -- it defines
# paths and stats them, nothing more. Checking R packages means starting R, about
# 1-2 s a time, for an answer that only changes when somebody installs software.
# Paying that per array task would be a slow way to learn nothing new.
#
# So: paths are asserted continuously (config.local.sh), toolchain is verified
# once (here). Both fail loudly; neither guesses.

set -uo pipefail          # NOT -e: we want every failure reported, not just the first

QUIET=0
[[ "${1:-}" == "-q" ]] && QUIET=1

_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -r "$_root/config.local.sh" ]]; then
    echo "FAIL  config.local.sh missing at $_root" >&2
    echo "      cp config.example.sh config.local.sh and edit it. See README.md." >&2
    exit 1
fi
# shellcheck disable=SC1090
source "$_root/config.local.sh" || { echo "FAIL  config.local.sh did not load cleanly" >&2; exit 1; }

MISSING=(); WARNED=(); OK=0
say()  { [[ $QUIET -eq 1 ]] || printf "%s\n" "$*"; }
pass() { OK=$((OK+1)); [[ $QUIET -eq 1 ]] || printf "  \033[32mok\033[0m    %-34s %s\n" "$1" "${2:-}"; }
fail() { MISSING+=("$1${2:+ -- $2}"); [[ $QUIET -eq 1 ]] || printf "  \033[31mMISS\033[0m  %-34s %s\n" "$1" "${2:-}"; }
# warn() is for things that block ONE optional deliverable rather than the
# pipeline. It never touches the exit status -- see the R packages section.
warn() { WARNED+=("$1${2:+ -- $2}"); [[ $QUIET -eq 1 ]] || printf "  \033[33mwarn\033[0m  %-34s %s\n" "$1" "${2:-}"; }

# --- conda envs ---------------------------------------------------------------
say; say "conda environments"
if [[ -r "$CONDA_SH" ]]; then
    pass "conda.sh" "$CONDA_SH"
    # shellcheck disable=SC1090
    set +u; source "$CONDA_SH"; set -u
    for e in "$ENV_VISOR" "$ENV_MANTA" "$ENV_SVTYPER" "$ENV_GATK" "$ENV_FACET" "$ENV_SVCLONE"; do
        [[ -d "$e/conda-meta" ]] && pass "env: $e" || fail "env: $e" "run tools/setup_rockfish_tool_envs.sh"
    done
else
    fail "conda.sh" "$CONDA_SH not readable"
fi

# --- binaries, each in the env that should provide it -------------------------
# name:env:install-hint   (env "-" means "must be on PATH already")
say; say "executables"
BINS=(
  "VISOR:$ENV_VISOR:conda install -c bioconda visor"
  "samtools:$ENV_VISOR:conda install -c bioconda samtools"
  "minimap2:$ENV_VISOR:conda install -c bioconda minimap2"
  "python3:$ENV_VISOR:-"
  "configManta.py:$ENV_MANTA:conda install -c bioconda manta"
  "svtyper:$ENV_SVTYPER:conda install -c bioconda svtyper"
)
for spec in "${BINS[@]}"; do
    IFS=: read -r bin env hint <<<"$spec"
    if [[ -x "$env/bin/$bin" ]] \
       || conda run --prefix "$env" which "$bin" >/dev/null 2>&1; then
        pass "$bin" "($env)"
    else
        fail "$bin" "in env $env: $hint"
    fi
done

# convertInversion.py is exposed in Manta's bin directory, and the chrX
# calling stage cannot represent inversions without it -- 2,190 of the 4,935
# simulated SVs are inversions, so a missing copy is not a cosmetic gap.
_manta_prefix="$ENV_MANTA"
if [[ -x "$_manta_prefix/bin/convertInversion.py" ]]; then
    pass "convertInversion.py" "($ENV_MANTA/bin)"
else
    fail "convertInversion.py" "expected in $ENV_MANTA/bin; inversions will not be called"
fi

# --- python modules -----------------------------------------------------------
say; say "python modules ($ENV_VISOR)"
for m in pysam; do
    if "$ENV_VISOR/bin/python" -c "import $m" >/dev/null 2>&1; then
        pass "python: $m"
    else
        fail "python: $m" "install $m in $ENV_VISOR"
    fi
done

# --- R packages ---------------------------------------------------------------
# Two groups, and the difference between them decides the exit status.
#
# REQUIRED -- SVCFit's Imports, from its DESCRIPTION. R CMD INSTALL refuses to
# build without every one of them, so a gap blocks installing SVCFit
# outright, and with it the scoring. Hard failure.
#
# FIGURE-ONLY -- what the figure scripts additionally need
# to knit Figure 1A/B, including the packages its three sourced prostate helpers
# pull in (data.table, gridExtra, gtools, plyr). Nothing in the simulation or
# scoring path imports any of these, so a gap is a WARNING and does not fail the
# gate. A machine set up only to run the pipeline is a legitimate configuration,
# and if it could never show a green result people would learn to ignore a red
# one -- which costs more than the gap being reported.
#
# Both groups are resolved in ONE Rscript call: starting R costs 1-2 s and the
# answer is the same either way, so the group is a tag on each output line.
say; say "R packages (RLIB=$RLIB)"
if [[ -x "$RSCRIPT" ]]; then
    _rout="$("$RSCRIPT" -e '
      lib <- Sys.getenv("RLIB"); if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
      req <- c("dplyr","stringr","GenomicRanges","IRanges","S4Vectors","tidyr",
               "ggplot2","reticulate","RColorBrewer","igraph")
      fig <- c("devtools","knitr","rmarkdown","purrr","ggpubr","rstatix","ccube",
               "data.table","gridExtra","gtools","plyr")
      for (p in req) cat("req", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING", "\n")
      for (p in fig) cat("fig", p, if (requireNamespace(p, quietly=TRUE)) "OK" else "MISSING", "\n")
    ' 2>/dev/null)"
    if [[ -z "$_rout" ]]; then
        fail "R packages" "Rscript produced no output; is $RSCRIPT healthy?"
    else
        _fig_header=0
        while read -r grp p st; do
            [[ -z "$p" ]] && continue
            # ccube is not on CRAN; the others are.
            if [[ "$p" == "ccube" ]]; then _hint="remotes::install_github('keyuan/ccube')"
            else _hint="install.packages('$p', lib='$RLIB')"; fi

            if [[ "$grp" == "fig" && $_fig_header -eq 0 ]]; then
                _fig_header=1
                say; say "R packages for figure rendering only"
            fi

            if [[ "$st" == "OK" ]]; then pass "R: $p"
            elif [[ "$grp" == "fig" ]]; then warn "R: $p" "$_hint"
            else fail "R: $p" "$_hint"; fi
        done <<<"$_rout"
    fi
else
    fail "Rscript" "cannot check R packages without it"
fi

# SVclone launches the R interpreter from its own legacy environment. Checking
# these packages in the core SVCFit library is insufficient because PATH order
# can otherwise select the wrong Rscript only after an expensive array starts.
say; say "SVclone R packages ($ENV_SVCLONE)"
if [[ -x "$ENV_SVCLONE/bin/Rscript" ]]; then
    for p in ccube foreach doParallel; do
        if "$ENV_SVCLONE/bin/Rscript" -e \
            "quit(status=if(requireNamespace('$p',quietly=TRUE)) 0L else 1L)" \
            >/dev/null 2>&1; then
            pass "SVclone R: $p"
        else
            fail "SVclone R: $p" "rebuild with tools/setup_rockfish_tool_envs.sh"
        fi
    done
else
    fail "SVclone Rscript" "$ENV_SVCLONE/bin/Rscript is not executable"
fi

# --- reticulate -> sklearn (stage 6 clustering only) --------------------------
# Checked through reticulate rather than by importing sklearn directly, because
# the two can disagree and only this one matters: cluster_data() has reticulate
# EMBED the interpreter, which needs a libpython to link against. An interpreter
# that imports sklearn perfectly well from a shell can still fail here.
#
# A warning, not a failure: only stage 6 needs it, and it is the last stage.
say; say "reticulate -> sklearn (stage 6 only)"
if [[ -n "${RETICULATE_PYTHON:-}" && -x "${RETICULATE_PYTHON:-}" ]]; then
    if RETICULATE_PYTHON="$RETICULATE_PYTHON" "$RSCRIPT" -e '
          lib <- Sys.getenv("RLIB"); if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
          q <- tryCatch({ reticulate::import("sklearn"); 0L }, error = function(e) 1L)
          quit(status = q)' >/dev/null 2>&1; then
        pass "sklearn via reticulate" "$(basename "$(dirname "$(dirname "$RETICULATE_PYTHON")")") env"
    else
        warn "sklearn via reticulate" "RETICULATE_PYTHON=$RETICULATE_PYTHON cannot import sklearn; stage 6 will fail"
    fi
else
    warn "sklearn via reticulate" "RETICULATE_PYTHON unset or not executable; see config.example.sh"
fi

# --- reference data -----------------------------------------------------------
# .fai is checked separately: Manta and samtools both need the index, and its
# absence surfaces late and confusingly rather than at startup.
say; say "reference data"
for v in REF_CHRX REF_AUTO; do
    f="${!v:-}"
    [[ -n "$f" && -s "$f" ]] && pass "$v" "$(basename "$f")" || fail "$v" "${f:-unset}"
    [[ -n "$f" && -s "$f.fai" ]] && pass "$v.fai" || fail "$v.fai" "samtools faidx $f"
done

# --- SLURM --------------------------------------------------------------------
say; say "cluster"
command -v sbatch >/dev/null 2>&1 && pass "sbatch" "$(command -v sbatch)" \
    || fail "sbatch" "not on PATH; check SLURM_BIN in config.local.sh"

# --- verdict ------------------------------------------------------------------
# Warnings are printed whether or not anything failed, but only MISSING sets the
# exit status. -q stays silent and returns that status, so the setup step in
# README.md still means "the pipelines can run here".
say
if ((${#WARNED[@]} > 0)); then
    say "${#WARNED[@]} optional dependenc$( ((${#WARNED[@]}==1)) && echo y || echo ies) absent:"
    for w in "${WARNED[@]}"; do say "  - $w"; done
    say
    say "These are needed ONLY for figure rendering. The simulation, scoring"
    say "and calling stages do not use them, and the gate below ignores them."
    say
fi

if ((${#MISSING[@]} == 0)); then
    say "ALL DEPENDENCIES PRESENT ($OK checks passed)"
    exit 0
fi
say "MISSING ${#MISSING[@]} of $((OK + ${#MISSING[@]})) required dependencies:"
for m in "${MISSING[@]}"; do say "  - $m"; done
say
say "Fix these before running the pipelines. R packages install into RLIB"
say "($RLIB) rather than the conda env, so the env is never mutated."
exit 1
