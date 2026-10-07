#!/usr/bin/env bash
# Install the pinned SVCFit checkout into an isolated library and run its tests.

set -Eeuo pipefail

: "${RUN_ROOT:?RUN_ROOT is required}"
: "${RUN_LIBRARY:?RUN_LIBRARY is required}"
: "${WORKFLOW_REPO:?WORKFLOW_REPO is required}"
: "${EXPECTED_WORKFLOW_COMMIT:?EXPECTED_WORKFLOW_COMMIT is required}"
: "${EXPECTED_SVCFIT_COMMIT:?EXPECTED_SVCFIT_COMMIT is required}"
: "${SVCFIT_PKG_DIR:?SVCFIT_PKG_DIR is required}"
: "${RSCRIPT:?RSCRIPT is required}"
: "${ENV_SVCFIT:?ENV_SVCFIT is required}"
: "${RLIB:?RLIB is required}"
: "${RETICULATE_PYTHON:?RETICULATE_PYTHON is required}"

write_failure() {
    local status=$1
    rm -f "$RUN_ROOT/PASS"
    printf 'failed_at=%s\nexit_status=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$status" > "$RUN_ROOT/FAIL"
    exit "$status"
}
fail() {
    local status=$?
    write_failure "$status"
}
die() {
    printf 'ERROR: %s\n' "$1" >&2
    write_failure 1
}
trap fail ERR

for path in "$RUN_ROOT" "$RUN_LIBRARY" "$WORKFLOW_REPO" "$SVCFIT_PKG_DIR" "$ENV_SVCFIT"; do
    [[ -d "$path" ]] || die "required directory is absent: $path"
done
for path in "$RSCRIPT" "$ENV_SVCFIT/bin/R" "$RETICULATE_PYTHON"; do
    [[ -x "$path" ]] || die "required executable is absent: $path"
done

observed_workflow_commit=$(git -C "$WORKFLOW_REPO" rev-parse HEAD)
observed_svcfit_commit=$(git -C "$SVCFIT_PKG_DIR" rev-parse HEAD)
[[ "$observed_workflow_commit" == "$EXPECTED_WORKFLOW_COMMIT" ]] || {
    die "workflow commit is $observed_workflow_commit; expected $EXPECTED_WORKFLOW_COMMIT"
}
[[ "$observed_svcfit_commit" == "$EXPECTED_SVCFIT_COMMIT" ]] || {
    die "SVCFit commit is $observed_svcfit_commit; expected $EXPECTED_SVCFIT_COMMIT"
}
[[ -z "$(git -C "$WORKFLOW_REPO" status --porcelain)" ]] || {
    die 'workflow checkout is dirty'
}
[[ -z "$(git -C "$SVCFIT_PKG_DIR" status --porcelain)" ]] || {
    die 'SVCFit checkout is dirty'
}

mkdir -p "$RUN_ROOT/cache" "$RUN_ROOT/r-user" "$RUN_ROOT/tmp"
export R_LIBS_USER="$RUN_LIBRARY:$RLIB"
export R_USER="$RUN_ROOT/r-user"
export XDG_CACHE_HOME="$RUN_ROOT/cache"
export TMPDIR="$RUN_ROOT/tmp"

{
    printf 'host\t%s\n' "$(hostname -f)"
    printf 'started_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-not-slurm}"
    printf 'workflow_commit\t%s\n' "$observed_workflow_commit"
    printf 'svcfit_commit\t%s\n' "$observed_svcfit_commit"
    printf 'rscript\t%s\n' "$RSCRIPT"
    printf 'r_library\t%s\n' "$RUN_LIBRARY"
    printf 'reticulate_python\t%s\n' "$RETICULATE_PYTHON"
} > "$RUN_ROOT/provenance/worker-environment.tsv"

"$ENV_SVCFIT/bin/R" CMD INSTALL --install-tests --library="$RUN_LIBRARY" \
    "$SVCFIT_PKG_DIR" > "$RUN_ROOT/install.log" 2>&1

RUN_LIBRARY="$RUN_LIBRARY" "$RSCRIPT" -e '
run_lib <- normalizePath(Sys.getenv("RUN_LIBRARY"), mustWork=TRUE)
installed <- normalizePath(find.package("SVCFit"), mustWork=TRUE)
if (!startsWith(installed, paste0(run_lib, .Platform$file.sep)))
    stop("SVCFit was not loaded from the run-local library: ", installed)
testthat::test_package(
    "SVCFit",
    reporter="summary",
    stop_on_failure=TRUE,
    stop_on_warning=FALSE
)
' > "$RUN_ROOT/testthat.log" 2>&1

"$RSCRIPT" -e 'sessionInfo()' > "$RUN_ROOT/provenance/sessionInfo.txt" 2>&1
(
    cd "$RUN_ROOT"
    sha256sum install.log testthat.log provenance/SVCFit.commit.txt \
        provenance/workflow.commit.txt provenance/worker-environment.tsv \
        > provenance/output-files.sha256
)

rm -f "$RUN_ROOT/FAIL"
printf 'passed_at=%s\nsvcfit_commit=%s\nscope=installed-package-tests\n' \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$observed_svcfit_commit" > "$RUN_ROOT/PASS"
trap - ERR
printf 'SVCFit software self-evaluation: PASS\n'
