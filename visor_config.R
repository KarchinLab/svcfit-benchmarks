## R-side loader for config.local.sh. The counterpart to the walk-up block in
## chrx_common.sh / pipeline_common.sh, and it resolves the config the same way,
## so R and bash entry points cannot disagree about where a path points.
##
## WHY THIS EXISTS SEPARATELY. The shell scripts get their paths by sourcing
## config.local.sh directly. An .Rmd knitted from RStudio or `Rscript -e
## rmarkdown::render(...)` has no such shell, so Sys.getenv() alone returns ""
## and the document would fall back to whatever the author last hardcoded. That
## is the failure this repo exists to prevent, so this file sources the config
## in a subshell and imports the result rather than trusting the environment.
##
## Usage:
##   source("<repo>/visor_config.R")
##   cfg <- visor_config(c("REPLICATES_DIR", "MC_HACK_BASE"))
##   cfg$REPLICATES_DIR

## Walk up from `start` looking for config.local.sh. $VISOR_CONFIG overrides,
## matching the bash side (used when a SLURM task is copied to a node-local
## spool directory and the tree is no longer above it).
.visor_find_config <- function(start = NULL) {
    override <- Sys.getenv("VISOR_CONFIG", unset = "")
    if (nzchar(override)) {
        if (!file.exists(override))
            stop("VISOR_CONFIG is set to '", override, "', which does not exist.", call. = FALSE)
        return(normalizePath(override))
    }

    if (is.null(start)) start <- .visor_this_dir()
    d <- normalizePath(start, mustWork = FALSE)
    repeat {
        cand <- file.path(d, "config.local.sh")
        if (file.exists(cand)) return(normalizePath(cand))
        parent <- dirname(d)
        if (parent == d) break          # reached /
        d <- parent
    }
    stop("config.local.sh not found by walking up from '", start, "'.\n",
         "       cp config.example.sh config.local.sh at the repo root and edit it.\n",
         "       See README.md.", call. = FALSE)
}

## Best guess at "the directory of the file being executed", which R does not
## expose uniformly. knitr sets its own root during a render; otherwise fall
## back to the working directory. Only ever used as a starting point for the
## walk up, so being one directory off is harmless as long as it is inside the
## tree.
.visor_this_dir <- function() {
    if (requireNamespace("knitr", quietly = TRUE) && !is.null(knitr::current_input())) {
        inp <- knitr::current_input(dir = TRUE)
        if (!is.null(inp) && nzchar(inp)) return(dirname(normalizePath(inp, mustWork = FALSE)))
    }
    getwd()
}

## Source config.local.sh in a subshell and return the named variables.
##
## Sourced rather than parsed: the file is real bash, including interpolation,
## a _visor_require assert that must be allowed to fire, and a PATH guard.
## set_libpaths: prepend $RLIB to .libPaths() as a side effect.
##
## On by default, and that is the point. RLIB is where this project's R packages
## actually live -- the CRAN dependencies and the figure-only packages both -- and
## it is NOT the conda env's library, so an R process that has not been told about
## it cannot see any of them. A bare Rscript, or an .Rmd knitted from RStudio,
## would otherwise report installed packages as missing.
visor_config <- function(vars, start = NULL, set_libpaths = TRUE) {
    stopifnot(is.character(vars), length(vars) > 0)
    cfg <- .visor_find_config(start)

    ## $1 is the config, the rest are variable names. Sourcing writes nothing to
    ## stdout on success (verified: `bash -c 'source ./config.local.sh'` is
    ## silent), but redirect anyway so a future `echo` in the config cannot be
    ## mistaken for a value. Errors stay on stderr and surface below.
    ## RLIB is always queried, whether or not the caller asked for it, so
    ## set_libpaths works without every call site having to know to request it.
    ## It is dropped from the return value below unless it was actually asked for.
    query <- unique(c(vars, if (set_libpaths) "RLIB"))

    script <- 'source "$1" >/dev/null || exit 1
               shift
               for v in "$@"; do printf "%s\t%s\n" "$v" "${!v-}"; done'
    out <- suppressWarnings(
        system2("bash", c("-c", shQuote(script), "_", shQuote(cfg), shQuote(query)),
                stdout = TRUE, stderr = TRUE)
    )
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0)
        stop("sourcing ", cfg, " failed:\n", paste(out, collapse = "\n"), call. = FALSE)

    kv <- strsplit(out[grepl("\t", out)], "\t", fixed = TRUE)
    res <- setNames(vapply(kv, function(p) if (length(p) > 1) p[2] else "", character(1)),
                    vapply(kv, `[`, character(1), 1))

    ## Hard-fail on anything unset or missing, for the reason config.local.sh
    ## gives in its own assert: a default is what turns "wrong machine" into a
    ## plausible wrong number instead of an error.
    missing <- vars[!vars %in% names(res) | !nzchar(res[vars])]
    if (length(missing))
        stop("config.local.sh does not set: ", paste(missing, collapse = ", "),
             "\n       Add them to ", cfg, " (see config.example.sh).", call. = FALSE)
    absent <- vars[!file.exists(res[vars])]
    if (length(absent))
        stop("config.local.sh points at paths that do not exist on this machine:\n",
             paste0("       ", absent, " -> ", res[absent], collapse = "\n"), call. = FALSE)

    ## Prepend, not replace: the conda env's own library must stay reachable, or
    ## base packages and anything installed there disappears. Silent when RLIB is
    ## unset or absent -- this is a convenience, and a machine that keeps its
    ## packages in the env default is a legitimate setup. What is NOT acceptable
    ## is the reverse, a library that exists and is never consulted.
    if (set_libpaths) {
        rlib <- res[["RLIB"]]
        if (!is.null(rlib) && nzchar(rlib) && dir.exists(rlib) && !rlib %in% .libPaths())
            .libPaths(c(rlib, .libPaths()))
    }

    as.list(res[vars])
}
