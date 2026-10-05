# FACETS setup

The calling workflows require three FACETS components:

1. the `facets` R package;
2. the `snp-pileup` executable used to count alleles from normal/tumour BAMs;
3. a command-line R driver that turns the count matrix into `.bed` and `.RData`
   outputs with the historical workflow contract.

The repository now supplies a reproducible path for all three. It includes a
portable conda environment specification, the project-specific driver at
`tools/facet.R`, and an installer that compiles the official MSKCC
`snp-pileup` source shipped with the R package.

## Install

From the repository root:

```bash
tools/setup_facets.sh
```

The script creates the `facet` conda environment when needed, builds
`snp-pileup`, installs a runnable copy of the tracked `tools/facet.R` driver,
and prints the exact `config.local.sh` settings. By default these installed
tools are placed beside this checkout under `01_software/facets`, not committed
to Git.

To choose another in-layout destination or rebuild after an operating-system
or conda change:

```bash
tools/setup_facets.sh --dest /path/to/svcfit/01_software/facets
tools/setup_facets.sh --rebuild
```

Then run:

```bash
./check_deps.sh
```

The dependency check fails clearly if the environment, R package, executable,
or bundled driver is missing.

## Upstream provenance

- FACETS R package: <https://github.com/mskcc/facets>
- Official `snp-pileup` build instructions:
  <https://github.com/mskcc/facets/blob/master/inst/extcode/README.txt>
- Optional modern wrapper suite: <https://github.com/mskcc/facets-suite>

The workflow retains its established FACETS parameters and output format. It
does not switch to `facets-suite`; that project is linked only as an optional
upstream wrapper. The repository-owned driver makes the exact analysis adapter
auditable and prevents a hidden dependency on `~/facet/facet.R`.
