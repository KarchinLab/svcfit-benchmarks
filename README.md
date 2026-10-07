# SVCFit benchmark workflows

Workflows that produced the benchmark results in:

Liu Y, Lai J, Yang Y, Markowski MC, Antonarakis ES, De Marzo AM, Yegnasubramanian S, Wood LD, Sena LA, Karchin R. Longitudinal structural variant phylogenies define tumor evolution under therapeutic selection pressure in metastatic prostate cancer. npj Precision Oncology (in revision).

The SVCFit R package is a separate repository, https://github.com/KarchinLab/SVCFit. All reported results correspond to SVCFit commit `7f32d81f3dd0eee0f2b8623e2775aae70ce0e917`.

This repository covers four benchmarks, from simulation to the event-level outputs:

| Folder | Benchmark | Paper |
|---|---|---|
| `visor_replicates/` | Autosomal VISOR accuracy benchmark: 75 scenarios, 5 purities, 30 read replicates; SVCFit, SVclone and assisted SVclone | Figure 3A and 3B, Supplementary Figures S3 to S7 and S12, Tables S7 to S14 |
| `visor_chrX/` | Hemizygous chromosome X simulation: 45 conditions, 30 read replicates | Figure 3B (hemizygous), Supplementary Note S9.4 |
| `prostate_mixture/` | Prostate cancer mixtures: 11 mixtures, 30 read replicates | Figure 3C, Table S15 |
| `tree_eval/eval_package/Phylogeny_benchmark/` | Longitudinal phylogeny simulations: 100 planned cases (S1 scenario) | Supplementary Note S3, Figures S8 to S11, Tables S3 and S4 |

The figures and the statistics tables are regenerated from these event-level outputs by the scripts in the Mendeley Data deposit (Data and scripts for "Longitudinal structural variant phylogenies...", Version 7, DOI to be added). That deposit also contains the event-level outputs themselves, so the paper's numbers can be checked without rerunning these workflows.

The clinical (COMBAT) analysis is not part of this repository.

## Requirements

- Linux with SLURM. All stages are written as SLURM jobs or arrays; set account, partition and QoS through `config.local.sh` (`SLURM_ACCOUNT`, `SLURM_PARTITION`, `SLURM_QOS`) or your site defaults.
- Conda. Environment specifications are in `tools/environments/` and `prostate_mixture/prostate_replicates/conda_envs/`. `tools/setup_rockfish_tool_envs.sh` builds them; `tools/setup_rockfish_runtime.sh` builds the core runtime (R 4.4.3, SVCFit, scikit-learn through reticulate). The scripts are named for the cluster the runs used and take every path from configuration.
- FACETS: run `tools/setup_facets.sh` (see `FACETS_SETUP.md`).
- Tool versions are listed in Supplementary Table S1 of the paper (VISOR 1.1.3, Manta 1.6.0, SVtyper 0.7.1, GATK 4.6.2.0, FACETS 0.6.2, DNAcopy 1.80.0, SVclone 1.1.2, scikit-learn 1.9.0).
- An SVCFit checkout at commit `7f32d81` (`SVCFIT_PKG_DIR`) and an installed SVCFit library built from it. `tools/run_svcfit_self_evaluation.sh` installs the library into a fresh directory, runs the package tests and records the commit in `provenance/SVCFit.commit.txt`. The autosomal VISOR stage loads SVCFit from the checkout (`devtools::load_all`); the chromosome X, prostate and phylogeny stages load the installed library (`SVCFIT_R_LIB` with `SVCFIT_R_LIB_COMMIT_FILE`, or `--rlib` with `--rlib-commit-file` for prostate) and refuse one built from another commit. The phylogeny launcher only warns if `SVCFIT_R_LIB` is unset and then uses whatever SVCFit is on `RLIB`, so set it.

## Configuration

```bash
cp config.example.sh config.local.sh   # edit PROJECT_ROOT, SVCFIT_RUNTIME_LOADER and the data paths
./check_deps.sh
```

Every script finds `config.local.sh` by walking up from its own location, or from `$VISOR_CONFIG`. `visor_config.R` is the R-side reader. The expected layout under `PROJECT_ROOT` is `01_software/` (this repository and SVCFit), `02_data/` (references and simulated data), `03_analysis/` (runs) and `04_qc/`.

Set `VISOR_ROOT` to this checkout; the default in `config.example.sh` names an older directory (`svcfit_workflows`), and the prostate stage sources its helper scripts from `$VISOR_ROOT`.

Several submit scripts record provenance and refuse to run from a modified checkout. They compare `git rev-parse HEAD` with `EXPECTED_WORKFLOW_COMMIT` and `EXPECTED_SVCFIT_COMMIT`; set the first to the commit of this repository you are running and the second to `7f32d81f3dd0eee0f2b8623e2775aae70ce0e917`. The autosomal VISOR and prostate submit scripts read these from the environment and do not source `config.local.sh` themselves, so export them (or `source config.local.sh`) first. Submit scripts are dry runs unless given `--submit` (the phylogeny launcher submits unless given `--dry-run`).

## Input data

- Reference genomes: chromosomes 1 and 2 (`REF_AUTO`) for the autosomal and phylogeny simulations, GRCh38 chromosomes 22 and X (`REF_CHRX`) for the chromosome X simulation, and hs37d5 (`PROSTATE_REF`) for the prostate mixtures.
- Chromosome X simulation design: `visor_chrX/resources/beds/` (planted SVs and copy-number changes, excluded regions, depth windows) and `visor_chrX/truth/` (condition table, per-clone SV BED files and the ground-truth table).
- Autosomal and phylogeny simulation design (both benchmarks use the same planted SVs), in the Mendeley deposit:
  - `VISOR_benchmark.zip`: `ground_truth/sv_beds/` (per-clone SV BED files; the scoring truth and the HACk input, `MC_HACK_BASE`), `input_data/beds/` (haplotype, SV and germline SNP BED files), `input_data/snp_vcfs/` (germline SNP VCFs for chromosomes 1 and 2, from dbSNP build 157) and `depth/HACk.random.bed`;
  - `Phylogeny_benchmark.zip`: `input_data/hack/` (`TREE_EVAL_TRUTH_DIR`), `input_data/resource/`, `input_data/norm_short/` and `input_data/reference/chr1-2.fa` (GRCh38 chromosomes 1 and 2, with `.fai` and `.dict`; build the BWA index with `bwa index`).
- The prostate mixtures are built from the prostate cancer whole-genome sequencing data of Cmero et al. (Nat. Commun. 2020), which is controlled-access; request access from its data controller. This repository and the Mendeley deposit contain no reads, germline variants or other controlled data from those samples.

## Simulation design (`simulation_design/scripts/genome_setup/`)

These scripts generated the design files deposited in Mendeley; the deposited files are what the reported runs used, so steps 1 and 2 are optional. Step 3 is required to simulate the autosomal benchmark from scratch: the clone genomes it builds are not deposited, and `00_visor_shorts.sh` reads them from `MC_HACK_BASE` as `e1/` to `e5/` (each with `c2/` and `c3/` holding `h1.fa` and `h2.fa`), `snp1/` and `snp1_1/`, next to the per-clone BED files. `run_m_hack.sh` writes them under `truth/fastas/clone_genomes/`, so move or link them into `MC_HACK_BASE`.

1. SNP background: `download_snp.sh` (dbSNP build 157), `prepare_dbsnp_vcfs.sh`, then `filter_snp.sh` (heterozygous SNPs near the SVs) and `filter_snp_beds.sh` (SNPs within 1 kb of breakpoints).
2. SVs and copy-number changes: `make_sv.sh` (random 100-kb SVs on chromosome 1 with translocation partners on chromosome 2, using VISOR's `randomregion.r`; set `VISOR_HOME`), `modify_sv.R`, `make_cnv.R`, and `sv2subclone.R` (per-clone BED files).
3. Clone genomes: `run_hack.sh` and `run_m_hack.sh` (VISOR HACk for the five SV-CNV configurations).

The scripts expect `resources/beds/`, `resources/snp_vcfs/` and `truth/sv_beds/` under `simulation_design/`; these correspond to `input_data/beds/`, `input_data/snp_vcfs/` and `ground_truth/sv_beds/` in the Mendeley `VISOR_benchmark.zip`. Set `REF` (or `REF_DIR` for `make_sv.sh`) to the GRCh38 chromosome 1 and 2 FASTA.

## Autosomal VISOR benchmark (`visor_replicates/`)

1. Simulate tumor BAMs: `00_visor_shorts.sh` (SLURM array of 2,250 tasks: 15 purity-mixture conditions x 5 configurations x 30 replicates; replicate-specific seeds). It needs the HACk clone genomes in `MC_HACK_BASE` (see Simulation design) and the matched normal BAMs `NORM_SHORT_DIR/normal/sim.srt.bam` and `NORM_SHORT_DIR/o_normal/norm.bam`. The 10% purity conditions (tasks with `task_id % 15` in 0 to 2) were simulated separately by `submit_10pct.sh`, which runs stages 00 to 03 for p10 only with its own seed offset; do not also run them through `00_visor_shorts.sh`, which would give them different seeds.
2. Call and characterize: `submit_all.sh` submits `01_manta_svtyper.sh` (Manta and SVtyper), `02_snp_pipeline.sh` (germline heterozygous SNPs and phasing) and `03_facet.sh` (FACETS) with dependencies. It checks for the 2,250 BAMs before submitting anything. It does not submit `04_svclone.sh`; the SVclone arms are step 3.
3. SVclone, two arms, from the existing calls:
   - assisted SVclone (true cellular fraction resolves multiplicity): `submit_downstream_only.sh`;
   - SVclone (FACETS purity and copy number, no truth): `submit_fair_downstream_only.sh`.
4. Build the condition ledger and extract the SVclone estimates: `Rscript 06_build_three_arm_manifest.R --help` and `Rscript 07_extract_svclone_arms.R --help` give the arguments.
5. SVCFit: `submit_three_arm_svcfit.sh --manifest <root>/three_arm_condition_manifest.tsv --output-root <root> --rscript "$SVCFIT_R" --svcfit-package "$SVCFIT_PKG_DIR" --truth-dir "$MC_HACK_BASE" --rlib <R library with devtools> --array 0-2249%25 --submit` (`08_svcfit_three_arm_array.sh`; the default `--array` is a four-task test). Use one analysis root `<root>` for steps 4 to 6: it holds `three_arm_condition_manifest.tsv` (step 4), `svclone_native_events_long.rds` (step 4) and `svcfit_events/` (this step).
6. Join, summarize and compare, with `--analysis-root <root>`:
   - `09_join_three_arms.R` (default tolerance), then `09_join_three_arms.R --tolerance 25 --output-prefix primary_tol025_ --rds-only` (the primary matching that 11 reads), then `--summary-only` runs at tolerance 0, 25 and 50 with prefixes `sensitivity_tol000_`, `sensitivity_tol025_` and `sensitivity_tol050_`;
   - `10_summarize_three_arms.R`, then `11_compare_three_arms.R`.

   09 caches the SVCFit events in `<root>/svcfit_events_long.rds` and reuses the cache if it exists. Never copy that file into a new root, or 09 will keep the old SVCFit values.

## Chromosome X simulation (`visor_chrX/scripts/`)

1. Build the haplotypes (chromosome 22 diploid, chromosome X single-copy) and plant the SVs: `01_chrx_hack.sh`; ground truth: `04_chrx_ground_truth.sh` (the condition table and truth table it writes are also in `visor_chrX/truth/`).
2. Build the shared matched normal once: `sbatch --export=ALL,COV=50 scripts/02_chrx_shorts.sh normal` (a single job, not an array). Every replicate links to it, and `15_chrx_replicate_chain.sh` refuses to start without it.
3. Simulate, segment, call and score all 30 replicates at 50x: `COV=50 16_chrx_replicates_all.sh` submits `15_chrx_replicate_chain.sh` for replicates 1 to 30 (`02_chrx_shorts.sh`, then `12_chrx_calls_array.sbatch` and `06_chrx_segment_array.sbatch`, then `17_chrx_score_one.sbatch`, which runs the 09 join and step 13). 17 takes SVCFit from `SVCFIT_R_LIB`, so set it before submitting. `18_chrx_replicate_sweep.sh` reconciles replicates against what is on disk and resubmits gaps.
4. SVclone, per replicate: `21_chrx_svclone.sbatch` (an array over the 45 conditions; it calls `20_chrx_svclone_input.R`), then `22_chrx_score_svclone.R` (with `COV` and `REP` set), which writes `scoring_rep<N>/chrx_svclone_scores_c50.tsv`.
5. Summarize: `19_chrx_bootstrap.R` (SVCFit) and `23_chrx_compare_bootstrap.R` (SVCFit against SVclone), with `CHRX_DIR` and `COV=50` set. To rescore SVCFit alone against a pinned library without touching the step 3 outputs, use `25_chrx_rescore_svcfit.sh --submit` (needs `EXPECTED_WORKFLOW_COMMIT`, `EXPECTED_SVCFIT_COMMIT`, `SVCFIT_R_LIB` and `SVCFIT_R_LIB_COMMIT_FILE`; `RUN_ROOT` and `N_REPS` are optional). It reruns step 13 for every replicate into a new root, then 19 and 23 there. Its downstream job also copies `chrx_svclone_fair_scores_c50.tsv` from each `scoring_rep<N>/`, and no script in this repository writes that file.

Set `COV=50` for every stage; it is required. The scripts read the design files from `$CHRX_DIR`, so copy `visor_chrX/resources/` and `visor_chrX/truth/` into `$CHRX_DIR` before step 1. `15` and `18` change to `$CHRX_DIR` and submit `scripts/<name>`, so also link the scripts there: `ln -s "$PWD/visor_chrX/scripts" "$CHRX_DIR/scripts"`.

## Prostate mixtures (`prostate_mixture/scripts/`)

1. `submit_all.sh`: source-BAM filtering (`00a`), mixture construction (`00`), Manta and SVtyper (`01`), germline SNPs (`02`, `03`), FACETS (`04`) and assisted SVclone (`05`).
2. SVclone (FACETS inputs, no truth): `submit_svclone_robust_only.sh`.
3. Chromosome X copy number from depth: `04a_chrx_depth_segmentation.sh` (SLURM array, 330 tasks; it runs `helper/chrx_depth_segmentation.R`), then, after it finishes, the segment-to-SV join `04b_chrx_segment_sv_join.sh` (also 330 tasks). 4m and 5m have no chromosome X and are skipped.
4. SVCFit: `submit_svcfit_correction.sh --output-root <out> --data-root "$PROSTATE_DATA_DIR" --rscript "$SVCFIT_R" --svcfit-package "$SVCFIT_PKG_DIR" --rlib <SVCFit library> --rlib-commit-file <its SVCFit.commit.txt> --array 0-29%10 --submit` (`09_svcfit_correction_array.sh` runs `06a_run_svcfit_chrx.R` for the nine three-cluster mixtures and for 4m/5m, writing `parts/svcfit_chrx_rep<N>.tsv` and `parts/svcfit_45_rep<N>.tsv`). Submit it from inside this checkout: `06a` evaluates the `setup` and `helpers` chunks of `06_svcfit_replicates_shared_sv.Rmd`, which finds `visor_config.R` by walking up from the working directory.
5. Combine the parts into `<out>/svcfit_chrx_all.tsv` and `<out>/svcfit_45_all.tsv`: the header once, then replicates 1 to 30 in numeric order without headers.
6. Three-way comparison on the shared event set: `08_prostate_three_arm.R --svcfit-chrx <out>/svcfit_chrx_all.tsv --svcfit-45 <out>/svcfit_45_all.tsv --assisted-root ... --fair-root ... --truth-dir ... --output-dir ... --svclone-rlib ... --bootstrap 2000 --seed 20260920` (`--assisted-root` is the step 1 replicate tree, `--fair-root` the step 2 result root).

## Phylogeny benchmark (`tree_eval/eval_package/Phylogeny_benchmark/scripts/`)

1. Simulate and call: `run_all.sh` submits `longi_short.sh` (VISOR SHORtS, two timepoints) and `longi_calling.sh` (Manta, SVtyper, SNPs, FACETS) for scenarios S1 to S4, purities 10% to 80% and simulation replicates BOOT 0 to 4. Only S1 at 20% to 80% is reported. These stages read the HACk files from `Phylogeny_benchmark/data/hack/` and write to `Phylogeny_benchmark/outputs/`, both inside the checkout, not from `TREE_EVAL_TRUTH_DIR` or to `TREE_EVAL_LONGITUDINAL`. Put `input_data/hack/` there before running, and move or link `outputs/` to `TREE_EVAL_LONGITUDINAL` afterwards. The `longi_svcfit.sh` and evaluation jobs that `run_all.sh` also submits do not get the variables `longi_svcfit.sh` requires (`OUTPUT_DIR`, `INPUT_DIR`, `SVCFIT_REPO`, ...), so they fail; step 2 replaces them.
2. Reported results: `run_svcfit_and_evaluate.sh --scope full` re-runs SVCFit, clustering and tree reconstruction at SVCFit `7f32d81` on the S1 simulations in `TREE_EVAL_LONGITUDINAL` (20% to 80% purity, 5 configurations, 5 replicates; 100 cases) and evaluates them with `evaluate_downstream.Rmd`. It needs `EXPECTED_WORKFLOW_COMMIT` plus `SVCFIT_R_LIB` and `SVCFIT_R_LIB_COMMIT_FILE`, and writes to `03_analysis/tree_eval/runs/<id>/` unless `RUN_ROOT` is set. `--scope smoke` runs one case (S1, 20% purity, BOOT 0, configuration 1); `--dry-run` checks without submitting.
3. Covariate analyses (Supplementary Note S3.5): `A2_coverage_correlation.R` and `a2_coverage_correlation/`.

## License

MIT; see `LICENSE`. Third-party tools are installed from their own distributions under their own licenses.

## Citation

See `CITATION.cff`.
