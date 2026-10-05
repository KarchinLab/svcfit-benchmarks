# SVCFit benchmark workflows

Workflows that produced the benchmark results in:

Liu Y, Lai J, Yang Y, Markowski MC, Antonarakis ES, De Marzo AM, Yegnasubramanian S, Wood LD, Sena LA, Karchin R. Longitudinal structural variant phylogenies define tumor evolution under therapeutic selection pressure in metastatic prostate cancer. npj Precision Oncology (in revision).

The SVCFit R package is a separate repository, https://github.com/KarchinLab/SVCFit. All reported results use SVCFit commit `e0d7e0b8d704caa3cfb227bc3dbf1ee99d66fac7`.

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
- An installed SVCFit library at commit `e0d7e0b`. `tools/run_svcfit_self_evaluation.sh` installs it into a fresh library, runs the package tests and records the commit; later stages take that library through `SVCFIT_RLIB` / `SVCFIT_R_LIB` and refuse a library built from any other commit.

## Configuration

```bash
cp config.example.sh config.local.sh   # edit PROJECT_ROOT, SVCFIT_RUNTIME_LOADER and the data paths
./check_deps.sh
```

Every script finds `config.local.sh` by walking up from its own location, or from `$VISOR_CONFIG`. `visor_config.R` is the R-side reader. The expected layout under `PROJECT_ROOT` is `01_software/` (this repository and SVCFit), `02_data/` (references and simulated data), `03_analysis/` (runs) and `04_qc/`.

Several submit scripts record provenance and refuse to run from a modified checkout. They compare `git rev-parse HEAD` with `EXPECTED_WORKFLOW_COMMIT` and `EXPECTED_SVCFIT_COMMIT`; set the first to the commit of this repository you are running and the second to `e0d7e0b8d704caa3cfb227bc3dbf1ee99d66fac7`. Submit scripts are dry runs unless given `--submit`.

## Input data

- Reference genomes: chromosomes 1 and 2 (`REF_AUTO`) for the autosomal and phylogeny simulations, GRCh38 chromosomes 22 and X (`REF_CHRX`) for the chromosome X simulation, and hs37d5 (`PROSTATE_REF`) for the prostate mixtures.
- Chromosome X simulation design: `visor_chrX/resources/beds/` (planted SVs and copy-number changes, excluded regions, depth windows) and `visor_chrX/truth/` (condition table, per-clone SV BED files and the ground-truth table).
- Autosomal and phylogeny simulation design (both benchmarks use the same planted SVs), in the Mendeley deposit:
  - `VISOR_benchmark.zip`: `ground_truth/sv_beds/` (per-clone SV BED files; the scoring truth and the HACk input, `MC_HACK_BASE`), `input_data/beds/` (haplotype, SV and germline SNP BED files), `input_data/snp_vcfs/` (germline SNP VCFs for chromosomes 1 and 2, from dbSNP build 157) and `depth/HACk.random.bed`;
  - `Phylogeny_benchmark.zip`: `input_data/hack/` (`TREE_EVAL_TRUTH_DIR`), `input_data/resource/`, `input_data/norm_short/` and `input_data/reference/chr1-2.fa` (GRCh38 chromosomes 1 and 2, with `.fai` and `.dict`; build the BWA index with `bwa index`).
- The prostate mixtures are built from the prostate cancer whole-genome sequencing data of Cmero et al. (Nat. Commun. 2020), which is controlled-access; request access from its data controller. This repository and the Mendeley deposit contain no reads, germline variants or other controlled data from those samples.

## Simulation design (`simulation_design/scripts/genome_setup/`)

These scripts generated the design files deposited in Mendeley; the deposited files are what the reported runs used, so rerunning them is optional.

1. SNP background: `download_snp.sh` (dbSNP build 157), `prepare_dbsnp_vcfs.sh`, then `filter_snp.sh` (heterozygous SNPs near the SVs) and `filter_snp_beds.sh` (SNPs within 1 kb of breakpoints).
2. SVs and copy-number changes: `make_sv.sh` (random 100-kb SVs on chromosome 1 with translocation partners on chromosome 2, using VISOR's `randomregion.r`; set `VISOR_HOME`), `modify_sv.R`, `make_cnv.R`, and `sv2subclone.R` (per-clone BED files).
3. Clone genomes: `run_hack.sh` and `run_m_hack.sh` (VISOR HACk for the five SV-CNV configurations).

The scripts expect `resources/beds/`, `resources/snp_vcfs/` and `truth/sv_beds/` under `simulation_design/`; these correspond to `input_data/beds/`, `input_data/snp_vcfs/` and `ground_truth/sv_beds/` in the Mendeley `VISOR_benchmark.zip`. Set `REF` (or `REF_DIR` for `make_sv.sh`) to the GRCh38 chromosome 1 and 2 FASTA.

## Autosomal VISOR benchmark (`visor_replicates/`)

1. Simulate tumor BAMs: `00_visor_shorts.sh` (SLURM array over 75 conditions x 30 replicates; replicate-specific seeds). The 10% purity conditions use `submit_10pct.sh`.
2. Call and characterize: `submit_all.sh` submits `01_manta_svtyper.sh` (Manta and SVtyper), `02_snp_pipeline.sh` (germline heterozygous SNPs and phasing) and `03_facet.sh` (FACETS) with dependencies.
3. SVclone, two arms, from the existing calls:
   - assisted SVclone (true cellular fraction resolves multiplicity): `submit_downstream_only.sh`;
   - SVclone (FACETS purity and copy number, no truth): `submit_fair_downstream_only.sh`.
4. Build the condition ledger and extract the SVclone estimates: `Rscript 06_build_three_arm_manifest.R --help` and `Rscript 07_extract_svclone_arms.R --help` give the arguments.
5. SVCFit: `submit_three_arm_svcfit.sh --manifest ... --svcfit-package ... --rlib ... --submit` (2,250 tasks; `08_svcfit_three_arm_array.sh`).
6. Join, summarize and compare: `09_join_three_arms.R`, `10_summarize_three_arms.R`, `11_compare_three_arms.R` (each with `--help`).

## Chromosome X simulation (`visor_chrX/scripts/`)

1. Build the haplotypes (chromosome 22 diploid, chromosome X single-copy) and plant the SVs: `01_chrx_hack.sh`; ground truth: `04_chrx_ground_truth.sh`.
2. Simulate, segment, call and score all 30 replicates at 50x: `16_chrx_replicates_all.sh` submits `15_chrx_replicate_chain.sh` per replicate (`02_chrx_shorts.sh`, `06_chrx_segment_array.sbatch`, `12_chrx_calls_array.sbatch`, `17_chrx_score_one.sbatch`). `18_chrx_replicate_sweep.sh` reconciles replicates against what is on disk and resubmits gaps.
3. SVclone: `20_chrx_svclone_input.R`, `21_chrx_svclone.sbatch`, `22_chrx_score_svclone.R`.
4. Score SVCFit with the pinned library and summarize: `25_chrx_rescore_svcfit.sh --submit` runs step 13 for every replicate, then `19_chrx_bootstrap.R` and `23_chrx_compare_bootstrap.R`.

Set `COV=50` for every stage; it is required. The scripts read the design files from `$CHRX_DIR`, so copy `visor_chrX/resources/` and `visor_chrX/truth/` into `$CHRX_DIR` before step 1.

## Prostate mixtures (`prostate_mixture/scripts/`)

1. `submit_all.sh`: source-BAM filtering (`00a`), mixture construction (`00`), Manta and SVtyper (`01`), germline SNPs (`02`, `03`), FACETS (`04`) and assisted SVclone (`05`).
2. SVclone (FACETS inputs, no truth): `submit_svclone_robust_only.sh`.
3. Chromosome X copy number from depth: `04a_chrx_depth_segmentation.sh`, then the segment-to-SV join `04b_chrx_segment_sv_join.sh`.
4. SVCFit: `submit_svcfit_correction.sh --output-root ... --data-root ... --svcfit-package ... --rlib ... --submit` (30 tasks; `09_svcfit_correction_array.sh` runs `06a_run_svcfit_chrx.R` for the nine three-cluster mixtures and for 4m/5m).
5. Three-way comparison on the shared event set: `Rscript 08_prostate_three_arm.R --help`.

## Phylogeny benchmark (`tree_eval/eval_package/Phylogeny_benchmark/scripts/`)

1. Simulate and call: `run_all.sh` submits `longi_short.sh` (VISOR SHORtS, two timepoints), `longi_calling.sh` (Manta, SVtyper, SNPs, FACETS) and `longi_svcfit.sh` for every scenario, purity and simulation replicate (BOOT 0 to 4).
2. Reported results: `run_svcfit_and_evaluate.sh --scope full` re-runs SVCFit, clustering and tree reconstruction at SVCFit `e0d7e0b` on the S1 simulations (20% to 80% purity, 5 configurations, 5 replicates; 100 cases) and evaluates them with `evaluate_downstream.Rmd`. `--scope smoke` runs one case.
3. Covariate analyses (Supplementary Note S3.5): `A2_coverage_correlation.R` and `a2_coverage_correlation/`.

## License

MIT; see `LICENSE`. Third-party tools are installed from their own distributions under their own licenses.

## Citation

See `CITATION.cff`.
