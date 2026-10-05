#!/bin/bash
# Manta SV calling + SVtyper genotyping for ONE chrX simulation condition.
#
# Supplies sv_alt (AO) and sv_ref (RO) -- the read counts calc_svcf() turns into VAF.
#
# This mirrors visor_replicates/01_manta_svtyper.sh step for step, deliberately,
# so chrX numbers stay comparable to the autosomal ones.
# Differences from that script, all forced and all recorded:
#
#   - absolute environment prefixes come from config.local.sh
#     (ENV_MANTA / ENV_SVTYPER). Same Manta 1.6.0, same SVtyper 0.7.1.
#   - reference is chr22-X.fa, not chr1-2.fa.
#   - one shared normal BAM (short/normal), since the simulation has exactly one.
#
# Kept identical on purpose, because they change the numbers:
#   - somatic mode: --normalBam AND --tumorBam
#   - convertInversion.py over somaticSV.vcf.gz -- without it Manta reports
#     inversions as BND pairs and parse_sv_info's INV branch never fires. 2,190
#     of the 4,935 simulated SVs are inversions, so this is not optional here.
#   - the CIPOS/CIEND = +/-100 rewrite before SVtyper
#
# Usage: ./11_chrx_manta_svtyper.sh <exp> <cond> [threads]
#        ./11_chrx_manta_svtyper.sh e1 c25p80m50 4

set -euo pipefail

exp="${1:?usage: $0 <exp> <cond> [threads]}"
cond="${2:?usage: $0 <exp> <cond> [threads]}"
threads="${3:-4}"

_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_DIR}/chrx_common.sh"      # resolves BASE_DIR, REF and the conda env names
                                     # from config.local.sh -- no literal paths here
tbam="${SHORT_DIR}/${exp}/${cond}/short.out/sim.srt.bam"
nbam="${SHORT_DIR}/normal_c${cov}/short.out/sim.srt.bam"   # cov-named; see 02_chrx_shorts.sh
name="${exp}_${cond}"

# cond is a free-form argument here, but the normal is chosen by $COV. Called
# directly as `11_... e1 c50p80m50` with COV unset, this would genotype a 50x
# tumour against the 25x normal: cn_bar is a tumour-to-normal depth ratio, so
# every cellular fraction on that contig would come out wrong by a factor of two
# with nothing looking broken. Refuse instead.
cond_cov="${cond#c}"; cond_cov="${cond_cov%%p*}"
[[ "$cond_cov" == "$cov" ]] || {
    echo "ERROR: condition '${cond}' is ${cond_cov}x but COV=${cov}, so the normal would be" >&2
    echo "       normal_c${cov} -- a tumour/normal coverage mismatch. Re-run with COV=${cond_cov}." >&2
    exit 1; }

# chrx_common.sh already rep-tags $CALLS_DIR; an explicit OUT_DIR is not tagged,
# so tag that one here. Branching keeps the tag from being applied twice.
if [[ -n "${OUT_DIR:-}" ]]; then
    out="${OUT_DIR}${rep_tag}/${name}"
else
    out="${CALLS_DIR}/${name}"
fi
manta_dir="${out}/manta"
svtyp_dir="${out}/svtyp"

for f in "$REF" "$REF.fai" "$tbam" "$nbam"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
mkdir -p "$manta_dir" "$svtyp_dir"

source "$CONDA_SH"

# ---- Manta -----------------------------------------------------------------
set +u   # the manta env's openjdk activate script reads $JAVA_HOME unguarded
conda activate "$ENV_MANTA"
set -u
echo "[$(date +%T)] manta $($MANTA_CONFIG --version 2>&1 | head -1)  ${name}"

"$ENV_MANTA/bin/python" "$MANTA_CONFIG" \
    --normalBam "$nbam" \
    --tumorBam  "$tbam" \
    --referenceFasta "$REF" \
    --runDir "$manta_dir" > "${out}/configManta.log" 2>&1

"$ENV_MANTA/bin/python" "$manta_dir/runWorkflow.py" -j "$threads" > "${out}/runWorkflow.log" 2>&1

# Manta writes this, and convertInversion.py below reads it. Asserted so that a
# missing input fails here rather than reading as an SVtyper problem.
[[ -s "$manta_dir/results/variants/somaticSV.vcf.gz" ]] || {
    echo "ERROR: ${name}: Manta produced no somaticSV.vcf.gz; see ${out}/runWorkflow.log" >&2
    exit 1; }

# convertInversion.py needs samtools; it lives in the visor env, not on PATH here.
samtools_bin="$SAMTOOLS"
convert_py="$MANTA_CONVERT"
[[ -r "$convert_py" ]] || {
    echo "ERROR: ${name}: convertInversion.py not found at $convert_py" >&2
    echo "       configured Manta prefix: ${ENV_MANTA}" >&2
    exit 1; }

"$ENV_MANTA/bin/python" "$convert_py" \
    "$samtools_bin" "$REF" \
    "$manta_dir/results/variants/somaticSV.vcf.gz" \
    > "$manta_dir/results/variants/${name}.vcf"

# convertInversion is what makes Manta's inversions parseable at all -- 2,190 of
# the 4,935 simulated SVs are inversions -- so an empty result here is not a
# quiet degradation, it removes half the truth set from scoring.
[[ -s "$manta_dir/results/variants/${name}.vcf" ]] || {
    echo "ERROR: ${name}: convertInversion.py produced no output" >&2; exit 1; }

conda deactivate

# ---- SVtyper ---------------------------------------------------------------
set +u; conda activate "$ENV_SVTYPER"; set -u
echo "[$(date +%T)] svtyper $(svtyper --version 2>&1 | head -1)  ${name}"

# Manta's own CIPOS/CIEND are tight; the submitted pipeline widens them to
# +/-100 so SVtyper gathers evidence over a comparable window. Same rewrite.
awk 'BEGIN{FS=OFS="\t"}
     /^#/ { print; next }
     { info=$8
       gsub(/CIPOS=[^;]+;?/,"",info)
       gsub(/CIEND=[^;]+;?/,"",info)
       $8="CIPOS=-100,100;CIEND=-100,100;"info
       print }' \
    "$manta_dir/results/variants/${name}.vcf" \
    > "$svtyp_dir/${name}_w_CI.vcf"

svtyper \
    -i "$svtyp_dir/${name}_w_CI.vcf" \
    -B "$tbam" \
    -l "$svtyp_dir/${name}.bam.json" \
    > "$svtyp_dir/svt_${name}.vcf"

conda deactivate

# POST-CONDITION. A task can exit 0 without writing a VCF, so the output is asserted
# here; `grep -vc || true` below cannot detect a missing file on its own.
#
# So refuse to exit 0 without the artefact the caller will look for. An empty
# call set is a legitimate outcome for a low-purity condition, so the assertion
# is existence and a header, not a record count.
[[ -s "$svtyp_dir/svt_${name}.vcf" ]] || {
    echo "ERROR: ${name}: svtyper produced no $svtyp_dir/svt_${name}.vcf" >&2; exit 1; }
grep -q '^#' "$svtyp_dir/svt_${name}.vcf" || {
    echo "ERROR: ${name}: $svtyp_dir/svt_${name}.vcf has no VCF header -- truncated output" >&2
    exit 1; }

n=$(grep -vc '^#' "$svtyp_dir/svt_${name}.vcf" || true)
echo "[$(date +%T)] DONE ${name}: ${n} genotyped records -> $svtyp_dir/svt_${name}.vcf"
