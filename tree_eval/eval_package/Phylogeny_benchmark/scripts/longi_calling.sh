#!/bin/bash
#SBATCH --job-name=longi_pipe
#SBATCH --output=log/longi_pipe_%A_%a.out
#SBATCH --error=log/longi_pipe_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=3G
#SBATCH --time=12:00:00
#SBATCH --array=0-4

###############################################################################
# run_longi_pipeline.sh
#
# Runs Manta -> SVtyper -> SNP calling -> phasing -> FACETS on BOTH timepoints
# of a longitudinal simulation scenario.
#
# Depends on: run_longi_short.sh having completed successfully.
#
# Usage:
#   Set PPUR, SCENARIO, COV below, then submit:
#     sbatch run_longi_pipeline.sh
###############################################################################

set -e

thread=$SLURM_CPUS_PER_TASK
cov=50

########## DEFAULTS (override via --export from run_all.sh) ##########
ppur=${ppur:-60}
SCENARIO=${SCENARIO:-S1}
BOOT=${BOOT:-0}
######################################################################

exp_lst=(exp1 exp2 exp3 exp4 exp5)
nor_lst=(normal o_normal o_normal o_normal o_normal)
# normal uses sim.srt.bam; o_normal uses norm.bam (matches visor_replicates/pipeline_common.sh)
nor_bam=(sim.srt.bam norm.bam norm.bam norm.bam norm.bam)
samp_name=${exp_lst[$SLURM_ARRAY_TASK_ID]}
normal_name=${nor_lst[$SLURM_ARRAY_TASK_ID]}

# --- PATHS ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
# shellcheck disable=SC1091
source "${VISOR_CONFIG:-$REPO_ROOT/config.local.sh}"
PUB_DIR="$(dirname "$SCRIPT_DIR")"
out_base=$PUB_DIR/outputs
ref=$PUB_DIR/data/reference/chr1-2.fa
samtool_path="$SAMTOOLS"
scrp_dir=$SCRIPT_DIR

if [[ $BOOT -eq 0 ]]; then
  longi_dir=$out_base/${SCENARIO}/c${cov}p${ppur}/$samp_name
else
  longi_dir=$out_base/${SCENARIO}/c${cov}p${ppur}/b${BOOT}/$samp_name
fi
nbam=$PUB_DIR/data/norm_short/$normal_name/${nor_bam[$SLURM_ARRAY_TASK_ID]}

source "$CONDA_SH"

[[ ! -s $nbam ]] && { echo "ERROR: normal BAM missing or empty: $nbam" >&2; exit 1; }

# Process each timepoint
for tp in t1 t2; do

    tbam=$longi_dir/$tp/short.out/sim.srt.bam
    [[ ! -s $tbam ]] && { echo "ERROR: tumor BAM missing or empty for $tp: $tbam" >&2; exit 1; }
    manta_dir=$longi_dir/$tp/manta
    svtyp_dir=$longi_dir/$tp/svtyp
    snp_dir=$longi_dir/$tp/SNP
    fac_dir=$longi_dir/$tp/facet

    echo "============================================"
    echo "Processing $SCENARIO / $samp_name / $tp"
    echo "============================================"

    ###################### MANTA ######################
    conda activate "$ENV_MANTA"

    rm -rf $manta_dir
    echo "Running Manta ($tp)..."

    "$ENV_MANTA/bin/python" "$MANTA_CONFIG" \
        --normalBam $nbam \
        --tumorBam $tbam \
        --referenceFasta $ref \
        --runDir $manta_dir

    "$ENV_MANTA/bin/python" $manta_dir/runWorkflow.py -j $thread

    "$ENV_MANTA/bin/python" "$MANTA_CONVERT" \
        $samtool_path $ref \
        $manta_dir/results/variants/somaticSV.vcf.gz \
        > $manta_dir/results/variants/${samp_name}_${tp}.vcf

    conda deactivate

    ##################### SVTYPER #####################
    conda activate "$ENV_SVTYPER"

    echo "Running SVtyper ($tp)..."
    mkdir -p $svtyp_dir
    rm -f $svtyp_dir/*

    # Add CIPOS/CIEND tags for SVtyper compatibility
    awk 'BEGIN{FS=OFS="\t"}
         /^#/ { print; next }
         {
           info = $8
           gsub(/CIPOS=[^;]+;?/,"",info)
           gsub(/CIEND=[^;]+;?/,"",info)
           $8 = "CIPOS=-100,100;CIEND=-100,100;" info
           print
         }' $manta_dir/results/variants/${samp_name}_${tp}.vcf \
         > $svtyp_dir/${samp_name}_${tp}_w_CI.vcf

    svtyper \
        -i $svtyp_dir/${samp_name}_${tp}_w_CI.vcf \
        -B $tbam \
        -l $svtyp_dir/${samp_name}_${tp}.bam.json \
        > $svtyp_dir/svt_${samp_name}_${tp}.vcf

    conda deactivate

    ####################### SNP #######################
    conda activate "$ENV_GATK"

    echo "Running SNP calling ($tp)..."
    mkdir -p $snp_dir
    rm -f $snp_dir/*

    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R $ref \
        -I $tbam \
        -O $snp_dir/SNP.vcf.gz

    bcftools view -v snps -g het -Oz -o $snp_dir/het_snp.vcf.gz $snp_dir/SNP.vcf.gz
    tabix -p vcf $snp_dir/het_snp.vcf.gz

    conda deactivate

    ###################### PHASE ######################
    conda activate "$ENV_VISOR"

    echo "Running phasing ($tp)..."

    # Get SV genomic ranges
    Rscript $scrp_dir/get_sv_range.R \
        -s ${samp_name}_${tp} \
        -o $snp_dir \
        -v $svtyp_dir/svt_${samp_name}_${tp}.vcf

    # Extract het SNPs near SVs
    bcftools view \
        -R $snp_dir/${samp_name}_${tp}.bed \
        $snp_dir/het_snp.vcf.gz \
        -Ov -o $snp_dir/het_near_sv_${samp_name}_${tp}.vcf

    # Get het SNPs on SV-supporting reads
    awk 'BEGIN { OFS="\t" } !/^#/ { print $1, $2-1, $2 }' \
        $snp_dir/het_near_sv_${samp_name}_${tp}.vcf \
        > $snp_dir/pos_${samp_name}_${tp}.bed

    samtools view -f 1 -F 2 -b $tbam > $snp_dir/sup_${samp_name}_${tp}.bam
    samtools index $snp_dir/sup_${samp_name}_${tp}.bam

    bcftools mpileup -f $ref -a DP,AD -A \
        -R $snp_dir/pos_${samp_name}_${tp}.bed \
        $snp_dir/sup_${samp_name}_${tp}.bam \
        -Ov > $snp_dir/het_on_sv_${samp_name}_${tp}.vcf

    rm $snp_dir/sup_${samp_name}_${tp}.bam*

    conda deactivate

    ###################### FACETS ######################
    conda activate "$ENV_FACET"

    echo "Running FACETS ($tp)..."
    mkdir -p $fac_dir
    rm -f $fac_dir/*

    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

    "$FACET_SNP_PILEUP" -g \
        -q15 -Q20 -P100 -r25,0 \
        $snp_dir/SNP.vcf.gz \
        $fac_dir/tmp_${samp_name}_${tp}.csv \
        $nbam \
        $tbam

    zcat $fac_dir/tmp_${samp_name}_${tp}.csv.gz \
        | sed '1!s/^chr//' \
        | gzip -c > $fac_dir/${samp_name}_${tp}.csv.gz

    Rscript "$FACET_R_SCRIPT" \
        -s ${samp_name}_${tp} \
        -p $fac_dir

    conda deactivate

    ##################### CLEANUP #####################
    echo "Cleaning up intermediates ($tp)..."

    #rm -f $svtyp_dir/${samp_name}_${tp}_w_CI.vcf
    #rm -f $svtyp_dir/${samp_name}_${tp}.bam.json

    #rm -f $snp_dir/SNP.vcf.gz $snp_dir/SNP.vcf.gz.tbi
    #rm -f $snp_dir/het_snp.vcf.gz $snp_dir/het_snp.vcf.gz.tbi
    #rm -f $snp_dir/${samp_name}_${tp}.bed
    #rm -f $snp_dir/het_near_sv_${samp_name}_${tp}.vcf
    #rm -f $snp_dir/pos_${samp_name}_${tp}.bed

    #rm -f $fac_dir/tmp_${samp_name}_${tp}.csv.gz

    # Manta: keep results/ and config files (runWorkflow.py*), delete only workspace
    #rm -rf $manta_dir/workspace

    echo "Done: $tp"

done

echo "=== Pipeline complete for $SCENARIO / $samp_name ==="
