#!/bin/bash -l
#
# Calls an R script to apply a model given a model file,
# an input, an output, and optionally a threshold

module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load R/3.2.5
module load gcc/4.8.1 openblas python/2.7.2 python-packages/2
module load tabix/0.2.6

readonly DEFTHRESH=0.66

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ ! -f "$1" ] || [ ! -f "$2" ]
then
    >&2 echo "$0 - Apply an ensemble model (provided) to the VCF"
    >&2 echo "Usage: $0 modelfile input.vcf raw_output.vcf output.vcf [thresh=${DEFTHRESH}]"
    exit 
fi

echo "$0 $1 $2 $3 $4 $5"

readonly MODEL=$1
readonly INPUTVCF=$2
readonly RAWOUTPUTVCF=$3
readonly OUTPUTVCF=$4
readonly THRESH=${5:-${DEFTHRESH}}

nin=$( zcat ${INPUTVCF} | grep -vc "^#" )

# new header
cat <<EOF > ${OUTPUTVCF}
##fileformat=VCFv4.1
##INFO=<ID=NumCallers,Number=1,Type=Integer,Description="Number of callers that made this call">
##INFO=<ID=Callers,Number=.,Type=String,Description="Callers that made this call">
##INFO=<ID=1000genomes_AF,Number=1,Type=Float,Description="Thousand Genomes phase 3 occurance fraction if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">
##INFO=<ID=1000genomes_ID,Number=1,Type=String,Description="Thousand Genomes phase 3 ID if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">
##INFO=<ID=VAF,Number=1,Type=Float,Description="VAF from mutect read filter if available">
##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count from mutect read filter if available">
##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour alt count from mutect read filter if available">
##INFO=<ID=cosmic,Number=1,Type=String,Description="(first) cosmic ID if found, COSMICv76">
##INFO=<ID=dbsnp,Number=1,Type=String,Description="(first) dbSNP ID if found, build 147, All_20160408.vcf.gz">
##INFO=<ID=repeat_masker,Number=1,Type=String,Description="Repeat masker region if in one">
##FILTER=<ID=LOWSUPPORT,Description="Not enough support in consensus model">
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
EOF

Rscript --vanilla scripts/filter_calls_by_model.R $MODEL $INPUTVCF $RAWOUTPUTVCF $THRESH \
    && ./scripts/clean_indel_calls.py <( grep -v "##contig" $RAWOUTPUTVCF ) | grep -v "^#" >> $OUTPUTVCF \
    && bgzip -f ${OUTPUTVCF} \
    && tabix -p vcf ${OUTPUTVCF}.gz \
    && nin=$( zcat ${INPUTVCF} | grep -v "^#" | grep -vc "Callers=broad;" ) \
    && nout=$( zcat ${OUTPUTVCF}.gz | grep -vc "^#" ) \
    && if [ $nin -ne $nout ]; then echo "TRUNCATED: ${INPUTVCF} ${OUTPUTVCF}"; fi
