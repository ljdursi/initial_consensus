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

readonly MODEL=$1
readonly INPUTVCF=$2
readonly RAWOUTPUTVCF=$3
readonly OUTPUTVCF=$4
readonly THRESH=${5:-${DEFTHRESH}}

nin=$( zcat ${INPUTVCF} | grep -vc "^#" )

Rscript --vanilla scripts/filter_calls_by_model.R $MODEL $INPUTVCF $RAWOUTPUTVCF $THRESH \
    && ./scripts/clean_indel_calls.py <( grep -v "##contig" $RAWOUTPUTVCF ) -o $OUTPUTVCF \
    && bgzip -f ${OUTPUTVCF} \
    && tabix -p vcf ${OUTPUTVCF}.gz \
    && nin=$( zcat ${INPUTVCF} | grep -v "^#" | grep -vc "Callers=broad;" ) \
    && nout=$( zcat ${OUTPUTVCF}.gz | grep -vc "^#" ) \
    && if [ $nin -ne $nout ]; then echo "TRUNCATED: ${INPUTVCF} ${OUTPUTVCF}"; fi
