#!/bin/bash -l

module purge
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load python/2.7.2
module load gcc/4.8.1 openblas python-packages/2

readonly SAMPLE=$1
readonly VARIANT=$2

readonly USAGE="Usage: $0 sample-id snv_mnv|indel|indel_normed"

if [ -z ${SAMPLE} ] || [ -z ${VARIANT} ]
then
    >&2 echo $USAGE
    >&2 echo "  combines the merged VCF file with annotations from all callers"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

if [ ${VARIANT} != "snv_mnv" ] && [ ${VARIANT} != "indel" ] && [ ${VARIANT} != "indel_normed" ]
then
    >&2 echo $USAGE
    >&2 echo " Invalid VARIANT type ${VARIANT}"
    >&2 echo "invocation: $0 $1 $2"
    exit
fi

readonly merged=./annotated/${VARIANT}/${SAMPLE}.annotated.${VARIANT}.vcf.gz

if [ ! -f $merged ]
then
    >&2 echo $USAGE
    >&2 echo "Invalid sample ${SAMPLE}"
    exit
fi

readonly ANNOTATED_DIR=./input_vcfs/oxog_v1/${VARIANT}
readonly broad=${ANNOTATED_DIR}/${SAMPLE}_annotated_broad_SNV.vcf.gz
readonly dkfz=${ANNOTATED_DIR}/${SAMPLE}_annotated_dkfz_embl_SNV.vcf.gz
readonly muse=${ANNOTATED_DIR}/${SAMPLE}_annotated_muse_SNV.vcf.gz

readonly sanger=${ANNOTATED_DIR}/${SAMPLE}_annotated_sanger_SNV.vcf.gz

if [ ! -f $broad ] || [ ! -f $dkfz ] || [ ! -f $sanger ]
then
    >&2 echo $USAGE
    >&2 echo "  Input files not found: one of "
    >&2 echo "  $broad"
    >&2 echo "  $dkfz"
    >&2 echo "  $muse"
    >&2 echo "  $sanger"
    exit
fi

readonly OUTPUTDIR=./annotated/caller_vafs/${VARIANT}
mkdir -p ${OUTPUTDIR}
readonly outputfile=${OUTPUTDIR}/${SAMPLE}.annotated.${VARIANT}.vcf

module load python/2.7.2
module load gcc/4.8.1 openblas python-packages/2

function addheaders {
    sed -e '/#CHROM/i\
##INFO=<ID=VAFs,Number=.,Type=Float,Description="Variant Allele Fractions identified by callers, in order of callers in Callers record">\
##INFO=<ID=DPs,Number=.,Type=Integer,Description="Total tumour depths by callers, in order of callers in Callers record">\
##INFO=<ID=RCs,Number=.,Type=Integer,Description="Total alt-supporting read counts by callers, in order of callers in Callers record">\
##INFO=<ID=medianVAF,Number=1,Type=Float,Description="median VAF identified by callers">\
##INFO=<ID=weightedmeanVAF,Number=1,Type=Float,Description="sum of all seen alt read counts / sum of all depths">\
##INFO=<ID=closest_depth,Number=1,Type=Integer,Description="depth from caller giving VAF closest to wmVAF">\
##INFO=<ID=closest_readcount,Number=1,Type=Integer,Description="readcount from caller giving VAF closest to wmVAF">'
}

if [ ${VARIANT} == "snv_mnv" ]
then
    python ./scripts/annotate_caller_vafs.py $merged \
        -b ${broad} \
        -d ${dkfz} \
        -s <( zcat ${sanger} | sed '/=$/d' ) \
        -m ${muse} \
        | addheaders \
            > ${outputfile}
else
    python ./scripts/annotate_caller_vafs.py $merged \
        -b ${filenames[0]} \
        -d ${filenames[1]} \
        -s <( zcat ${filenames[2]} | sed '/=$/d' ) \
        -i \
        | addheaders \
            > ${outputfile}
fi
