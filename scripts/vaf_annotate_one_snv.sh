#!/bin/bash
module load vcfanno
module load tabix/0.2.6

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       annotates one merged snv tumor vcf"
    exit 1
}

readonly VARIANT=snv_mnv
readonly INDIR=dbsnp_annotated/${VARIANT}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly input_file=${INDIR}/${ID}.annotated.${VARIANT}.vcf

if [ ! -f $input_file ] 
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

readonly OUTDIR=annotated/${VARIANT}
mkdir -p $OUTDIR
readonly output_file=${OUTDIR}/${variant}/${ID}.annotated.${VARIANT}.vcf
sed -e "s/@@SAMPLE@@/${ID}/" annotation/vaf.annotations.conf.template > annotation/vaf.${ID}.conf

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/vaf.${ID}.conf ${input_file} > ${output_file}
bgzip ${output_file}
tabix -p vcf ${output_file}.gz
rm annotation/vaf.${ID}.conf
