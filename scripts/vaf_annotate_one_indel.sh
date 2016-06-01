#!/bin/bash
module load vcfanno
module load tabix/0.2.6

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       annotates one merged snv tumor vcf"
    exit 1
}

readonly VARIANT=indel
readonly INDIR=dbsnp_annotated/${VARIANT}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly input_file=${INDIR}/${ID}.annotated.${VARIANT}.vcf
readonly nbroad=$( grep -c broad ${input_file} )

if [ ! -f $input_file ] 
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

readonly OUTDIR=annotated/${VARIANT}
mkdir -p $OUTDIR
readonly output_file=${OUTDIR}/${variant}/${ID}.annotated.${VARIANT}.vcf
if [[ $nbroad -eq 0 ]]
then
    template_file=annotation/vaf.indel.annotations-nobroad.conf.template 
else
    template_file=annotation/vaf.indel.annotations.conf.template 
fi

sed -e "s/@@SAMPLE@@/${ID}/" $template_file > annotation/vaf.indel.${ID}.conf

if [[ -f $output_file ]] && [[ $outfile -nt $input_file ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/vaf.indel.${ID}.conf ${input_file} > ${output_file}
bgzip -f ${output_file}
tabix -p vcf ${output_file}.gz
rm annotation/vaf.indel.${ID}.conf
