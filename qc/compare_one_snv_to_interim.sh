#!/bin/bash

readonly sample=$1

if [ -z $sample ]
then
    >&2 echo "Usage: $0 sample_id"
    exit 1
fi


readonly VARIANT=snv_mnv

readonly NEWDIR=/oicr/data/pancanxfer/consensus/initial_consensus/annotated/${VARIANT}
readonly OLDDIR=/.mounts/labs/simpsonlab/users/jdursi/interim_consensus/results/${VARIANT}
readonly MERGEDIR=./qc/compare-to-interim/${VARIANT}

mkdir -p ${MERGEDIR}

oldfile=${OLDDIR}/${sample}.merged.somatic.${VARIANT}.vcf
newfile=${NEWDIR}/${sample}.annotated.${VARIANT}.vcf.gz

if [ ! -f $newfile ] 
then
    echo "${sample}: new file not present: ${newfile}"
    exit 0
fi
if [ ! -f $oldfile ]
then
    echo "${sample}: old file not present ${oldfile}"
    exit 0
fi
echo -n ${sample}
newcalls=$( zgrep -cv "^#" $newfile )
oldcalls=$( grep -cv "^#" $oldfile )

if [ "$newcalls" -ne "$oldcalls" ] 
then
    echo ": mismatch number of calls (new: ${newcalls} old: ${oldcalls} )"
else
    echo ""
fi
if [ ${newcalls} -eq 0 ] || [ ${oldcalls} -eq 0 ] 
then
    exit 0
fi

outfile=${MERGEDIR}/${sample}.comparison.${VARIANT}.vcf
mergevcf -l interim,new_mutect_vafs \
    ${oldfile} ${newfile} \
    --ncallers \
    -o ${outfile}
