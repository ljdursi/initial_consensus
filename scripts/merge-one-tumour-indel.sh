#!/bin/bash -l

export PATH=/u/jdursi/.local/bin:${PATH}
module load python/2.7.2
module load gcc/4.8.1 openblas python-packages/2
module load tabix/0.2.6

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       merges indel VCFs for sample $0"
    exit 1
}

readonly VARIANT=indel

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

readonly INDIR=input_vcfs/processed/${VARIANT}

readonly broadfile=${INDIR}/${ID}.broad.${VARIANT}.vcf.gz
readonly dkfzfile=${INDIR}/${ID}.dkfz.${VARIANT}.vcf.gz
readonly sangerfile=${INDIR}/${ID}.sanger.${VARIANT}.vcf.gz
readonly smufinfile=${INDIR}/${ID}.smufin.${VARIANT}.vcf.gz  

usetwo=0
if [ ! -f $smufinfile ] || [ ! -f $broadfile ]
then
    >&2 echo "files missing: one of ${smufinfile} ${broadfile}.  Merging just two core callers"
    usetwo=1
fi

if [ ! -f $dkfzfile ] || [ ! -f $sangerfile ]
then
    >&2 echo "manditory files missing: one of ${dkfzfile} ${sangerfile} not found"
    usage
fi

readonly OUTDIR=merged/${VARIANT}
mkdir -p ${OUTDIR}

newest=$dkfzfile
if [[ $sangerfile -nt $newest ]]; then newest=$sangerfile; fi
if [[ $usetwo == 0 ]] && [[ $broadfile -nt $newest ]]; then newest=$broadfile; fi
if [[ $usetwo == 0 ]] && [[ $smufinfile -nt $newest ]]; then newest=$smufinfile; fi

readonly outfile=${OUTDIR}/${ID}.merged.${VARIANT}.vcf

if [[ -f $outfile ]] && [[ $outfile -nt $newest ]]
then
    >&2 echo "$0: ${outfile} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

if [ $usetwo == 0 ]
then
    mergevcf -l broad,dkfz,sanger,smufin \
        ${broadfile} ${dkfzfile} ${sangerfile} ${smufinfile} \
        --ncallers \
        -o ${outfile}.tmp
else
    mergevcf -l dkfz,sanger \
        ${dkfzfile} ${sangerfile} \
        --ncallers \
        -o ${outfile}.tmp
fi

grep "^#" ${outfile}.tmp > ${outfile}
grep -v "^#" ${outfile}.tmp \
    | grep -v "Callers=smufin;" \
    | sort -k1,1d -k2,2n \
    >> ${outfile}
rm ${outfile}.tmp

rm -f ${outfile}.gz ${outfile}.gz.tbi
bgzip ${outfile}
tabix -p vcf ${outfile}.gz
