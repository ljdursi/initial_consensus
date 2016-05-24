#!/bin/bash

mkdir -p miscalls

for mergedfile in annotated/snv_mnv/*gz
do
    base=$( basename $mergedfile )
    sample=$( echo $base | cut -f 1 -d . )

    vaffile=vafs/${sample}.mutect_counts.vcf.gz

    if [ ! -f ${vaffile} ]
    then
        continue
    fi

    outfile=miscalls/${sample}.miscalls.txt

    zgrep -v "^#" ${mergedfile} \
        | grep -v "^#" \
        | grep -v VAF \
        | awk '{printf "%s_%s	%s\n", $1, $2, $0}' \
        | sort \
        > ${outfile}.tmp

    join <( zgrep -v "^#" ${vaffile} | awk '{printf "%s_%s	%s\n", $1, $2, $0}' | sort ) ${outfile}.tmp \
        | cut -f 2- -d' ' \
        > ${outfile}

    rm ${outfile}.tmp
done
