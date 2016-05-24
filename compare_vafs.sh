#!/bin/bash

echo "broad, ncallers, oxogVAF, medianVAF, wmeanVAF"
ls annotated/caller_vafs/snv_mnv/ \
    | grep -vf /.mounts/labs/simpsonlab/users/jdursi/interim_consensus_v2/conest/bad-samples.txt \
    | xargs -n1 -I{} cat annotated/caller_vafs/snv_mnv/{} \
    | grep -v "^#" \
    | grep -v "LOWSUPPORT" \
    | grep VAF \
    | cut -f 8 \
    | awk -F\; 'BEGIN {cols["NumCallers"]=1; cols["VAF"]=2; cols["medianVAF"]=3; cols["weightedmeanVAF"]=4; } \
    { for(i=1; i<=NF; i++) { split($i, kv, "="); if (kv[1] == "Callers") {printf "%d,",index(kv[2], "broad");} if (kv[1] in cols) {printf "%s,", kv[2]; } } print "";}' \
    | sed -e 's/,$//' 
