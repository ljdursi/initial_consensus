#!/bin/bash -l

module purge
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

function usage {
    >&2 echo "usage: $0 sample-id "
    >&2 echo "       merges and annotates indel VCFs for a single tumour"
    exit 1
}

readonly ID=$1
if [ -z $ID ]
then
    >&2 echo "argument missing"
    usage
fi

./scripts/merge-one-tumour-indel.sh $ID
./scripts/dbsnp_annotate_one_indel.sh $ID
./scripts/vaf_annotate_one_indel.sh $ID
