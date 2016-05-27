#!/bin/bash
module load tabix/0.2.6
module load parallel

readonly INDIR=annotated/snv_mnv
readonly OUTDIR=filtered/snv_mnv

mkdir -p ${OUTDIR}

function copy_and_reheader {
    local file=$1
    local newdir=filtered/snv_mnv

    base=$( basename $file )
    nogz=$( basename $base .gz )

    echo "${nogz}"
    cat > ${newdir}/${nogz} <<EOF
##fileformat=VCFv4.1
##INFO=<ID=NumCallers,Number=1,Type=Integer,Description="Number of callers that made this call">
##INFO=<ID=Callers,Number=.,Type=String,Description="Callers that made this call">
##INFO=<ID=1000genomes_AF,Number=A,Type=Float,Description="Thousand Genomes phase 3 occurance fraction if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">
##INFO=<ID=1000genomes_ID,Number=1,Type=String,Description="Thousand Genomes phase 3 ID if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">
##INFO=<ID=VAF,Number=1,Type=Float,Description="VAF from mutect read filter if available">
##INFO=<ID=t_alt_count,Number=1,Type=Integer,Description="Tumour alt count from mutect read filter if available">
##INFO=<ID=t_ref_count,Number=1,Type=Integer,Description="Tumour alt count from mutect read filter if available">
##INFO=<ID=cosmic,Number=1,Type=String,Description="(first) cosmic ID if found, COSMICv76">
##INFO=<ID=dbsnp,Number=1,Type=String,Description="(first) dbSNP ID if found, build 147, All_20160408.vcf.gz">
##INFO=<ID=repeat_masker,Number=1,Type=String,Description="Repeat masker region if in one">
##FILTER=<ID=LOWSUPPORT,Description="Not called by enough callers in ensemble">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
    zgrep -v "^#" $file >> ${newdir}/${nogz}
    bgzip ${newdir}/${nogz}
    tabix -p vcf ${newdir}/${nogz}.gz
}

export -f copy_and_reheader

find ${INDIR} -type f -name "*.gz" \
    | grep -vf blacklist.txt \
    | parallel -j 2 copy_and_reheader {}
