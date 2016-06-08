#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import sys

def round_dig(value, ndigits=2):
    scale=10**ndigits
    return int(num*scale)*1./scale

def main():
    parser = argparse.ArgumentParser(description='Annotate merged vcf with VAF information where available')
    parser.add_argument('inputvcf', type=argparse.FileType('r'), default=sys.stdin, help="Merged VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    args = parser.parse_args()

    reader = vcf.Reader(args.inputvcf)
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        new_info = {}
        # skip broad private calls
        if record.INFO['Callers'] == ['broad']:
            continue

        # skip calls that are defined to be SVs
        varlen = abs(len(record.REF) - len(record.ALT[0]))
        if varlen >= 100:
            continue

        # copy some records over directly
        for item in ['dbsnp', 'cosmic', 'Callers', 'NumCallers', 'repeat_masker', '1000genomes_AF', '1000genomes_ID']:
            if item in record.INFO:
                new_info[item] = record.INFO[item]

        if 'dbsnp_VP' in record.INFO:
            qualbyte = int(record.INFO['dbsnp_VP'],16) & 255
            somatic = (qualbyte & 2**5) > 0
            if somatic:
                new_info['dbsnp_somatic'] = True

        if ('TumorVarDepth' in record.INFO) and (record.INFO['TumorVarDepth'] > 0):
            new_info['t_vaf'] = record.INFO['TumorVAF']
            new_info['t_alt_count'] = record.INFO['TumorVarDepth']
            new_info['t_ref_count'] = record.INFO['TumorTotalDepth']-record.INFO['TumorVarDepth']

        record.INFO = new_info
        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
