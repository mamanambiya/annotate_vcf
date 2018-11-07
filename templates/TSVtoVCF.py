#!/usr/bin/env python2.7

import argparse,sys
import time as t

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", help="")
parser.add_argument("--outVCF", help="")
args = parser.parse_args()

def toVCF(TSV, VCF):
    '''
    :param TSV:
    :param VCF:
    :return:
    '''
    out = open(VCF, 'w')
    for line in open(TSV):
        if line.startswith('##'):
            out.writelines(line)
        else:
            if 'CHROM' in line:
                head = line.strip().split(' ')
                out.writelines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            else:
                line = line.strip().split('\t')
                if len(line) > 1:
                    for i in range(4,len(line)):
                        line[i] = head[i]+'='+line[i]
                    out.writelines('\t'.join([line[0],line[1],'.',line[2],line[3],'.','.',';'.join(line[4:])])+'\n')
    out.close()
toVCF(args.inTSV, args.outVCF)