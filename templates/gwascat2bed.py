#!/usr/bin/env python2.7
'''

'''
import argparse,sys

parser = argparse.ArgumentParser()
parser.add_argument("--toBed", default="${toBed}", help="BED file to extract columns")
parser.add_argument("--outBed", default="${outBed}", help="BED output of extracted columns")
args = parser.parse_args()

def tobed(inBed, outBed):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    unmapped = 0
    mapped = 0
    out = open(args.outBed, 'w')
    for line in open(args.toBed):
        line = line.strip().split('\\t')
        chr = line[11]
        pos = line[12]
        if chr == '' or pos == '':
            unmapped += 1
        else:
            mapped += 1
        if 'x' in chr:
            chr = chr.split('x')[0]
            pos = pos.split('x')[0]
        if ';' in chr: chr = chr.split(';')[0]
        if ';' in pos: pos = pos.split(';')[0]
        out.writelines(str(chr)+" "+str(pos)+" "+str(pos)+"\\n")
    out.close()
    return mapped, unmapped

tobed(args.toBed, args.outBed)