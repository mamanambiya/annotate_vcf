#!/usr/bin/env python2.7
'''

'''
import argparse,sys

parser = argparse.ArgumentParser()
parser.add_argument("--mappedBed", default="${mappedBed}", help="BED mapped")
parser.add_argument("--gwascat_b38", default="${gwascat_b38}", help="GWAS Catalog file b38")
parser.add_argument("--gwascat_b37", default="${gwascat_b37}", help="")
args = parser.parse_args()

def gwascat38to37(gwascat, mappedBed, outCat):
    '''
    :param mappedBed:
    :param outCat:
    :return:
    '''
    datas = {}
    for line in open(mappedBed):
        line = line.strip().split('->')
        try:
            map1 = line[0].strip().split()[0]+'-'+line[0].strip().split()[1]
            map2 = line[1].strip().split()[0]+'-'+line[1].strip().split()[1]
        except:
            print "Failed map", '\\t'.join(line)
        datas[map1] = map2
    out = open(outCat, 'w')
    for line in open(gwascat):
        line = line.strip().split('\\t')
        chr = str(line[11])
        pos = str(line[12])
        if 'x' in chr:
            chr = chr.split('x')[0]
            pos = pos.split('x')[0]
        if ';' in chr: chr = chr.split(';')[0]
        if ';' in pos: pos = pos.split(';')[0]
        if chr+'-'+pos in datas:
            chr1 = datas[chr+'-'+pos].split('-')[0]
            pos1 = datas[chr+'-'+pos].split('-')[1]
            line[11] = chr1
            line[12] = pos1
        out.writelines('\\t'.join(line)+'\\n')
    out.close()
    return

gwascat38to37(args.gwascat_b38, args.mappedBed, args.gwascat_b37)