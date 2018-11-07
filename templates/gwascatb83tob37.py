#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
'''

'''
import argparse,sys

parser = argparse.ArgumentParser()
parser.add_argument("--toBed", help="BED file to extract columns")
parser.add_argument("--outBed", help="BED output of extracted columns")
parser.add_argument("--mappedBed", help="BED mapped")
parser.add_argument("--gwascat", help="GWAS Catalog file")
parser.add_argument("--outb37", help="")
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
        line = line.strip().split('\t')
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
        out.writelines(str(chr)+" "+str(pos)+" "+str(pos)+"\n")
    out.close()
    return mapped, unmapped

def gwascat38to37(gwascat, mappedBed, outCat):
    '''
    :param mappedBed:
    :param outCat:
    :return:
    '''
    datas = {}
    for line in open(args.mappedBed):
        line = line.strip().split('->')
        try:
            map1 = line[0].strip().split()[0]+'-'+line[0].strip().split()[1]
            map2 = line[1].strip().split()[0]+'-'+line[1].strip().split()[1]
        except:
            print "Failed map", '\t'.join(line)
        datas[map1] = map2
    out = open(args.outb37, 'w')
    for line in open(args.gwascat):
        line = line.strip().split('\t')
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
        out.writelines('\t'.join(line)+'\n')
    out.close()
    return

if(args.toBed):
    tobed(args.toBed, args.outBed)
if(args.mappedBed):
    gwascat38to37(args.gwascat, args.mappedBed, args.outb37)