#!/usr/bin/env python2.7
'''

'''
import argparse,sys,time,gzip

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", default="${inTSV}", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--outTSV", default="${outTSV}", help="TSV annotation file")
parser.add_argument("--chrm", default="${chrm}", help="chromosome to consider")
args = parser.parse_args()

def mergeMAFannotations(inTSV, chrm='', outTSV=''):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    outTSV_out = open(outTSV, 'w')
    print "Reading ",inTSV
    for line in open(inTSV):
        if "#CHRM" in line or '#chr' in line:
            line = line.strip().split('\\t')
            header = line[1:]
            outTSV_out.writelines('\\t'.join(['#CHRM:POS'] + header) + '\\n')
        else:
            line = line.strip().split('\\t')
            try:
                chr = line[0].split(':')[0]
                pos = line[0].split(':')[1]
            except:
                pass
            if str(chr) == str(chrm):
                chr_pos = chr + ':' + pos
                maf = line[1:]
                maf_ = []
                for i in range(len(header)):
                    try:
                        maf_i = float(maf[i])
                        if maf_i <= 0.5:
                            maf_.append(header[i] + '=' + str(maf_i))
                        else:
                            maf_.append(header[i] + '=' + str(1 - maf_i))
                    except:
                        print maf
                outTSV_out.writelines('\\t'.join([chr_pos] + maf_) + '\\n')
    outTSV_out.close()

if __name__ == '__main__':
    if args.inTSV != '' and args.chrm != '' and args.outTSV != '':
        mergeMAFannotations(args.inTSV, args.chrm, args.outTSV)