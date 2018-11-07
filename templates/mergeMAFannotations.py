#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--outTSV", help="TSV annotation file")
parser.add_argument("--chrm", help="chromosome to consider")
args = parser.parse_args()


def readTSV(inTSV_list, chrm='', outTSV=''):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    print chrm
    inTSV_list = inTSV_list.split(';')
    if outTSV != '':
        outTSV_out = open(outTSV,'w')
    datas_tsv = {}
    datas_tsv['CHRM:POS'] = []
    for tsv in inTSV_list:
        for line in open(tsv):
            if "#CHRM" in line or '#chr' in line:
                line = line.strip().split('\t')
                datas_tsv['CHRM:POS'] += line[1:]
                header = line[1:]
            else:
                line = line.strip().split('\t')

                try:
                    chr = line[0].split(':')[0]
                    pos = line[0].split(':')[1]
                except:
                    pos = line[0]
                    chr = line[0]
                maf = line[1:]
                maf_ = []
                if chrm == '':
                    for i in range(len(header)):
                        try:
                            maf_.append(header[i] + '=' + maf[i])
                        except:
                            print datas_tsv['CHRM:POS'], maf
                    if not chr + ':' + pos in datas_tsv:
                        datas_tsv[chr + ':' + pos] = []
                    datas_tsv[chr + ':' + pos] += maf_
                else:
                    if str(chrm) == str(chr):
                        for i in range(len(header)):
                            try:
                                maf_.append(header[i]+'='+maf[i])
                            except:
                                print datas_tsv['CHRM:POS'], maf
                        if not chr+':'+pos in datas_tsv:
                            datas_tsv[chr+':'+pos] = []
                        datas_tsv[chr+':'+pos] += maf_

    if outTSV != '':
        # outTSV_out.writelines('\t'.join(['#CHRM:POS'] + datas_tsv['CHRM:POS']) + '\n')
        for data in datas_tsv:
            if 'CHRM:POS' not in data:
                outTSV_out.writelines('\t'.join([data]+datas_tsv[data])+'\n')
        outTSV_out.close()

elif args.inTSV:
    if not args.chrm:
        args.chrm = ''
    readTSV(args.inTSV, args.chrm, args.outTSV)