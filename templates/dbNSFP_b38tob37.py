#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
'''

'''
import argparse,sys,time


parser = argparse.ArgumentParser()
parser.add_argument("--inDB", help="")
parser.add_argument("--outDB", help="")
args = parser.parse_args()

def dbnsfp38to37(inDB, outDB):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    datas = {}
    out = open(outDB, 'w')
    for line in open(inDB):
        if 'hg19_pos' in line:
            out.writelines(line)
        else:
            line = line.strip().split('\t')
            pos38 = line[1]
            pos37 = line[8]
            if pos37 != '.' and pos37 != '':
                line[1] = pos37
                line[8] = pos38
                datas[pos37] = '\t'.join(line)+"\n"
    data = [str(pos) for pos in sorted([int(pos) for pos in datas])]
    for pos in data:
        out.writelines(datas[pos])
    out.close()

if(args.inDB):
    dbnsfp38to37(args.inDB, args.outDB)