#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--cosmicVCF", help="")
parser.add_argument("--outVCF", help="")
args = parser.parse_args()

def fix_cosmic(inVCF, outVCF):
    '''
    :param inVCF:
    :param outVCF:
    :return:
    '''
    eff = ['GENE', 'STRAND', 'CDS', 'AA', 'CNT', 'SNP']
    out = open(outVCF, 'w')
    for line in open(inVCF):
        if not line.startswith("#"):
            line = line.strip().split('\t')
            ann = line[-1].split(';')
            if len(ann) != 6:
                ann.append('N/A')
                line[-1] = ';'.join(ann)
            line = '\t'.join(line) + '\n'
        out.writelines(line)
    out.close()

if(args.cosmicVCF):
    fix_cosmic(args.cosmicVCF, args.outVCF)