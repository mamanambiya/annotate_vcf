#!/home/mamana/miniconda3/envs/ngs_py27/bin/python

import argparse,sys,time,random

parser = argparse.ArgumentParser()
parser.add_argument("--sampleIN", help="")
parser.add_argument("--sampleOUT", help="")
parser.add_argument("--dataset", help="")
parser.add_argument("--sample_size", help="")
args = parser.parse_args()

def reduce_sample_size(sampleIN, sampleOUT, dataset='', sample_size=25):
    '''
    '''
    out = open(sampleOUT, 'w')
    samples = open(sampleIN).readlines()
    random.shuffle(samples)
    for sample in samples[:int(sample_size)]:
        out.writelines(sample)
    out.close()

reduce_sample_size(args.sampleIN, args.sampleOUT, args.dataset, args.sample_size)
