#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
# Author: Jeffrey M Kidd, Modified my Mamana M
# 2 September 2011
# add-ANC-to-vcf.py
# adds ancestral annotation based on ensemble takes from genome data archive
# you'll have to do your own filtering based on qual etc.

import sys
import os
import genomedata
import math

from genomedata import Genome
from optparse import OptionParser


USAGE = """
add-ANC-to-vcf.py --in <vcf file to process> --out <new VCF file name> -g <in/out is gziped>
                          --genomedata <path to genome data archieve with GERP scores>

Adds ancestral state SNPs in VCF, based on values in genomedata archieve (in 'anc' track).

Use -g if input VCF is gzipped, output file will also be gzipped.

Note: current version assumes all variants in VCF are SNPs.

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inVCF', help = 'input VCF file')
parser.add_option('--out',dest='outVCF', help = 'output VCF file')
parser.add_option('-g',dest='isGzip',  action='store_true', default = False, help = 'output VCF file')
parser.add_option('--genomedata',dest='genomedata', help = 'genomedata archive with GERP scores')


(options, args) = parser.parse_args()

if options.inVCF is None:
    parser.error('input VCF not given')
if options.outVCF is None:
    parser.error('output VCF not given')
if options.genomedata is None:
    parser.error('genomedata archive not given')

###############################################################################


# try to open up the genome data archieve
try:
    genome = Genome(options.genomedata)
except:
    print "ERROR!! Couldn't open the genomedata archive:  " + options.genomedata + "\n"
    sys.exit(1)

#setup file open/close with or without gzip
# if options.isGzip is True:
#     try:
#         gc = 'gunzip -c ' + options.inVCF
#         inFile = os.popen(gc, 'r')
#     except:
#         print "ERROR!! Couldn't open the file" + options.inVCF + " (with gzip)\n"
#         sys.exit(1)
#
#     try:
#         gc = 'gzip > ' + options.outVCF
#         outFile = os.popen(gc, 'w')
#     except:
#         print "ERROR!! Couldn't open the output file" + options.outVCF + " (with gzip)\n"
#         sys.exit(1)
# else:
inFile = open(options.inVCF, 'r')
outFile = open(options.outVCF, 'w')

# Set current chrom as something that isn't a chrom
currentChrom = 'notAChrom'
INFO = False
for line in inFile:
    # read through VCF file up to the chrom line, we will then add addtional info fields
    if line.startswith("#"):
        if line.startswith("##INFO") and INFO == False:
            # at this point, line is the 'header' line of the VCF.  Output header for the GERP info line
            outFile.writelines('##INFO=<ID=AA,Number=1,Type=Character,Description="ancestral state from ensemble">\n')
            INFO = True
        outFile.writelines(line)
    else:
        # rest of the VCF file should now just be the variants
        if line == "":
            break
        line = line.decode('utf-8').strip().split('\t')
        if line[0] != currentChrom:
            chrom = genome["chr"+line[0]]
            currentChrom = "chr"+line[0]
        pos = int(line[1]) - 1 #switch to zero based indexing
        score = chrom[pos,'anc']
        # check to see if there is a GERP score for the position, if not output line and continue
        # We should probably check ot see if the variant is not a SNP, as GERP isn't well defined
        # for non-SNP variants
        if math.isnan(score):
            anc = '.'
        else:
            anc = chr(score)

        if line[7] == '.':
            line[7] = 'AA=%s' % (anc.upper())
        else :
            rsField = ';AA=%s' % (anc.upper())
            line[7] += rsField
        outFile.writelines('\t'.join(line)+'\n')

genome.close()
inFile.close()
outFile.close()
