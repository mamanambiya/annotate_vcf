#!/home/mamana/miniconda3/envs/ngs_py27/bin/python
'''

'''
import argparse,sys,time

parser = argparse.ArgumentParser()
parser.add_argument("--inTSV", help="One or many TSV annotation files, if many use ';' as separator")
parser.add_argument("--outTSV", help="TSV annotation file")
parser.add_argument("--inVCF", help="VCF file to annotate")
parser.add_argument("--chrm", help="chromosome to consider")
parser.add_argument("--outVCF", help="Annotated VCF file")
args = parser.parse_args()

annot = {'AGVP_AF':'AGVP', 'SAHGP_AF':'SAHGP', 'TRYPANOGEN_AF':'TRYPANOGEN', 'KG_AC':'KG', 'KG_AF':'KG', 'KG_AFR_AC':'KG', 'KG_AFR_AF':'KG', 'KG_EUR_AC':'KG', 'KG_EUR_AF':'KG', 'KG_AMR_AC':'KG', 'KG_AMR_AF':'KG', 'KG_EAS_AC':'KG', 'KG_EAS_AF':'KG', 'KG_SAS_AC':'KG', 'KG_SAS_AF':'KG', 'TWINSUK_AC':'TWINSUK', 'TWINSUK_AF':'TWINSUK', 'ALSPAC_AC':'ALSPAC', 'ALSPAC_AF':'ALSPAC', 'ESP6500_AA_AC':'ESP6500', 'ESP6500_AA_AF':'ESP6500', 'ESP6500_EA_AC':'ESP6500', 'ESP6500_EA_AF':'ESP6500', 'ExAC_AC':'ExAC', 'ExAC_AF':'ExAC', 'ExAC_Adj_AC':'ExAC', 'ExAC_Adj_AF':'ExAC', 'ExAC_AFR_AC':'ExAC', 'ExAC_AFR_AF':'ExAC', 'ExAC_AMR_AC':'ExAC', 'ExAC_AMR_AF':'ExAC', 'ExAC_EAS_AC':'ExAC', 'ExAC_EAS_AF':'ExAC', 'ExAC_FIN_AC':'ExAC', 'ExAC_FIN_AF':'ExAC', 'ExAC_NFE_AC':'ExAC', 'ExAC_NFE_AF':'ExAC', 'ExAC_SAS_AC':'ExAC', 'ExAC_SAS_AF':'ExAC', 'ExAC_nonTCGA_AC':'ExAC', 'ExAC_nonTCGA_AF':'ExAC', 'ExAC_nonTCGA_Adj_AC':'ExAC', 'ExAC_nonTCGA_Adj_AF':'ExAC', 'ExAC_nonTCGA_AFR_AC':'ExAC', 'ExAC_nonTCGA_AFR_AF':'ExAC', 'ExAC_nonTCGA_AMR_AC':'ExAC', 'ExAC_nonTCGA_AMR_AF':'ExAC', 'ExAC_nonTCGA_EAS_AC':'ExAC', 'ExAC_nonTCGA_EAS_AF':'ExAC', 'ExAC_nonTCGA_FIN_AC':'ExAC', 'ExAC_nonTCGA_FIN_AF':'ExAC', 'ExAC_nonTCGA_NFE_AC':'ExAC', 'ExAC_nonTCGA_NFE_AF':'ExAC', 'ExAC_nonTCGA_SAS_AC':'ExAC', 'ExAC_nonTCGA_SAS_AF':'ExAC', 'ExAC_nonpsych_AC':'ExAC', 'ExAC_nonpsych_AF':'ExAC', 'ExAC_nonpsych_Adj_AC':'ExAC', 'ExAC_nonpsych_Adj_AF':'ExAC', 'ExAC_nonpsych_AFR_AC':'ExAC', 'ExAC_nonpsych_AFR_AF':'ExAC', 'ExAC_nonpsych_AMR_AC':'ExAC', 'ExAC_nonpsych_AMR_AF':'ExAC', 'ExAC_nonpsych_EAS_AC':'ExAC', 'ExAC_nonpsych_EAS_AF':'ExAC', 'ExAC_nonpsych_FIN_AC':'ExAC', 'ExAC_nonpsych_FIN_AF':'ExAC', 'ExAC_nonpsych_NFE_AC':'ExAC', 'ExAC_nonpsych_NFE_AF':'ExAC', 'ExAC_nonpsych_SAS_AC':'ExAC', 'ExAC_nonpsych_SAS_AF':'ExAC', 'gnomAD_exomes_AC':'gnomAD', 'gnomAD_exomes_AN':'gnomAD', 'gnomAD_exomes_AF':'gnomAD', 'gnomAD_exomes_AFR_AC':'gnomAD', 'gnomAD_exomes_AFR_AN':'gnomAD', 'gnomAD_exomes_AFR_AF':'gnomAD', 'gnomAD_exomes_AMR_AC':'gnomAD', 'gnomAD_exomes_AMR_AN':'gnomAD', 'gnomAD_exomes_AMR_AF':'gnomAD', 'gnomAD_exomes_ASJ_AC':'gnomAD', 'gnomAD_exomes_ASJ_AN':'gnomAD', 'gnomAD_exomes_ASJ_AF':'gnomAD', 'gnomAD_exomes_EAS_AC':'gnomAD', 'gnomAD_exomes_EAS_AN':'gnomAD', 'gnomAD_exomes_EAS_AF':'gnomAD', 'gnomAD_exomes_FIN_AC':'gnomAD', 'gnomAD_exomes_FIN_AN':'gnomAD', 'gnomAD_exomes_FIN_AF':'gnomAD', 'gnomAD_exomes_NFE_AC':'gnomAD', 'gnomAD_exomes_NFE_AN':'gnomAD', 'gnomAD_exomes_NFE_AF':'gnomAD', 'gnomAD_exomes_SAS_AC':'gnomAD', 'gnomAD_exomes_SAS_AN':'gnomAD', 'gnomAD_exomes_SAS_AF':'gnomAD', 'gnomAD_exomes_OTH_AC':'gnomAD', 'gnomAD_exomes_OTH_AN':'gnomAD', 'gnomAD_exomes_OTH_AF':'gnomAD', 'gnomAD_AC':'gnomAD', 'gnomAD_AN':'gnomAD', 'gnomAD_AF':'gnomAD', 'gnomAD_AFR_AC':'gnomAD', 'gnomAD_AFR_AN':'gnomAD', 'gnomAD_AFR_AF':'gnomAD', 'gnomAD_AMR_AC':'gnomAD', 'gnomAD_AMR_AN':'gnomAD', 'gnomAD_AMR_AF':'gnomAD', 'gnomAD_ASJ_AC':'gnomAD', 'gnomAD_ASJ_AN':'gnomAD', 'gnomAD_ASJ_AF':'gnomAD', 'gnomAD_EAS_AC':'gnomAD', 'gnomAD_EAS_AN':'gnomAD', 'gnomAD_EAS_AF':'gnomAD', 'gnomAD_FIN_AC':'gnomAD', 'gnomAD_FIN_AN':'gnomAD', 'gnomAD_FIN_AF':'gnomAD', 'gnomAD_NFE_AC':'gnomAD', 'gnomAD_NFE_AN':'gnomAD', 'gnomAD_NFE_AF':'gnomAD', 'gnomAD_OTH_AC':'gnomAD', 'gnomAD_OTH_AN':'gnomAD', 'gnomAD_OTH_AF':'gnomAD'}

annots = {"KG":"Allele Frequency from Thousand Genomes Project (1KG)",
          "gnomAD":"Allele Frequency from the Genome Aggregation Database (gnomAD)",
          "ExAC":"Allele Frequency from the Exome Aggregation Consortium (ExAC)",
          "ESP6500":"Allele Frequency from the NHLBI GO Exome Sequencing Project (ESP)",
          "ALSPAC":"Allele Frequency from the Avon Longitudinal Study of Parents and Children (ALSPAC)",
          "TWINSUK":"Allele Frequency from the TwinsUK project",
          "AGVP": "Allele Frequency among the African Genome Variation Project (AGVP)",
          "SAHGP": "Allele Frequency from the Southern African Human Genome Project (SAHGP)",
          "TRYPANOGEN": "Allele Frequency from H3Africa Trypanogen Project"}

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
            line = line.strip().split('\t')
            header = line[1:]
            outTSV_out.writelines('\t'.join(['#CHRM:POS'] + header) + '\n')
        else:
            line = line.strip().split('\t')
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
                outTSV_out.writelines('\t'.join([chr_pos] + maf_) + '\n')
    outTSV_out.close()


def readTSV(inTSV_list):
    '''
    :param inBed:
    :param outBed:
    :return:
    '''
    datas_tsv = {}
    inTSV_list = inTSV_list.split(';')
    for inTSV in inTSV_list:
        for line in open(inTSV):
            line = line.strip().split('\t')
            chr_pos = line[0]
            maf = ';'.join(line[1:])
            # if chr_pos.startswith('#CHRM'):
            #     if chr_pos not in datas_tsv:
            #         datas_tsv[chr_pos] = maf.split(';')
            #     else:
            #         datas_tsv[chr_pos] += maf.split(';')
            # else:
            if chr_pos not in datas_tsv:
                datas_tsv[chr_pos] = maf
            else:
                datas_tsv[chr_pos] += ';' + maf
    return datas_tsv

def readVCF(inTSV, inVCF, outVCF):
    '''
    :param inVCF:
    :return:
    '''
    out = open(outVCF, 'w')
    ann = readTSV(inTSV)
    for line in open(inVCF):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            chr = line[0]
            pos = line[1]
            id = chr+':'+pos
            if id in ann:
                line[7] += ';'+ ann[id]
            out.writelines('\t'.join(line)+'\n')
        else:
            head = ann['#CHRM:POS'].split(';')
            if line.startswith('#CHROM'):
                for info in head:
                    out.writelines("##INFO=<ID="+info+",Number=A,Type=Float,Description=\""+annots[annot[info]]+"\">\n")
            out.writelines(line)

if args.inVCF:
    readVCF(args.inTSV, args.inVCF, args.outVCF)
elif args.inTSV and not args.inVCF:
    mergeMAFannotations(args.inTSV, args.chrm, args.outTSV)