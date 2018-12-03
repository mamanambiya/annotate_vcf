#!/usr/bin/env bash

HOMEDIR="/home/mamana/annotate_vcf"
OUTDIR="/spaces/mamana/annotate_vcf"

nextflow pull mypandos/annotate_vcf
### Nextflowscript here
mkdir -p ${OUTDIR}/LOG
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/annotate_vcf.log \
    -C ${HOMEDIR}/WITS/nextflow.config \
    run mypandos/annotate_vcf/main.nf \
    -w ${OUTDIR}/work \
    -resume \
    -profile slurm
