#!/usr/bin/env bash

HOMEDIR="/home/mamana/GAPW/main"
OUTDIR="/spaces/mamana/GAPW"

nextflow pull mypandos/annotate_vcf
### Nextflowscript here
mkdir -p ${OUTDIR}/LOG
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/annotate_vcf.log \
    run mypandos/annotate_vcf/main.nf \
    -w ${OUTDIR}/work \
    -resume \
    -profile pbs
