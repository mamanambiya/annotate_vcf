#!/usr/bin/env bash

HOMEDIR="/home/mamana/GAPW/main"
OUTDIR="/spaces/mamana/GAPW"

### Nextflowscript here
mdkir -p ${OUTDIR}/LOG
cd ${OUTDIR}
nextflow -log ${OUTDIR}/LOG/annotate_vcf.log \
    run mypandos \
    -w ${OUTDIR}/work \
    -resume \
    -profile pbs
