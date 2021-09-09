#!/bin/bash

FEATURE="genes"

python /ceph/opt/enhanced_fastqtl/python/run_FastQTL_threaded.py \
       --covariates ../../../covariates/_m/${FEATURE}.combined_covariates.txt \
       --maf_threshold 0.01 --chunks 250 --threads 32 \
       --window 0.5e6 --ma_sample_threshold 10 \
       ../../../vcf/_m/genotypes_chr.vcf.gz \
       ../../_m/${FEATURE}.expression.bed.gz \
       Brainseq_LIBD
