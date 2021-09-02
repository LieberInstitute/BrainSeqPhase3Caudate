#!/bin/bash

FEATURE="genes"

python /ceph/opt/enhanced_fastqtl/python/run_FastQTL_threaded.py \
       --covariates ../../../covariates/_m/${FEATURE}.combined_covariates.txt \
       --permute 1000 10000 --ma_sample_threshold 10 --threads 32 \
       --maf_threshold 0.01 --chunks 250 --window 0.5e6 \
       ../../../vcf/_m/genotypes_chr.vcf.gz \
       ../../_m/${FEATURE}.expression.bed.gz \
       Brainseq_LIBD
