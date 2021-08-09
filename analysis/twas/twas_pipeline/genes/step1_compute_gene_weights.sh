#!/bin/bash

SIZE=1

ln -sfn . ./output
mkdir log_files

export END=`grep -v Br ../../_m/genes/twas_gene_expression.txt | wc -l`

## Full run
for ii in `seq 100 $SIZE $END`; do
    qsub -V -l mem_free=26.0G,h_vmem=30G,h_fsize=50G -N cw_${ii} -o \
         log_files/cw_gene_summary.out -e log_files/cw_gene_error.out \
         -cwd -b y "bash ../_h/parallel_compute_weights.sh --batch_number $ii \
                         --batch_size $SIZE --threads 1 --feature genes"
done
