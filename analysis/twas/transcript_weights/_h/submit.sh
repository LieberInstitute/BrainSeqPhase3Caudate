#!/bin/bash

SIZE=1

# ln -s . ./output

export end=`grep -v Br ../../_m/transcripts/twas_gene_expression.txt | wc -l`

# for ii in `seq 1 $SIZE $end`; do
#     qsub -V -l mem_free=26.0G,h_vmem=30G,h_fsize=50G -N cw_${ii} -o cw_${ii}_summary.out -e cw_${ii}_error.out \
# 	 -cwd -b y "bash ../_h/parallel_compute_weights.sh --batch_number $ii --batch_size $SIZE --threads 1"
# done

# for ii in `cat ../_h/num_error.txt`; do
#     qsub -V -l mem_free=26.0G,h_vmem=30G,h_fsize=50G -N cw_${ii} -o cw_${ii}_summary.out -e cw_${ii}_error.out \
# 	 -cwd -b y "bash ../_h/parallel_compute_weights.sh --batch_number $ii --batch_size $SIZE --threads 1"
# done

# for ii in `cat ../_h/stalled.txt`; do
#     qsub -V -l mem_free=32.0G,h_vmem=35G,h_fsize=50G -N cw_${ii} -o cw_${ii}_summary.out -e cw_${ii}_error.out \
# 	 -cwd -b y "bash ../_h/parallel_compute_weights.sh --batch_number $ii --batch_size $SIZE --threads 1"
# done

for ii in 1; do
    qsub -V -l mem_free=32.0G,h_vmem=35G,h_fsize=50G -N cw_${ii} -o cw_${ii}_summary.out -e cw_${ii}_error.out \
	 -cwd -b y "bash ../_h/parallel_compute_weights.sh --batch_number $ii --batch_size $SIZE --threads 1"
done

#bash ../_h/generate_summary_stats.sh 2> Phase3_Caudate.err 1> Phase3_Caudate.profile
#python ../_h/gen_positions.py 
