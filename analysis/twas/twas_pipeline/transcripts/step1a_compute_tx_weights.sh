#!/bin/bash
#$ -cwd -V
#$ -t 1-50000
#$ -l mem_free=32.0G,h_vmem=35G,h_fsize=50G
#$ -N cw_tx_part1
#$ -e log_files/cw_tx_error.out
#$ -o log_files/cw_tx_summary.out

SIZE=1

ln -sfn . ./output

mkdir log_files

bash ../_h/parallel_compute_weights.sh --batch_number $SGE_TASK_ID \
     --batch_size $SIZE --threads 1 --feature transcripts
