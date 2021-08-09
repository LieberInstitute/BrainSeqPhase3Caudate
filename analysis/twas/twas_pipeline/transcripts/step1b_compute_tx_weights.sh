#!/bin/bash
#$ -cwd -V
#$ -t 50001-102191:1
#$ -l mem_free=32.0G,h_vmem=35G,h_fsize=50G
#$ -N cw_tx_part2
#$ -e log_files/cw_tx_error.out
#$ -o log_files/cw_tx_summary.out

SIZE=1

bash ../_h/parallel_compute_weights.sh --batch_number $SGE_TASK_ID \
     --batch_size $SIZE --threads 1 --feature transcripts
