#!/bin/bash
#$ -cwd -V
#$ -t 100000-153770:1
#$ -l mem_free=22.0G,h_vmem=25G,h_fsize=50G
#$ -N cw_jxn_part3
#$ -e log_files/cw_error.out
#$ -o log_files/cw_summary.out

SIZE=1

bash ../_h/parallel_compute_weights.sh --batch_number $SGE_TASK_ID \
     --batch_size $SIZE --threads 1 --feature junctions
