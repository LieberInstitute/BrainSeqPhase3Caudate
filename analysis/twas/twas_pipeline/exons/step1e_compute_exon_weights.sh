#!/bin/bash
#$ -cwd -V
#$ -t 300001-354440:1
#$ -l mem_free=22.0G,h_vmem=25G,h_fsize=50G
#$ -N cw_exon_part5
#$ -e log_files/cw_exon_error.out
#$ -o log_files/cw_exon_summary.out

SIZE=1

bash ../_h/parallel_compute_weights.sh --batch_number $SGE_TASK_ID \
     --batch_size $SIZE --threads 1 --feature exons
