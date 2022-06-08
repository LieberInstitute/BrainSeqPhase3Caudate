#!/bin/bash

FEATURE="junctions"

# Node 5
Rscript ../_h/parallel_mashr_junctions.R --feature $FEATURE --run_chunk \
        --chunk_size 1500 --threads 4 1>$FEATURE/output.log 2>$FEATURE/error.log

# # Node 6 or 7
# Rscript ../_h/parallel_mashr_junctions.R --feature $FEATURE --run_chunk \
#         --chunk_size 1500 --threads 2 1>$FEATURE/output.log 2>$FEATURE/error.log
