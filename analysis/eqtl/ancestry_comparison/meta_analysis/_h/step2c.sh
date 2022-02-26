#!/bin/bash

FEATURE="exons"

# # Node 5
# Rscript ../_h/parallel_mashr.R --feature $FEATURE --run_chunk \
#         --chunk_size 1000 --threads 8 1>$FEATURE/output.log 2>$FEATURE/error.log

# Node 6 or 7
Rscript ../_h/parallel_mashr.R --feature $FEATURE --run_chunk \
        --chunk_size 1000 --threads 4 1>$FEATURE/output.log 2>$FEATURE/error.log
