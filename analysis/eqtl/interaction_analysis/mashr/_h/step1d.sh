#!/bin/bash

FEATURE="junctions"

Rscript ../_h/parallel_mashr_junctions.R --feature $FEATURE --run_chunk \
        --chunk_size 500 --threads 2 1>$FEATURE/output.log 2>$FEATURE/error.log
