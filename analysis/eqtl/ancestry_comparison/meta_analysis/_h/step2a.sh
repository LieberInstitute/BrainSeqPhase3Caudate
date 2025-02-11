#!/bin/bash

FEATURE="genes"

Rscript ../_h/parallel_mashr.R --feature $FEATURE --run_chunk \
        --chunk_size 250 --threads 8 1>$FEATURE/output.log 2>$FEATURE/error.log
