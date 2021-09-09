#!/bin/bash

FEATURE="transcripts"

Rscript ../_h/parallel_mashr.R --feature $FEATURE --run_chunk \
        --chunk_size 2500 --threads 8 1>output.log 2>error.log
