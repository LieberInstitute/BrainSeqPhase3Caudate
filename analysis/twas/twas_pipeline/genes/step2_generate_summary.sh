#!/bin/bash

bash ../_h/generate_summary_stats.sh 2> Phase3_Caudate.err 1> \
     Phase3_Caudate.profile
python ../_h/gen_positions.py
