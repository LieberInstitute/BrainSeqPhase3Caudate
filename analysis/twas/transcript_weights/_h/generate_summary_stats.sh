#!/bin/bash

FUSION="/ceph/opt/fusion_twas-master"

# ls WEIGHTS/*RDat > Phase3_Caudate.list
find WEIGHTS/ -type f -name "*RDat" -exec ls {} + > Phase3_Caudate.list

Rscript $FUSION/utils/FUSION.profile_wgt.R Phase3_Caudate.list
