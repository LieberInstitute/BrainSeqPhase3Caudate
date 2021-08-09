#!/bin/bash

FUSION="../_h"
WEIGHTS="../../_m"
PGC2="../../../../../inputs/sz_gwas/pgc2_clozuk/summary_statistics/_m"
LDREF="../../../../../inputs/genotypes/ld_references/reference_hg38/_m/LDREF_hg38"

python ../_h/get_glist.py

mkdir sig_analysis

for chr in `seq 1 22`;do
    Rscript $FUSION/FUSION.assoc_test.R \
            --weights $WEIGHTS/Phase3_Caudate.pos \
            --sumstats $PGC2/clozuk_pgc2.meta.reformatted.sumstats_hg38 \
            --ref_ld_chr $LDREF/1000G.EUR. \
            --out PGC2.SCZ.Caudate.${chr}.dat \
            --weights_dir $WEIGHTS \
            --chr $chr --perm 100

    Rscript $FUSION/FUSION.post_process.R \
            --sumstats $PGC2/clozuk_pgc2.meta.reformatted.sumstats_hg38 \
            --input PGC2.SCZ.Caudate.${chr}.dat \
            --out PGC2.SCZ.Caudate.${chr}.analysis \
            --ref_ld_chr $LDREF/1000G.EUR. \
            --chr ${chr} --verbose 2 --plot --plot_eqtl \
            --locus_win 100000 --report

    mv -v *.analysis* sig_analysis/
done

Rscript $FUSION/FUSION.post_process.R \
        --sumstats $PGC2/clozuk_pgc2.meta.reformatted.sumstats_hg38 \
        --input PGC2.SCZ.Caudate.6.dat.MHC \
        --out PGC2.SCZ.Caudate.6.analysis.MHC \
        --ref_ld_chr $LDREF/1000G.EUR. \
        --chr 6 --verbose 2 --plot --plot_eqtl \
        --locus_win 100000 --report

mv -v *.analysis* sig_analysis/
