# Annotate significant eqtls combining nominal and permutation eqtls
# See:
# https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/fastqtl.wdl

python /ceph/opt/enhanced_fastqtl/python/annotate_outputs.py \
    ../../fastqtl_permutation/_m/Brainseq_LIBD.genes.txt.gz \
    0.05 ../../gtf/_m/feature.gtf \
    --nominal_results ../../fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz
