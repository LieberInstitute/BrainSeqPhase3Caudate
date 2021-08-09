#!/bin/bash
#@Author: Kynon J. Benjamin
#@Edited from example script ($FUSION/examples/GEUV.compute_weights.sh)
show_help() {
    cat <<EOF
    Usage: ${0##*/} [-h] [-n BATCH NUMBER] [-s BATCH SIZE] [-t THREADS]
    Computer weights for TWAS analysis in batch of user defined size (default 100).

    -h|--help            display this help text
    -n|--batch_number N  initial batch number (default 0)
    -s|--batch_size N    size of batch run (default 100)
    -t|--threads N       number of threads (default 1)
    -f|--feature         feature to work on (default genes)
EOF
}

set -o errexit # stop execution on error

PLINK="/users/jbenjami/software/bin/plink"
FUSION="/dcl02/lieber/apaquola/opt/fusion_twas-master"
GCTA="/dcl02/lieber/apaquola/opt/fusion_twas-master/gcta64"
PRE_GENO="../../../../inputs/genotypes/byChrom/_m/LIBD_Brain_TopMed"
RSCRIPT="/dcl01/lieber/apaquola/trash/conda/env/my_svnR-3.5/R/3.5/bin/Rscript"
LDREF="../../../../inputs/genotypes/ld_references/reference_hg38/_m/LDREF_hg38"
GEMMA="/dcs04/lieber/ds2b/users/kynon/v4_phase3_paper/analysis/twas/gene_weights/_h/gemma-0.98.4"
CHR_INFO="/dcl01/lieber/apaquola/genomes/human/gencode26/GRCh38.p10.ALL/chromInfo/_m/chromInfo.txt"

# Reset all variables that might be set
NUM=0
SIZE=100
THREADS=1
FEATURE="genes"

# Initial output
while :; do
    case $1 in
        -h|=\?|--help)
            show_help
            exit
            ;;
        -n|--batch_number)
            if [ -n "$2" ]; then
                NUM=$2
                shift
            else
                printf 'ERROR: "--batch_number" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --batch_number=?*)
            NUM=${1#*=}
            ;;
        -s|--batch_size)
            if [ -n "$2" ]; then
                SIZE=$2
                shift
            else
                printf 'ERROR: "--batch_size" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --batch_size=?*)
            SIZE=${1#*=}
            ;;
        -t|--threads)
            if [ -n "$2" ]; then
                THREADS=$2
                shift
            else
                printf 'ERROR: "--threads" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --threads=?*)
            THREADS=${1#*=}
            ;;
        -f|--feature)
            if [ -n "$2" ]; then
                FEATURE=$2
                shift
            else
                printf 'ERROR: "--feature" requires a non-empty option argument.\n' >&2
                exit 1
            fi
            ;;
        --feature=?*)
            FEATURE=${1#*=}
            ;;
        --)
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)
            break
    esac
    shift
done

# --- BEGIN SCRIPT:
COVAR="../../_m/${FEATURE}/Brainseq_LIBD.cov"
PRE_GEXP="../../_m/${FEATURE}/twas_gene_expression.txt"
OUT_DIR="./WEIGHTS/"
mkdir -p $OUT_DIR/

BATCH_START=`expr $NUM + 1`
BATCH_END=`expr $NUM + $SIZE`
NR="${BATCH_START}"

mkdir --parents ./tmp/$NR
mkdir --parents ./hsq/$NR

# Loop through each gene expression phenotype in the batch
cat $PRE_GEXP | \
    awk -vs=$BATCH_START 'NR == s' |  \
    while read PARAM; do
        # Get the gene positions +/- 500 kb
        GNAME=`echo $PARAM | awk '{print $1}'`
        CHR=`echo $PARAM | awk '{print $2}'`
        P0=`echo $PARAM | awk '{print $3 - 0.5e6}' | awk '{if($1<0)print 0;else print}'`
        P1=`echo $PARAM | awk '{print $3 + 0.5e6}'`
        P1=`awk -vn=$CHR -ve=$P1 '$1=="chr"n {if(e > $2)print $2;else print e}' $CHR_INFO`

        PREFIX=`basename $PRE_GEXP .txt`
        OUT="./tmp/$NR/$PREFIX.$GNAME"

        echo $GNAME $CHR $P0 $P1
        # Pull out the current gene expression phenotype
        echo $PARAM | sed 's/ /\n/g' | tail -n+4 | paste $PRE_GEXP.ID - > $OUT.pheno
        # Get the locus genotypes for all samples and set current gene expression as the phenotype
        $PLINK --bfile $PRE_GENO.$CHR \
               --extract $LDREF/1000G.EUR.$CHR.bim \
               --from-bp $P0 --to-bp $P1 \
               --chr $CHR --mind 0.1 \
               --pheno $OUT.pheno \
               --keep $OUT.pheno \
               --allow-no-sex \
               --make-bed \
               --out $OUT
        # --force-intersect \

        # Process all samples together (for reference purposes only since this is mult-ethnic data)
        FINAL_OUT="$OUT_DIR/$PREFIX.$GNAME"

        $RSCRIPT $FUSION/FUSION.compute_weights.R \
                 --bfile $OUT \
                 --save_hsq \
                 --verbose 2 \
                 --hsq_p 0.01 \
                 --covar $COVAR \
                 --tmp $OUT.tmp \
                 --out $FINAL_OUT \
                 --PATH_gcta $GCTA \
                 --threads $THREADS \
                 --pheno $OUT.pheno \
                 --PATH_gemma $GEMMA \
                 --PATH_plink $PLINK \
                 --models lasso,top1,enet,blup #,bslmm
               # --tmpdir tmp/$NR/ \

        # Append heritability output to hsq file
        cat $FINAL_OUT.hsq >> hsq/$NR.hsq

        # Clean-up just in case
        rm -f $FINAL_OUT.hsq $OUT.tmp.*

        # Remove all intermediate files
        rm $OUT.*
    done
