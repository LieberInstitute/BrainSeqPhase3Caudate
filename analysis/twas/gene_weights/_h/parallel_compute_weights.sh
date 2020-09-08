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
EOF
}

set -o errexit # stop execution on error

COVAR="../../_m/genes/Brainseq_LIBD.cov"
PLINK="/users/jbenjami/software/bin/plink"
PRE_GEXP="../../_m/genes/twas_gene_expression.txt"
FUSION="/dcl02/lieber/apaquola/opt/fusion_twas-master"
GCTA="/dcl02/lieber/apaquola/opt/fusion_twas-master/gcta64"
GEMMA="/dcl02/lieber/apaquola/opt/fusion_twas-master/gemma-0.98.1"
RSCRIPT="/dcl01/lieber/apaquola/trash/conda/env/my_svnR-3.5/R/3.5/bin/Rscript"
CHR_INFO="/dcl01/lieber/apaquola/genomes/human/gencode26/GRCh38.p10.ALL/chromInfo/_m/chromInfo.txt"
PRE_GENO="/dcl02/lieber/apaquola/projects/phase3_paper/data_definition/machineLearning/jaffe_twas/byChrom/_m/Brainseq_LIBD"
LDREF="/dcl02/lieber/apaquola/projects/phase3_paper/data_definition/machineLearning/jaffe_twas/reference_hg38/_m/LDREF_hg38"

# Reset all variables that might be set
NUM=0
SIZE=100
THREADS=1

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
OUT_DIR="./WEIGHTS/"
mkdir -p $OUT_DIR/

BATCH_START=`expr $NUM + 2`
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
	       --make-bed \
	       --out $OUT

	# Process all samples together (for reference purposes only since this is mult-ethnic data)
	FINAL_OUT="$OUT_DIR/$PREFIX.$GNAME"

        $RSCRIPT $FUSION/FUSION.compute_weights.R \
	    	 --bfile $OUT \
	    	 --save_hsq \
	    	 --verbose 2 \
		 --hsq_p 0.01 \
	    	 --tmp $OUT.tmp \
		 --covar $COVAR \
		 --pheno $OUT.pheno \
	    	 --out $FINAL_OUT \
	    	 --PATH_gcta $GCTA \
		 --PATH_gemma $GEMMA \
		 --PATH_plink $PLINK \
		 --threads $THREADS \
	    	 --models lasso,top1,enet,blup #,bslmm
	    
	# Append heritability output to hsq file
	cat $FINAL_OUT.hsq >> hsq/$NR.hsq

	# Clean-up just in case
	rm -f $FINAL_OUT.hsq $OUT.tmp.*
	    
	# Remove all intermediate files
	rm $OUT.*
    done
