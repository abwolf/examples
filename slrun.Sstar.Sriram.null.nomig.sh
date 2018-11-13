#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=5G
#SBATCH --qos=1hr
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/tmp/abwolf/Sstar/null/slrun.Sstar.Sriram.null.nomig.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

mdl=$( echo Sriram)
seed=$( echo ${SLURM_ARRAY_TASK_ID}${RANDOM} )

n1=$( echo 0.0)
n2=$( echo 0.0)

m_AF_B=$( echo 0.0 )
m_B_AF=$( echo 0.0 )
m_AF_EU=$( echo 0.0 )
m_EU_AF=$( echo 0.0 )

eur=$( echo 1006 )   #1006
asn=$( echo 1008 )   #1008 ; 2040
len=$( awk 'BEGIN {print 1e6}' )

tag=$(echo "$mdl"_nonAfr_"$seed"_n1_"$n1"_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF")
dir=~/SimulatedDemographic/msprime/
sstardir=~/SimulatedDemographic/Sstar/

mkdir -p RegionFiles/
mkdir -p vcfs/
mkdir -p SstarSigFiles/
mkdir -p bedfiles/

echo mdl: $mdl	seed: $seed n1: $n1 n2: $n2 eur: $eur asn: $asn		tag: $tag
echo **RUN MSPRIME SIMULATION AND OUTPUT VCF**
cmd=$( echo " python $dir/bin/msprime.Admixture_Simulate.py \n
			-p nonAfr \n
			-o $mdl \n
			-s $seed \n
			-i 2 \n
			-n $n1 \n
			-d $n2 \n
			-e $eur \n
			-a $asn \n
			-r 216 \n
			--migration_AF_B $m_AF_B \n
			--migration_B_AF $m_B_AF \n
			--migration_AF_AS 0 \n
			--migration_AS_AF 0 \n
			--migration_AF_EU $m_AF_EU \n
			--migration_EU_AF $m_EU_AF \n
			--migration_EU_AS 0 \n
			--migration_AS_EU 0 \n
			-c vcf \n
			-l $len  \n
			| gzip -c - > $tag.vcf.gz")
echo -e $cmd

eval $( echo -e $cmd )

echo fin simulation
date


echo **MODIFY VCF , BGZIP , TABIX**
echo change chr
zcat $tag.vcf.gz | awk 'BEGIN {OFS="\t"} /^#/{print$0} !/^#/{$1="'$seed'" ; print $0}' | bgzip -c > $tag.mod.vcf.gz

echo tabix
tabix -p vcf $tag.mod.vcf.gz

echo remove original
rm $tag.vcf.gz
echo fin
date


echo **EXTRACT ARCHAIC DATA FROM VCF**
vcftools --gzvcf $tag.mod.vcf.gz \
	--keep $sstardir/bin/vcf_keep_archaic.txt \
	--recode --stdout \
	| bgzip -c \
	> $tag.mod.ARCHAIC.vcf.gz
tabix -p vcf $tag.mod.ARCHAIC.vcf.gz

echo fin
date

echo **CALCULATE WINDOWED SSTAR SCORES AND MATCH PCT -- updated Sstar version**
cmd=$( echo " ~/software/anaconda2/bin/python $sstardir/s_star_git/bin/windowed_calculations.py \n
	--vcf-has-illumina-chrnums \n
	-vcfz $tag.mod.vcf.gz \n
	-indf ./$mdl.popfile \n
	-target-pops EUR ASN \n
	-ref-pops AFR \n
	--archaic-vcf $tag.mod.ARCHAIC.vcf.gz \n
	-p 10 \n
	-s-star \n
	-winlen 50000 \n
	-winstep 10000 \n
	-no-pvalues \n
	-range 0 $len \n
	| gzip -c -	> $tag.windowcalc_out.gz" )

echo -e $cmd
eval $( echo -e $cmd )
echo ''
echo fin Sstar calc
date
