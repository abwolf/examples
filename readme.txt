# Directory for example files for step by step process to go from msprime simulated data to Neand called bedfiles and desert calls

# NULL DATA
# 1) Generate null simulated vcfs and windows S*-scores : slrun.Sstar.Sriram.null.nomig.sh

$ mkdir null/
$ cd null/
$ for i in {1..5}; do
    sbatch --array=1-1000 slrun.Sstar.Sriram.null.nomig.sh
  done
$ mv *.vcf.* vcfs/
$ mv *.windowcalc_out.gz RegionFiles/
$ ls -C RegionFiles/ | cut -f 3 -d '_' | sort - | uniq - > Sriram.chr_list

# 2) Generate match-pvalues and match-database from null data : write_config_file.sh , run_cetus.sh
# see https://github.com/lparsons/archaic_match & https://bitbucket.org/lance_parsons/abwolf_match_pvalue_simulations/ 
# for setting up match_pvalue_simulations environment and scripts

$ sh write_config_file.sh n1_0.0_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0
$ cd null/
$ source activate match_pvalue_simulations
	## set $configfile=="config_n1_0.0_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.yaml"
$ sh run_cetus.sh &> run_cetus.o


# 3) Generate S*-score ecdfs from null data : slrun.SstarECDFpvalueCalculation.sh

	## comment out "Calculating S*-pvalues" cmd and use only "Generating ECDFs" cmd
	## set $mdl=="Sriram"
$ cd null/
$ sbatch slrun.SstarECDFpvalueCalculation.sh

# A useable null/ directory, at minimum, needs a "null-XX.windows-50000-10000/match_percent_database.db" & an "SstarECDF_maxchrm_XX.RData.gz"

# ADMIXTURE DATA
# 4) Generate simulated vcfs and windowed S*-scores

$ mkdir n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0/
$ cd n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0/
$ sbatch --array=1-1000 slrun.Sstar.Sriram.test.mig.sh
$ mv *.vcf.* vcfs/
$ mv *.windowcalc_out.gz RegionFiles/
$ ls -C RegionFiles/ | cut -f 3 -d '_' | sort - | uniq - > Sriram.chr_list


# 5) Calculate match_pct p-values

$ sh write_config_file.sh n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0
$ cd n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0/
$ source activate match_pvalue_simulations
        ## set $configfile=="config_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.yaml"
$ sh run_cetus.sh &> run_cetus.o


# 6) Generate S*-pvalue for admix data and call Neanderthal haplotypes as bedfiles

	## Use only the "Calculating S*-pvalues" cmd
	## set $mdl=="Sriram
	## set $admix=="n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0"
$ cd n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0/ 
$ sbatch slrun.SstarECDFpvalueCalculation.sh


# 7) Sort and merge bedfiles, and separate EAS/ASN and EUR calls into separate bedfiles

$ for dir in $(ls -C | grep n1_0.1 ); do
    echo $dir
    for i in $(ls -C $dir/bedfiles ); do
        echo $i
        cat $dir/bedfiles/$i | grep -v msp_ID | sort-bed - | bedops --merge - | gzip -c - > $dir/bedfiles/$i.merged.gz
    done
done
$ for dir in $(ls -C | grep n1_0.1 ); do
    echo $dir
    for i in $(ls -C $dir/bedfiles/ | grep ".bed.merged.gz" ); do
        python separate_pop_bedfiles.py $dir/Sriram.popfile $dir/bedfiles/$i
    done
done


# 8) Calculate admixture proportion from S*-calls per population
# The pct_admixture is calculated as the (sum of introgressed bp)/(2 x #_individuals)/(#_replicates)/(simulated_chr_length)

$ for dir in $(ls -C | grep n1_); do 
    echo $dir
    zcat $dir/bedfiles/*.merged_ASN.gz | awk 'BEGIN {OFS="\t"} {sum_bp+=$3-$2} END {print "ASN: " sum_bp/1008/1000/1000000}'
    zcat $dir/bedfiles/*.merged_EUR.gz | awk 'BEGIN {OFS="\t"} {sum_bp+=$3-$2} END {print "EUR: " sum_bp/1006/1000/1000000}'
done > admixture_Sstar.txt


# CALCULATING DESERTS ; Requires simulatied chromosome be >=10Mb
# 9) Combine introgressed bedfiles and split chromosomes in to 5_to_10Mb windows

$ cd n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0/bedfiles/
$ cat *.bed.merged.gz > Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.gz
$ python split_chromosomes.temp.loop.py \
    Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.gz \			# infile name
    Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb.gz \	# outfile name
    10000000											# simulated_chr_length


# 10) Calculate number of windows with admixture proportion below significant threshold

$ zcat Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb.gz \
    | awk 'BEGIN {OFS="\t"} {print $0, "0.1", "0.0"}' \		# add columns for n1_admixture and n2_admixture level
    > Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb
$ Rscript toPct_Int.1_to_15Mb.R \
    Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb \	# infile name
    2014 \										# 2x #_individuals
    Sriram \										# model name
    0.000316 \										# admixture proprotion threshold
    Sriram.chr_list \									# list of chormosomes
    > Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb.prop_blw_thrhld
$ rm Sriram_nonAfr_n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0.bed.merged.5_to_10Mb


