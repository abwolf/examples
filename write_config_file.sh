
mdl=$( echo Sriram )
data_dir=/Genomics/akeylab/abwolf/SimulatedDemographic/Sstar/chr1_variable_ref/simulations
admix=$( echo $1 )
chromsize=$( echo 1000000)
configfile=$( echo config_$admix.yaml)

echo '
 window_size: "50000"
 window_step_size: "10000"
 informative_site_method: "derived_in_archaic"
 informative_site_range: "0"

 null_dataset:
     name: "null"
     directory: "'$data_dir'/'$mdl'/null"
     # VCF Info
     vcf_file_pattern: "vcfs/{sample}.mod.vcf.gz"
     chrom_sizes: "'$chromsize'"
     # Population info
     population_file: "'$mdl'.popfile"
     archaic_populations:
         - "Neand1"
     modern_populations:
         - "EUR"
         - "ASN"

 test_dataset:
     name: "'$mdl'_'$admix'"
     directory: "'$data_dir'/'$mdl'/'$admix'"
     # VCF Info
     vcf_file_pattern: "vcfs/{sample}.mod.vcf.gz"
     chrom_sizes: "'$chromsize'"
     # Population info
     population_file: "'$mdl'.popfile"
     archaic_populations:
         - "Neand1"
     modern_populations:
         - "EUR"
         - "ASN"
' \
> $configfile
