#!/bin/bash
#use output info.vcf.gz from prepare_vcf_mlcnv.sh to build a pon
#seperate by gender with provided gender info

Inputfileloc="/gatk/data/in"
outfolder="/gatk/data/out"
echo $cmdstr
cd $Inputfileloc
pri_loc=$(dirname "$(readlink -f "$0")")
#gvcf_file=$(ls *.vcf.gz | grep -v ALL.wgs | grep -v dbsnp)
#sv_file=$(ls ALL.wgs.*.vcf.gz)
#dbsnp_file=$(ls dbsnp*.vcf.gz)
#reference=$(ls *.fa* | grep -v ".fai")
#csv2jpg=$(ls *.py)
#bed=$(ls *.bed)
n_proc=$(grep -c ^processor /proc/cpuinfo)
#n_proc=$(awk -v var=$n_proc 'BEGIN{print var/2}')

#chk_pon = $(ls *.vcf.gz | grep pon.vcf.gz | wc -l)
#if [ $chk_pon -ge 1 ]; then pon_vcf=$(ls *.pon.vcf.gz); fi

ls *.info.vcf.gz > infovcf.list
echo "$thread_n cores"

addinfo=$(ls *.tsv)

#reg2cov_py=$(ls reg2cov*.py)
echo "n_proc:$n_proc; reference: $reference"
cat infovcf.list

for i in $(cat infovcf.list); do
    #echo $i 
    #to avoid merge issue with FORMAT:AD, change Number=R to .
    # bcftools view $i | sed 's/ID=AD,Number=R/ID=AD,Number=\./' | bcftools view - -Oz 
    bcftools index -t $i
    name=${i%%.*}
    gender=$(cat $addinfo| awk -v var=$name '{if($1 == var) print $2}')
    echo $name $gender
    if [ $(echo $gender | grep -w 'male' | wc -l) -ge 1 ]; then echo $i >> male.infovcf.list; fi
    if [ $(echo $gender | grep -w 'female' | wc -l) -ge 1  ]; then echo $i >> female.infovcf.list; fi 
done

for i in male female; do
    echo $i
    cat $i.infovcf.list
    echo "build pon with merge, sum the DP relDP and sGQ"
    #build pon with merge, sum the DP and GQ
    bcftools merge -m all -i DP:sum,relDP:sum,sGQ:sum -l $i.infovcf.list --force-samples | bcftools sort - -Oz > temp.pon.info.vcf.gz
    sample_num=$(bcftools query -l temp.pon.info.vcf.gz |wc -l)
    cp temp.pon.info.vcf.gz mlcnv_${i}_${sample_num}.pon.info.vcf.gz 
    bcftools index -t mlcnv_${i}_${sample_num}.pon.info.vcf.gz 
    cp mlcnv_${i}_${sample_num}.pon.info.vcf.gz  $outfolder/mlcnv_${i}_${sample_num}.pon.info.vcf.gz
    cp mlcnv_${i}_${sample_num}.pon.info.vcf.gz.tbi  $outfolder/mlcnv_${i}_${sample_num}.pon.info.vcf.gz.tbi
done
ls -LR $outfolder

