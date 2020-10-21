#!/bin/bash
#use output info.vcf.gz from prepare_vcf_mlcnv.sh to build a pon
#seperate by gender with provided gender info
#v3 use provided sv vcf for generating truthset

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
svvcf_chk=$(ls *.vcf.gz | grep ALL.wgs | wc -l)
if [ $svvcf_chk -ge 1 ]; then
    sv_file=$(ls *.vcf.gz | grep ALL.wgs)
fi

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


#if sv vcf (1000G) was provide and the sample exist in the list, add truth information to the info.vcf and info.csv
if [ $svvcf_chk -ge 1 ]; then
    for j in male female; do
        ponvcf=$(ls mlcnv_${j}_*.pon.info.vcf.gz)
        for i in $(cat $j.infovcf.list); do
            gvcf_name=${i%%.*}
            svname_chk=$(bcftools query -l $sv_file | grep $gvcf_name | wc -l)
            if [ $svname_chk -ge 1 ]; then
                echo "$gvcf_name exist in $sv_file as $j with $ponvcf"
                bcftools index $sv_file
                bcftools view -s $gvcf_name $sv_file \
                | bcftools norm -m- | bcftools query -f "%CHROM\t%POS\t%INFO/END[\t%GT]\t%ALT\t%INFO/SVTYPE\n" \
                | awk '{if($(NF-2) == "0|0") $NF="NOR"; else {if($(NF) == "CNV"){if($(NF-1) == "<CN0>") $NF="DEL"; else $NF="DUP"}}; print $0}' \
                | awk '{if($3 == ".") $3=$2; print $1"\t"$2"\t"$3"\t"$NF}' | awk '{if($4 >= 50) print}' \
                | bgzip > $gvcf_name.tab.gz
                tabix -s1 -b2 -e3 $gvcf_name.tab.gz
                echo '##INFO=<ID=CNV,Number=1,Type=String,Description="CNV events annotation">' > ann.hdr
                bcftools annotate -a $gvcf_name.tab.gz -h ann.hdr -c CHROM,FROM,TO,CNV $gvcf_name.info.vcf.gz -Oz > $gvcf_name.train.info.vcf.gz
                bcftools index -t -f $gvcf_name.train.info.vcf.gz 
                
                bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\t%majAF\t%CNV\n' $gvcf_name.train.info.vcf.gz | awk -F'\t' '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'> temp.info.bed
                pon_num=$(bcftools query -l $ponvcf | wc -l)
                bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\n' $ponvcf | awk -F'\t' -v var=$pon_num '{print $1"\t"($2-1)"\t"$2"\t"$3/var"\t"$4/var}' > temp.pon.bed 
                bedtools intersect -a temp.info.bed -b temp.pon.bed -wo | awk -F'\t' '{print $1","$3","$4/($11+0.0001)","$5/($12+0.0001)","$6","$7}' > $gvcf_name.train.info.csv
                cp $gvcf_name.train.info.csv $outfolder/$gvcf_name.train.info.csv
                cp $gvcf_name.train.info.vcf.gz $outfolder/$gvcf_name.train.info.vcf.gz
                rm temp*.bed
                
            fi
        done
    done
fi

ls -LR $outfolder

