#!/bin/bash
#prepare gvcf to vcf for cnv calling 

Inputfileloc="/gatk/data/in"
outfolder="/gatk/data/out"
echo $cmdstr
cd $Inputfileloc
pri_loc=$(dirname "$(readlink -f "$0")")
gvcf_file=$(ls *.vcf.gz | grep -v ALL.wgs | grep -v dbsnp)
dbsnp_file=$(ls dbsnp*.vcf.gz)
reference=$(ls *.fa* | grep -v ".fai")
#bed=$(ls *.bed)
n_proc=$(grep -c ^processor /proc/cpuinfo)
#n_proc=$(awk -v var=$n_proc 'BEGIN{print var/2}')
echo "$thread_n cores"
#reg2cov_py=$(ls reg2cov*.py)
norm_py=$(ls vcftsv_norm*.py)
echo "vcf:$gvcf_file; dbsnp:$dbsnp_file; n_proc:$n_proc; norm_py:$norm_py; reference: $reference"

unzip $(ls gatk*.zip)

vcf_subbed(){
    #subset vcf to certain chr or region
    # vcf_subset a.vcf.gz 1:1000-2000
    local gvcf_file=$1
    local bed_name=$2
    local ref=$3
    #local gvcf_name=$(echo $gvcf_file |sed "s/\.vcf.*$//")
    local out_name=$4
    #local m_proc=$(awk -v var=$n_proc 'BEGIN{printf "%.0f\n", var/8}')
    chk_idx=$(ls $gvcf_file | grep tbi |wc -l)
    if [ $chk_idx -eq 0 ]; then 
        #bcftools index -t -f --threads $n_proc $gvcf_file
        bcftools index -t -f $gvcf_file
    fi
    #bcftools index -t -f $gvcf_file
    #check if bed and vcf is consistence with chr or not 
    chk_vcf_chr=$(bcftools view -h $gvcf_file | grep "^##contig" | grep "chr" |wc -l)
    chk_bed=$(echo $bed_name | grep -v vcf | wc -l)
    if [ $chk_bed -ge 1 ]; then
        chk_bed_chr=$(cat $bed_name | head -n5 | grep "^chr" | wc -l)
        if [ $chk_vcf_chr -ge 1 ]; then 
            if [ $chk_bed_chr -lt 1 ]; then
                echo "add chr to bed to match vcf"
                sed -i 's/^/chr/' $bed_name
            fi
        else 
            if [ $chk_bed_chr -ge 1 ]; then
                echo "remove chr to bed to match vcf"
                sed -i 's/^chr//' $bed_name
            fi 
        fi
        head -n5 $bed_name 
    else
        chk_bed_chr=$(bcftools view -H $bed_name | head -n5 | grep "^chr" | wc -l) 
        if [ $chk_vcf_chr -ge 1 ]; then 
            if [ $chk_bed_chr -lt 1 ]; then
                echo "add chr to dbsnp to match vcf"
                cp $bed_name temp.vcf.gz
                bcftools view temp.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | bcftools view - -Oz > $bed_name
                bcftools index -t -f $bed_name
            fi
        else 
            if [ $chk_bed_chr -ge 1 ]; then
                echo "remove chr to dbsnp to match vcf"
                cp $bed_name temp.vcf.gz
                bcftools view temp.vcf.gz | awk '{gsub(/^chr/,""); print}' | bcftools view - -Oz > $bed_name
                bcftools index -t -f $bed_name
            fi 
        fi
        bcftools view -H $bed_name| head -n5
    fi
    
    #bcftools view -R $bed_name --threads $n_proc $gvcf_file | bcftools convert --gvcf2vcf -f $ref --threads $n_proc - | bcftools sort - -Oz -o $out_name
    #bcftools convert -R $bed_name --gvcf2vcf -f $ref --threads $n_proc $gvcf_file | bcftools sort -m ${m_proc}G - -Oz -o $out_name
    #bcftools convert -R $bed_name --gvcf2vcf -f $ref --threads $n_proc $gvcf_file -Oz -o temp.vcf.gz
    #bcftools convert -R $bed_name --gvcf2vcf -f $ref $gvcf_file -Oz -o temp.vcf.gz
    #bcftools view -R $bed_name $gvcf_file | bcftools convert --gvcf2vcf -f $ref - -Oz -o temp.vcf.gz

    /gatk/data/in/gatk*/gatk GenotypeGVCFs \
    -R $ref \
    -V $gvcf_file \
    -L $bed_name \
    --include-non-variant-sites \
    -O temp.vcf \
    -stand-call-conf 00

    echo "finish gvcf2vcf conversion"
    #bcftools sort temp.vcf.gz -Oz -o $out_name
    bcftools sort temp.vcf -Oz -o $out_name
    #bcftools index -t -f --threads $n_proc $out_name
    bcftools index -t -f $out_name
    rm temp.vcf*
    echo "subset: $out_name"
}

gvcf_name=$(echo $gvcf_file |sed "s/\..*$//")
echo $gvcf_name

#subset gvcf with dbsnp
bcftools query -f "%CHROM\t%POS\n" $dbsnp_file | awk -F'\t' '{print $1"\t"($2-1)"\t"$2}' > $dbsnp_file.bed
echo "vcf_subbed $gvcf_file $dbsnp_file.bed $reference $gvcf_name.bed.vcf.gz"
vcf_subbed $gvcf_file $dbsnp_file $reference $gvcf_name.bed.vcf.gz 
#wait $pid
echo "#check bed.vcf.gz first 10 lines"
bcftools view -H $gvcf_name.bed.vcf.gz | head -n10
cp $gvcf_name.bed.vcf.gz $outfolder/$gvcf_name.dbsnp.vcf.gz
echo "variants number $(bcftools view -H $gvcf_name.bed.vcf.gz | wc -l)"

#echo -e "CHROM\tPOS\tDP\tAF\tAD\tGQ" > $gvcf_name.tsv
bcftools query -f '%CHROM\t%POS\t%DP\t%AF[\t%AD][\t%GQ][,%RGQ]\n' $gvcf_name.bed.vcf.gz > $gvcf_name.tsv


python $norm_py $gvcf_name.tsv
#'CHROM','POS','relDP','maxAF','GQ'

cat $gvcf_name.norm.tsv | sed 1d | bgzip > $gvcf_name.tab.gz
tabix -s1 -b2 -e2 $gvcf_name.tab.gz
echo '##INFO=<ID=relDP,Number=1,Type=Integer,Description="Combined depth across samples, relative to each sample 100">' > ann.hdr
echo '##INFO=<ID=majAF,Number=1,Type=Float,Description="Major Allele Frequency of each samples">' >> ann.hdr
echo '##INFO=<ID=sGQ,Number=1,Type=Integer,Description="Genotype quality per samples,combine both GQ for variant site and RGQ for non-variant site">' >> ann.hdr
#annotate the vcf, change FORMAT:AD Number=R to . to avoid merging issues
bcftools annotate -a $gvcf_name.tab.gz -h ann.hdr -c CHROM,POS,relDP,majAF,sGQ $gvcf_name.bed.vcf.gz \
| sed 's/ID=AD,Number=R/ID=AD,Number=\./' | bcftools sort - -Oz > $gvcf_name.info.vcf.gz
bcftools index -t -f $gvcf_name.info.vcf.gz
rm $gvcf_name.tab.gz* $gvcf_name*.tsv 

cp $gvcf_name.info.vcf.gz $outfolder/$gvcf_name.info.vcf.gz
cp $gvcf_name.info.vcf.gz.tbi $outfolder/$gvcf_name.info.vcf.gz.tbi

ls $outfolder
#ls -LR $outfolder

