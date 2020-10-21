#!/bin/bash
#prepare gvcf to vcf for cnv calling 

#v2 if provide pon do normalization and generate csv; it can also start from eitehr the dbsnp.vcf.gz or the info.vcf.gz
#v3 create an additional train.csv if a sample is found in the 1000G vcf 

Inputfileloc="/gatk/data/in"
outfolder="/gatk/data/out"
echo $cmdstr
cd $Inputfileloc
pri_loc=$(dirname "$(readlink -f "$0")")
gvcf_file=$(ls *.vcf.gz | grep -v ALL.wgs | grep -v dbsnp)
dbsnp_file=$(ls dbsnp*.vcf.gz)
reference=$(ls *.fa* | grep -v ".fai")
#if cnvvcf exist , not need to convert gvcf to vcf again
cnvvcf_chk=$(ls *.vcf.gz | grep dbsnp.vcf.gz |wc -l)
if [ $cnvvcf_chk -ge 1 ]; then 
    cnvvcf=$(ls *.vcf.gz | grep dbsnp.vcf.gz)
else 
    unzip $(ls gatk*.zip)
fi

infovcf_chk=$(ls *.vcf.gz | grep info.vcf.gz | grep -v pon.info.vcf.gz | wc -l)
if [ $infovcf_chk -ge 1 ]; then 
    infovcf=$(ls *.vcf.gz | grep info.vcf.gz | grep -v pon.info.vcf.gz)
fi

ponvcf_chk=$(ls *.vcf.gz | grep pon.info.vcf.gz |wc -l)
if [ $ponvcf_chk -ge 1 ]; then 
    ponvcf=$(ls *.vcf.gz | grep pon.info.vcf.gz)
fi

svvcf_chk=$(ls *.vcf.gz | grep ALL.wgs | wc -l)
if [ $svvcf_chk -ge 1 ]; then
    sv_file=$(ls *.vcf.gz | grep ALL.wgs)
fi

#bed=$(ls *.bed)
n_proc=$(grep -c ^processor /proc/cpuinfo)
#n_proc=$(awk -v var=$n_proc 'BEGIN{print var/2}')
echo "$thread_n cores"
#reg2cov_py=$(ls reg2cov*.py)
norm_py=$(ls vcftsv_norm*.py)
echo "vcf:$gvcf_file; dbsnp:$dbsnp_file; n_proc:$n_proc; norm_py:$norm_py; reference: $reference; cnvvcf_chk:$cnvvcf_chk $cnvvcf; ponvcf_chk:$ponvcf_chk $ponvcf; infovcf_chk:$infovcf; svvcf_chk:$sv_file"


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

if [ $infovcf_chk -ge 1 ]; then 
    gvcf_name=$(echo $infovcf| sed 's/\.info.vcf.gz//')
    echo "$gvcf_name.info.vcf.gz exist"
else
    if [ $cnvvcf_chk -ge 1 ]; then
        gvcf_name=$(echo $cnvvcf| sed 's/\.dbsnp.vcf.gz//')
        echo "$gvcf_name.dbsnp.vcf.gz exist"
        cp $gvcf_name.dbsnp.vcf.gz $gvcf_name.bed.vcf.gz
    else
        #subset gvcf with dbsnp
        bcftools query -f "%CHROM\t%POS\n" $dbsnp_file | awk -F'\t' '{print $1"\t"($2-1)"\t"$2}' > $dbsnp_file.bed
        echo "vcf_subbed $gvcf_file $dbsnp_file.bed $reference $gvcf_name.bed.vcf.gz"
        vcf_subbed $gvcf_file $dbsnp_file $reference $gvcf_name.bed.vcf.gz 
        #wait $pid
        echo "#check bed.vcf.gz first 10 lines"
        bcftools view -H $gvcf_name.bed.vcf.gz | head -n10
        cp $gvcf_name.bed.vcf.gz $outfolder/$gvcf_name.dbsnp.vcf.gz
        echo "variants number $(bcftools view -H $gvcf_name.bed.vcf.gz | wc -l)"

    fi

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
fi

if [ $ponvcf_chk -ge 1 ]; then 
    bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\t%majAF\n' $gvcf_name.info.vcf.gz | awk -F'\t' '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5}'> temp.info.bed
    pon_num=$(bcftools query -l $ponvcf | wc -l)
    bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\n' $ponvcf | awk -F'\t' -v var=$pon_num '{print $1"\t"($2-1)"\t"$2"\t"$3/var"\t"$4/var}' > temp.pon.bed 
    bedtools intersect -a temp.info.bed -b temp.pon.bed -wo | awk -F'\t' '{print $1","$3","$4/($10+0.0001)","$5/($11+0.0001)","$6}' > $gvcf_name.info.csv
    #paste -d',' temp.info.csv temp.pon.csv | awk -F',' '{print $1","$2","$3/$6","$4/$7","$5}'
    cp $gvcf_name.info.csv $outfolder/$gvcf_name.info.csv
    #rm temp*.bed
fi

#if sv vcf (1000G) was provide and the sample exist in the list, add truth information to the info.vcf and info.csv
if [ $svvcf_chk -ge 1 ]; then
    svname_chk=$(bcftools query -l $sv_file | grep $gvcf_name | wc -l)
    if [ $svname_chk -ge 1 ]; then
        echo "$gvcf_name exist in $sv_file"
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
        
        if [ $ponvcf_chk -ge 1 ]; then 
            bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\t%majAF\t%CNV\n' $gvcf_name.train.info.vcf.gz | awk -F'\t' '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'> temp.info.bed
            #pon_num=$(bcftools query -l $ponvcf | wc -l)
            #bcftools query -f '%CHROM\t%POS\t%relDP\t%sGQ\n' $ponvcf | awk -F'\t' -v var=$pon_num '{print $1"\t"($2-1)"\t"$2"\t"$3/var"\t"$4/var}' > temp.pon.bed 
            bedtools intersect -a temp.info.bed -b temp.pon.bed -wo | awk -F'\t' '{print $1","$3","$4/($11+0.0001)","$5/($12+0.0001)","$6","$7}' > $gvcf_name.train.info.csv
            cp $gvcf_name.train.info.csv $outfolder/$gvcf_name.train.info.csv
            cp $gvcf_name.train.info.vcf.gz $outfolder/$gvcf_name.train.info.vcf.gz
            rm temp*.bed
        fi
    fi
fi

ls $outfolder
#ls -LR $outfolder

