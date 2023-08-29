#!/bin/bash
set -e

outvcf=$1
vcf=$2
bam=$3
ref_fasta=$4
mapq0perc=$5
outdir=$(dirname "$outvcf")

#grab sites that don't already have the MQ0 field
zgrep -v "^#" "$vcf" | grep -v "MQ0" | cut -f 1,2 | while read chr pos;do
    pysamstats --type mapq --chromosome $chr --start $pos --end $((pos+1)) "$bam"  | grep $pos | cut -f 1,2,5 >>$outdir/mapq0counts
done


if [[ ! -s $outdir/mapq0counts ]];then
    #no sites to process, just make a copy of the vcf. 
    #Have to handle both gzipped and unzipped vcfs
    if [[ ${vcf: -3} == ".gz" ]];then 
        gunzip -c $vcf > $outdir/mapq0.vcf
    else
        cp $vcf $outdir/mapq0.vcf
    fi
else 
    #need to add MQ0
    
    #does the file contain the MQ0 field already?
    mqcount=0
    case "$vcf" in
        *.gz | *.tgz )
            #gzipped vcf
            mqcount=$(gunzip -c "$vcf" | grep "^#" | grep -w MQ0 | wc -l)
            ;;
        *)
            #non-gzipped vcf
            mqcount=$(grep "^#" "$vcf" | grep -w MQ0 | wc -l)
            ;;
    esac
    
    if [[ $mqcount -gt 0 ]];then
        #already has mq0 set, we're all good
        vcf-info-annotator --overwrite -o $outdir/mapq0.vcf "$vcf" $outdir/mapq0counts MQ0 
    else
        #no mq0, need to set the header line as well
        vcf-info-annotator -o $outdir/mapq0.vcf -f Integer -d "Number of MAPQ == 0 reads covering this record" "$vcf" $outdir/mapq0counts MQ0
    fi
fi
#finally, set the filter tags on the vcf
#the multiplication by 1.0 is necessary to convert integers to floats before dividing in the JEXL expression 
gatk VariantFiltration --reference $ref_fasta --output $outdir/mapq_filtered.vcf.gz --variant $outdir/mapq0.vcf --filter-expression "((MQ0*1.0) / (DP*1.0)) > $mapq0perc" --filter-name "MAPQ0"
#java -Xmx6g -jar /opt/GenomeAnalysisTK.jar -T VariantFiltration -R $ref_fasta -o $outdir/mapq_filtered.vcf.gz --variant $outdir/mapq0.vcf --filterExpression "((MQ0*1.0) / (DP*1.0)) > $mapq0perc" --filterName "MAPQ0"
