#!/bin/bash

set -o pipefail
set -o errexit
set -o nounset

vcfs=$1
outdir=$9

#Modify formatting of vcfs from different tools
/usr/bin/python3 /usr/bin/modify_VCF.py -i ${vcfs} -o ${outdir}

#Outputs from modify_VCF.py
echo ${outdir}/0.mod.vcf > vcf_list.txt
echo ${outdir}/1.mod.vcf >> vcf_list.txt
echo ${outdir}/2.mod.vcf >> vcf_list.txt

#for i in $(echo ${1} | tr "," "\n"); do
#  echo "$i" >> vcf_list.txt
#done

#Run SURVIVOR
/usr/bin/SURVIVOR merge vcf_list.txt ${2} ${3} ${4} ${5} ${6} ${7} ${8}

#Modify format of output
/usr/bin/python3 /usr/bin/modify_SURVIVOR.py -i ${8} > mod.${8}

