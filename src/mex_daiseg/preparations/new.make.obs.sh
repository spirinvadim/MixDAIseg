#!/bin/bash

CHR=$1
panelfinal=$2

mex=$3
eu=$4
na=$5
af=$6
arch=$7
aa=$8
bed=$9




outtxt=${10}

echo ${arch}
bcftools query -f '%POS\t%REF\t%ALT\n' ${panelfinal}>positions.chr${CHR}.txt

bcftools query -S ${eu} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.eu.txt
echo '...'
bcftools query -S ${na} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.na.txt
echo '...'
bcftools query -S ${af} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's/|//g' |sed 's|/||g' |sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.af.txt

echo '...'
bcftools query -S ${arch} -f '[%GT]\n' ${panelfinal}|awk '{for(i=1;i<=length($0);i++){a[substr($0,i,1)]=1} for(i in a){printf("%s",i)} print "";delete a}'|sed 's|/||g' |sed 's/|//g'|sed -e 's/\.//g'|sed 's/^$/\./' >al.chr${CHR}.arch.txt

echo '...'
bcftools query -S ${mex}  -f '[%GT ]\n'  ${panelfinal} |sed  's/|/ /g'|sed 's|/| |g' >  obs.chr${CHR}.mex.txt

paste  positions.chr${CHR}.txt al.chr${CHR}.eu.txt al.chr${CHR}.na.txt al.chr${CHR}.af.txt al.chr${CHR}.arch.txt obs.chr${CHR}.mex.txt> 3.chr${CHR}.txt


printf '#POSITIONS\t#REF\t#ALT\t#ANCESTRAL\t#EU\t#NA\t#AF\t#ARCHAIC\t#MEX\n'>header.txt
cat header.txt 3.chr${CHR}.txt > temp.allels.ref.and.obs.chr${CHR}.txt

rm header.txt
paste positions.chr${CHR}.txt obs.chr${CHR}.mex.txt > obs.chr${CHR}.txt

rm al.*
rm positions.chr${CHR}.txt 
rm obs.chr${CHR}.mex.txt
rm 3.chr${CHR}.txt


python3 obs2.py ${CHR}  temp.allels.ref.and.obs.chr${CHR}.txt ${aa} ${bed} ${outtxt}

rm temp.allels.ref.and.obs.chr${CHR}.txt














