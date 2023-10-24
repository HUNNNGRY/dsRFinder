#!/bin/bash

rediTab=$1 # REDItools tbl results as input 
chrSize=$2 # chrsize
gnFa=$3 # gn/tx fa
binBed=$4 # gn/tx bin bed
EERNum=$5 # EER num cutoff per bin
binSize=$6 # length/size of bin
extSize=$7 # length/size of extended flank, usually half of binSize
minLen=$8 # min length/size cutoff of merged,extended EER cluster region
maxLen=$9 # max length/size cutoff of merged,extended EER cluster region
tmpDir=$10 # tmp dir

source /BioII/lulab_b/baopengfei/mambaforge/bin/activate REDItools
echo "start getEERcluster at `date`"

## filter EER-cluster bins
### count&filter
TableToGFF.py \
	-i ${rediTab} \
	-b 300000 -T ${tmpDir} \
	-o ${rediTab}.gff
gff2bed < ${rediTab}.gff \
	> ${rediTab}.bed
cut -f1-6 ${rediTab}.bed >  ${rediTab}.bed6


cat ${rediTab}.bed6 \
	| sort -k1,1 -k2,2n \
	| bedtools coverage -s -sorted -counts \
		-a ${binBed} -b - \
	| awk 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,$4,$7,$6}}' \
	| awk -v num=${EERNum} '($5>=num) {print $0}' \
		> ${rediTab}_filter_sort.bed6.count
#-g ${chrSize} \

## (optional: append editing sites position to Tab/bed)
#bedtools intersect ... 
#bedtools merge ... -c 2 -o collapse


## merge&extend EER-cluster bins
#note that extended region has potential editing sites, we neglect for current version
bedtools merge -s -d ${binSize} -c 4,5,6 -o distinct,sum,distinct -i ${rediTab}_filter_sort.bed6.count > ${rediTab}_filter_sort.bed6.count.merge
bedtools slop -s -l ${extSize} -r ${extSize} -g ${chrSize} \
	-i ${rediTab}_filter_sort.bed6.count.merge \
	> ${rediTab}_filter_sort.bed6.count.merge.ext

## filter length & add name
awk -v minL=${minLen} -v maxL=${maxLen} '(($3-$2)>=minL) && (($3-$2)<=maxL) {print $0}' ${rediTab}_filter_sort.bed6.count.merge.ext \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3"_"$6 "\t" $5 "\t" $6}' \
> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen


## RNAfold
### extract fa
#bedtools>=v2.30
bedtools getfasta -nameOnly -s \
	-fi ${gnFa} \
	-bed ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen \
	> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa
echo "end getEERcluster at `date`"

