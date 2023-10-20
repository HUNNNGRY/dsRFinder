------------------------------------
# prepare environment  & reference
------------------------------------
## prepare&activate conda environment 
ssh baopengfei@101.6.121.24 # hub-remote

```bash
cd /BioII/lulab_b/baopengfei/shared_reference/REDItools/testREDItools

### REDItools version>1.2.1 env
#20231010
git clone https://github.com/BioinfoUNIBA/REDItools 
mamba env create --name REDItools --file /BioII/lulab_b/baopengfei/gitsoft/REDItools/mamba_env.yml

### RNAfold env
mamba create -n py37 python=3.7 seaborn multiprocess

### IntaRNA v3.3.2 env
mamba create -n IntaRNA -c conda-forge -c bioconda intarna
```
## prepare genome/transcriptome bin bed
binSize=50
pre='/BioII/lulab_b/baopengfei/projects/WCHSU-FTC'

### genome
chrSize="${pre}/../exOmics/DNA-seq/genome/chrom.size"
binBed="${pre}/hg38.bins.${binSize}.bed"

bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}.tmp
#add strand
awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "-"}' ${binBed}.tmp >  ${binBed}.tmp2
cat ${binBed}.tmp ${binBed}.tmp2 | sort -k1,1 -k2,2n > ${binBed}


### transcriptome 
chrSize='${pre}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_sort_uniq_newTxID'
binBed='${pre}/exSeek-dev/genome/hg38/tbed/transcriptome_sort_uniq_newTxID.bins.${binSize}.bed6'

cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID  | grep -v "^chr" | grep -v "_miR_" | grep -v "^hsa____">  ${chrSize}
bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}


## prepare filtering reference
### download reference
#### SNP
option1: All_20180418.vcf.gz (too large): 
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz
option2: common_all_20180418.vcf.gz (might be better):
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz 
#convert chr naming from ensembl to UCSC
zcat common_all_20180418.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | gzip -c > common_all_20180418_chr.vcf.gz

#### seq_artifacts+germline_SNV
somatic-hg38_1000g_pon.hg38.vcf.gz: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38

### SNP & seq_artifacts+germline_SNV vcf2bed (only needed for 1st run)
```bash
SNPDir="/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/SNV_ref"
alias bcftools="/BioII/lulab_b/baopengfei/anaconda3/bin/bcftools"
alias vcf2bed="/BioII/lulab_b/baopengfei/gitsoft/bedops_v2.4.40/vcf2bed"
cd $SNPDir

#### SNP
#keep snp only
SNPpre="common_all_20180418_chr"
SNPbed="${SNPDir}/${SNPpre}_snp_sort_addStrandfilterChr.hg38.bed"
bcftools view -v snps ${SNPDir}/${SNPpre}.vcf.gz | perl -lane 'BEGIN {srand(1984)} if (/^#/) { print } elsif (length($F[3]) == 1) { if (rand(1) > 0.5) {print} }' | bgzip > ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz
gzip -dc ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz > ${SNPDir}/${SNPpre}_snp.hg38.vcf

#opt1: vcf2bed, add strand, rm non-auto-chr, rm dup rows 
#split large file like All_20180418
split -l 1000000  ${SNPDir}/${SNPpre}_snp.hg38.vcf  ${SNPDir}/${SNPpre}_snp.hg38.vcf.split 
for i in `ls  ${SNPDir}/${SNPpre}_snp.hg38.vcf.split*`; do echo $i; vcf2bed < $i > ${i}.bed; done
cat ${SNPDir}/${SNPpre}_snp.hg38.vcf.split*.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' | grep "^chr" | grep -vE "alt|random|decoy|chrUn" | uniq > ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.+

#opt2: vcf2bed, add strand, rm non-auto-chr, rm dup rows 
#small file like common_all_20180418
vcf2bed < ${SNPDir}/${SNPpre}_snp.hg38.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' | grep "^chr" | grep -vE "alt|random|decoy|chrUn" | uniq > ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.+


awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-"}' ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.+ > ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.-
cat ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.{+,-} | sort -k1,1 -k2,2n > ${SNPbed}
rm ${SNPDir}/${SNPpre}_snp_sort.hg38.bed.{+,-}
rm ${SNPDir}/${SNPpre}_snp.hg38.vcf


#### seq_artifacts+germline_SNV
#keep snp only
ARTpre="somatic-hg38_1000g_pon"
ARTbed="${SNPDir}/${ARTpre}_snp_sort_addStrandfilterChr.hg38.bed"
bcftools view -v snps ${SNPDir}/${ARTpre}.hg38.vcf.gz | perl -lane 'BEGIN {srand(1984)} if (/^#/) { print } elsif (length($F[3]) == 1) { if (rand(1) > 0.5) {print} }' | bgzip > ${SNPDir}/${ARTpre}_snp.hg38.vcf.gz
gzip -dc ${SNPDir}/${ARTpre}_snp.hg38.vcf.gz > ${SNPDir}/${ARTpre}_snp.hg38.vcf

#vcf2bed, add strand, rm non-auto-chr, rm dup rows
vcf2bed < ${SNPDir}/${ARTpre}_snp.hg38.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' | grep "^chr" | grep -vE "alt|random|decoy|chrUn" | uniq > ${ARTbed}_+
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-"}' ${ARTbed}_+ > ${ARTbed}_-
cat ${ARTbed}_{+,-} | sort -k1,1 -k2,2n > ${ARTbed}
rm ${ARTbed}_{+,-}
rm ${SNPDir}/${ARTpre}_snp.hg38.vcf
```


### convert bed to gff/gtf (gn-mode-only)
```bash
#### SNP
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${SNPbed} \
	-o ${SNPbed}.gtf
	
#### seq_artifacts+germline_SNV
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${ARTbed} \
	-o ${ARTbed}.gtf
```

### convert gtf to gff.gz (gn-mode-only)
```bash
tmpDir='tmp'
ulimit -n 4000

#### SNP
GFFtoTabix.py -i ${SNPbed}.gtf -b 300000 -t $tmpDir > $tmpDir/log1 2>&1 &
rm ${SNPbed}.gtf 

#### seq_artifacts+germline_SNV
GFFtoTabix.py -i ${ARTbed}.gtf -b 300000 -t $tmpDir > $tmpDir/log2 2>&1 &
rm ${ARTbed}.gtf 
```

### convert gn coordinate to tx coordinate (tx-mode-only)
```bash
refPre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38"

#### SNP
#split -d -l 200000  ${SNPbed}  ${SNPbed}. # total: 33630436

/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores ${cores} \
	-i ${SNPbed} \
	--type RNA --mode 1toN \
	--txRef ${refPre}/chrom_sizes/tx_gn_length_newTxID.txt \
	--bed12 ${refPre}/bed/11RNA_map_bed12_newTxID.bed \
	--blockRef ${refPre}/bed/11RNA_map_block_newTxID.txt \
	-o ${SNPbed}.tx.RNA > ${SNPbed}.tx.RNA.log 2>&1 &
#1/50 converted, need add 4DNA gn2tx

/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores ${cores} \
	-i ${SNPbed} \
	--type DNA --mode 1toN \
	--txRef ${refPre}/chrom_sizes/tx_gn_length_newTxID.txt \
	--bed12 ${refPre}/bed/long_DNA_newTxID.bed \
	-o ${SNPbed}.tx.DNA > ${SNPbed}.tx.DNA.log 2>&1 &
cat ${SNPbed}.tx.{RNA,DNA} | sort -k1,1 -k2,2n > ${SNPbed}.tx
rm  ${SNPbed}.tx.{RNA,DNA} 


#### seq_artifacts+germline_SNV
#split -d -l 200000  ${ARTbed}  ${ARTbed}. # total: 2262418

/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores 14 \
	-i ${ARTbed} \
	--type RNA --mode 1toN \
	--txRef ${refPre}/chrom_sizes/tx_gn_length_newTxID.txt \
	--bed12 ${refPre}/bed/11RNA_map_bed12_newTxID.bed \
	--blockRef ${refPre}/bed/11RNA_map_block_newTxID.txt \
	-o ${ARTbed}.tx.RNA > ${ARTbed}.tx.RNA.log 2>&1 &
/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores ${cores} \
	-i ${ARTbed} \
	--type DNA --mode 1toN \
	--txRef ${refPre}/chrom_sizes/tx_gn_length_newTxID.txt \
	--bed12 ${refPre}/bed/long_DNA_newTxID.bed \
	-o ${ARTbed}.tx.DNA > ${ARTbed}.tx.DNA.log 2>&1 &
cat ${ARTbed}.tx.{RNA,DNA} | sort -k1,1 -k2,2n > ${ARTbed}.tx
rm  ${ARTbed}.tx.{RNA,DNA} 
```

### convert bed to gff/gtf (tx-mode-only)
```bash
#### SNP
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${SNPbed}.tx \
	-o ${SNPbed}.tx.gtf
	
#### seq_artifacts+germline_SNV
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${ARTbed}.tx \
	-o ${ARTbed}.tx.gtf
```

### convert gtf to gff.gz (tx-mode-only)
```bash
#### SNP
GFFtoTabix.py -b 300000 -i ${SNPbed}.tx.gtf > $tmpDir/log1 2>&1 &
rm ${SNPbed}.tx.gtf

#### seq_artifacts+germline_SNV
GFFtoTabix.py -b 300000 -i ${ARTbed}.tx.gtf > $tmpDir/log2 2>&1 &
rm ${ARTbed}.tx.gtf
```




------------------------------------
dsRNAfinder (genome mode)
------------------------------------
## run EM on STAR bam (holding)
ssh hub
```bash
#current cfRNA-SEEK total pipeline only report 1 best align record for multi-align reads, if you need further quantification on repeats dsRNA, not just determination of dsRNA positions, we need remap bam using flag: --outSAMmultNmax -1 --outFilterMultimapNmax 100, and then run EM-rescue below:
rawBam="~/tmp/20230523-S048-1_BC11/hg38_v38_sortbyCoord.bam"
emDir="~/tmp/20230523-S048-1_BC11/CLAM"
source ~/anaconda3/bin/activate clam
python /BioII/lulab_b/baopengfei/gitsoft/CLAM-1.2.0/run/preprocessor.py \
	-i ${rawBam} \
	-o ${emDir}
#Unique reads = 20853688;  Multi reads = 3831753 (15.52 %)
samtools merge -f \
	${emDir}/merged.bam \
	${emDir}/multi.sorted.bam \
	${emDir}/unique.sorted.bam
samtools sort -@ 4 ${emDir}/merged.bam > ${emDir}/merged.sorted.bam
samtools index -@ 4 ${emDir}/merged.sorted.bam
rm ${emDir}/merged.bam
```

## run REDItools
```bash
source activate REDItools

cores=10
tmpDir='tmp'
sample_id='test'
inBam='/BioII/lulab_b/baopengfei/tmp/20230523-S048-1_BC11/hg38_v38_sortbyCoord.bam'
gnFa='/BioII/lulab_b/baopengfei/shared_reference/hg38/genome.fa'
outDir='output/SLE/REDItoolDnaRna'

mkdir -p $outDir
REDItoolDnaRna.py -t ${cores} \
	-i ${inBam} \
	-f ${gnFa} \
	-o ${outDir} \
	> ${outDir}/log 2>&1 &
#40285769 records
```

## test dsRNAfinder-downstream
```bash
### rm SNP in REDItools tab
#extract pid
pid=`head -n3 ${outDir}/log | grep "Analysis ID" | sed s/"Analysis ID: "/""/g`
rediTab="../output/GSE71008/REDItoolDnaRna/outTable_${pid}"
SNPgff='${SNPbed}.tx.gff.gz'
FilterTable.py \
	-i ${rediTab} \
	-s ${SNPgff} \
	-S peak \
	-o ${rediTab}_tmp \
	-E -p
#6/96877 rm


## filter seq_artifacts+germline_SNV sites
ARTgff='${ARTbed}.tx.gff.gz'
FilterTable.py \
	-i ${rediTab}_tmp \
	-s ${ARTgff} \
	-S peak \
	-o ${rediTab}_filter \
	-E -p
#?/? rm



## filter EER-cluster bins
### count&filter
TableToGFF.py \
	-i ${rediTab}_filter \
	-s -t \
	-o ${rediTab}_filter.gff
gff2bed < ${rediTab}_filter.gff \
	> ${rediTab}_filter_sort.bed
cut -f1-6 ${rediTab}_filter_sort.bed >  ${rediTab}_filter_sort.bed6

cat ${rediTab}_filter_sort.bed6 \
            | bedtools sort \
	   | sort -k1,1 -k2,2n \
            | bedtools coverage -s -sorted -counts \
		-g ${chrSize} \
                -a ${binBed} -b - \
            | awk 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,$4,$7,$6}}' \
	| awk -v num=${EERNum} '($5>=num) {print $0}' \
                > ${rediTab}_filter_sort.bed6.count

(### optional: append editing sites position to Tab/bed)
#bedtools intersect ... 
#bedtools merge ... -c 2 -o collapse


## merge&extend EER-cluster bins
#note that extended region has potential editing sites, we neglect for current version
extSize=25
bedtools merge -s -d ${binSize} -c 4,5,6 -o distinct,sum,distinct -i ${rediTab}_filter_sort.bed6.count > ${rediTab}_filter_sort.bed6.count.merge
bedtools slop -s -l ${extSize} -r ${extSize} -g ${chrSize} \
	-i ${rediTab}_filter_sort.bed6.count.merge \
	> ${rediTab}_filter_sort.bed6.count.merge.ext

## filter length & add name
minLen=50
maxLen=3000
awk -v minL=${minLen} -v maxL=${maxLen} '(($3-$2)>=minL) && (($3-$2)<=maxL) {print $0}' ${rediTab}_filter_sort.bed6.count.merge.ext \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3"_"$6 "\t" $5 "\t" $6}' \
> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen


## RNAfold
### extract fa
bedtools getfasta -nameOnly -s \
	-fi ${gnFa} \
	-bed ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen \
	> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa
	
### get MFE FDR
source activate py37
python3 \
	${pre}/exSeek-dev/scripts/rnafold_dinushuffle_parallel.py  \
	${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa 200 1234 ${inFile}.csv
rm ${inFile}_perm


## convert tx to gn (tx-only)
gnBed=${rediTab}_filter_sort.bed6.count.merge.ext.filterLen_gn.bed
mv ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen $gnBed

(## optional: filter tRNA,miRNA?)


## todo: annotation
/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/scripts/lulab/sequential.assign.long.sh
```

## todo: inter-molecular
```bash
source activate IntaRNA

IntaRNA --threads 4 -q query.fa -t target.fa > outpre
```




------------------------------------
dsRNAfinder (transcript mode)
------------------------------------
## modify REDItools 
vi /BioII/lulab_b/baopengfei/shared_reference/REDItools/REDItoolDnaRna.py
#special chr ID: ALR/Alpha__chr5___45962494____45970456_pos et al.
#lead to REDItoolDnaRna.py when creating tmp dir using chr name
#need modify REDItoolDnaRna.py for tx Bam:
        if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))

#delete empty dir
def remove_empty_folders(path_abs):
    walk = list(os.walk(path_abs))
    for path, _, _ in walk[::-1]:
        if len(os.listdir(path)) == 0:
            os.remove(path)
remove_empty_folders(outfolder)

#~844th row, add abelow to reduce iteration range
chr_list={}
min_count=10 
# bamfile = pysam.AlignmentFile(bamfile, "rb")
chromosomes = bamfile.references
# Iterate over each chromosome in the BAM file
for chrom in dicregions.keys():
    # Get the count of reads for the current chromosome
    count = bamfile.count(chrom)
    # If the count is greater than or equal to the minimum threshold, keep the chromosome
    if count >= min_count:
        # Iterate over each read in the current chromosome and write it to the output file
        chr_list[chrom] = count
# bamfile.close()
chrs=[x for x in dicregions.keys() if (x not in nochrs) and chr_list[x]]

cp /BioII/lulab_b/baopengfei/shared_reference/REDItools/REDItoolDnaRna.py /BioII/lulab_b/baopengfei/mambaforge/envs/REDItools/bin/REDItoolDnaRna.py
chmod 755 /BioII/lulab_b/baopengfei/mambaforge/envs/REDItools/bin/REDItoolDnaRna.py


## run REDItools
```bash
source activate REDItools

cores=8
tmpDir='tmp'
sample_id='REDItoolDnaRna'
inBam="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/tbam/NCpool/bam-EM/merge19_sort/merged.sorted.bam"
#inBam="./merge19_sort_subset01.bam"
#inBam="./testREDItools/tx_19.bam"
gnFa="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa"
#outDir="output/GSE71008_19/REDItoolDnaRna"
outDir="output/GSE71008_19/${sample_id}"
mkdir -p $outDir


REDItoolDnaRna.py -t ${cores} \
	-i ${inBam} \
	-f ${gnFa} \
	-o ${outDir} \
	> ${outDir}/log 2>&1 &
#?h, ? records
```

## test dsRNAfinder-downstream
```bash
### rm SNP in REDItools tab
#extract pid
pid=`head -n3 ${outDir}/log | grep "Analysis ID" | sed s/"Analysis ID: "/""/g`
rediTab="${outDir}/DnaRna_${pid}/outTable_${pid}"
SNPgff="${SNPbed}.tx.sorted.gff.gz"
FilterTable.py \
	-i ${rediTab} \
	-s ${SNPgff} \
	-S peak \
	-o ${rediTab}_tmp \
	-E -p
#6/96877 rm


## filter seq_artifacts+germline_SNV sites
ARTgff="${ARTbed}.tx.sorted.gff.gz"
FilterTable.py \
	-i ${rediTab}_tmp \
	-s ${ARTgff} \
	-S peak \
	-o ${rediTab}_filter \
	-E -p
#?/? rm



## filter EER-cluster bins
### count&filter
TableToGFF.py \
	-i ${rediTab}_filter \
	-s -t \
	-o ${rediTab}_filter.gff
gff2bed < ${rediTab}_filter.gff \
	> ${rediTab}_filter_sort.bed
cut -f1-6 ${rediTab}_filter_sort.bed >  ${rediTab}_filter_sort.bed6

cat ${rediTab}_filter_sort.bed6 \
            | bedtools sort \
	   | sort -k1,1 -k2,2n \
            | bedtools coverage -s -sorted -counts \
		-g ${chrSize} \
                -a ${binBed} -b - \
            | awk 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,$4,$7,$6}}' \
	| awk -v num=${EERNum} '($5>=num) {print $0}' \
                > ${rediTab}_filter_sort.bed6.count

(### optional: append editing sites position to Tab/bed)
#bedtools intersect ... 
#bedtools merge ... -c 2 -o collapse


## merge&extend EER-cluster bins
#note that extended region has potential editing sites, we neglect for current version
extSize=25
bedtools merge -s -d ${binSize} -c 4,5,6 -o distinct,sum,distinct -i ${rediTab}_filter_sort.bed6.count > ${rediTab}_filter_sort.bed6.count.merge
bedtools slop -s -l ${extSize} -r ${extSize} -g ${chrSize} \
	-i ${rediTab}_filter_sort.bed6.count.merge \
	> ${rediTab}_filter_sort.bed6.count.merge.ext

## filter length & add name
minLen=50
maxLen=3000
awk -v minL=${minLen} -v maxL=${maxLen} '(($3-$2)>=minL) && (($3-$2)<=maxL) {print $0}' ${rediTab}_filter_sort.bed6.count.merge.ext \
| awk '{print $1 "\t" $2 "\t" $3 "\t" $1":"$2"-"$3"_"$6 "\t" $5 "\t" $6}' \
> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen


## RNAfold
### extract fa
bedtools getfasta -nameOnly -s \
	-fi ${gnFa} \
	-bed ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen \
	> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa
	
### get MFE FDR
source activate py37 \
	${pre}/exSeek-dev/scripts/rnafold_dinushuffle_parallel.py  \
	${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa 200 1234 ${inFile}.csv
rm ${inFile}_perm


## convert tx to gn (tx-only)
txBed=${rediTab}_filter_sort.bed6.count.merge.ext.filterLen
gnBed=${rediTab}_filter_sort.bed6.count.merge.ext.filterLen_gn.bed
{{
  grep -v '^chr' ${txBed} | $pre/exSeek-dev/bin/tbed2gbed <(cat  $pre/exSeek-dev/genome/hg38/bed/{long_DNA,long_RNA,tRNA,pri_miRNA,piRNA,rRNA}_newTxID.bed) /dev/stdin /dev/stdout
  awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${txBed}
}} > ${gnBed}


(## optional: filter tRNA,miRNA?)
```

## todo: inter-molecular
```bash
source activate IntaRNA

IntaRNA --threads 4 -q query.fa -t target.fa > outpre
```

