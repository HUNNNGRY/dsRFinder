------------------------------------
# prepare environment  & reference
------------------------------------
(below only needed for 1st run)

ssh baopengfei@101.6.121.24 # hub-remote

## global config variable
```bash
pre='/BioII/lulab_b/baopengfei/projects/WCHSU-FTC'
SNPDir="/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/SNV_ref" # ~ARTDir
SNPpre="common_all_20180418_chr"
SNPbed="${SNPDir}/${SNPpre}_snp_sort_addStrandfilterChr.hg38.bed"
ARTpre="somatic-hg38_1000g_pon"
ARTbed="${SNPDir}/${ARTpre}_snp_sort_addStrandfilterChr.hg38.bed"
cores=6
tmpDir="tmp"
```

## prepare&activate conda environment 
### REDItools 
```bash
#v>1.2.1 (20231010)
cd /BioII/lulab_b/baopengfei/shared_reference/REDItools/testREDItools
git clone https://github.com/BioinfoUNIBA/REDItools 
mamba env create --name REDItools --file /BioII/lulab_b/baopengfei/gitsoft/REDItools/mamba_env.yml
```

### modify REDItools
```bash
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
#Iterate over each chromosome in the BAM file
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
```

	• raw: /BioII/lulab_b/baopengfei/gitsoft/REDItools/main/REDItoolDnaRna.py (github)
	• new: /BioII/lulab_b/baopengfei/mambaforge/envs/REDItools/bin/REDItoolDnaRna.py (self-optimized)


### RNAfold 
```bash
#v2.4.14
mamba create -n py37 python=3.7 seaborn multiprocess
```

### IntaRNA
```bash
#v3.3.2 
mamba create -n IntaRNA -c conda-forge -c bioconda intarna
```

## prepare genome/transcriptome bin bed
```bash
### genome
binSize=50
chrSize="${pre}/../exOmics/DNA-seq/genome/chrom.size"
binBed="${pre}/../exOmics/DNA-seq/genome/hg38.bins.${binSize}.bed"

bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}.tmp
#add strand
awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "-"}' ${binBed}.tmp >  ${binBed}.tmp2
cat ${binBed}.tmp ${binBed}.tmp2 | sort -k1,1 -k2,2n > ${binBed}
rm ${binBed}.tmp*

### transcriptome 
bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}
```

## prepare filtering reference
```bash
### download reference
#### SNP
option1: All_20180418.vcf.gz (too large): 
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/All_20180418.vcf.gz
option2: common_all_20180418.vcf.gz (this reduced set might be better):
wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz 
#convert chr naming from ensembl to UCSC
zcat common_all_20180418.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | gzip -c > $SNPDir/common_all_20180418_chr.vcf.gz

#### seq_artifacts+germline_SNV
somatic-hg38_1000g_pon.hg38.vcf.gz: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38
```


### SNP & seq_artifacts+germline_SNV vcf2bed
```bash
cd $SNPDir
alias bcftools="/BioII/lulab_b/baopengfei/anaconda3/bin/bcftools"
alias vcf2bed="/BioII/lulab_b/baopengfei/gitsoft/bedops_v2.4.40/vcf2bed"

#### SNP
#keep snp only
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
bcftools view -v snps ${SNPDir}/${ARTpre}.hg38.vcf.gz | perl -lane 'BEGIN {srand(1984)} if (/^#/) { print } elsif (length($F[3]) == 1) { if (rand(1) > 0.5) {print} }' | bgzip > ${SNPDir}/${ARTpre}_snp.hg38.vcf.gz
gzip -dc ${SNPDir}/${ARTpre}_snp.hg38.vcf.gz > ${SNPDir}/${ARTpre}_snp.hg38.vcf

#vcf2bed, add strand, rm non-auto-chr, rm dup rows
vcf2bed < ${SNPDir}/${ARTpre}_snp.hg38.vcf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' | grep "^chr" | grep -vE "alt|random|decoy|chrUn" | uniq > ${ARTbed}_+
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "-"}' ${ARTbed}_+ > ${ARTbed}_-
cat ${ARTbed}_{+,-} | sort -k1,1 -k2,2n > ${ARTbed}
rm ${ARTbed}_{+,-}
rm ${SNPDir}/${ARTpre}_snp.hg38.vcf
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

### convert bed to gtf
```bash
#gn
#### SNP
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${SNPbed} \
	-o ${SNPbed}.gtf
	
#### seq_artifacts+germline_SNV
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${ARTbed} \
	-o ${ARTbed}.gtf

#tx
#### SNP
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${SNPbed}.tx \
	-o ${SNPbed}.tx.gtf
	
#### seq_artifacts+germline_SNV
python2 /BioII/lulab_b/baopengfei/biosoft/bed2gtf.py \
	-i ${ARTbed}.tx \
	-o ${ARTbed}.tx.gtf
```

### convert gtf to gff.gz
```bash
ulimit -n 4000

#gn
#### SNP
GFFtoTabix.py -i ${SNPbed}.gtf -b 300000 -t $tmpDir > $tmpDir/gn_log1 2>&1 # 33630436
rm ${SNPbed}.gtf 

#### seq_artifacts+germline_SNV
GFFtoTabix.py -i ${ARTbed}.gtf -b 300000 -t $tmpDir > $tmpDir/gn_log2 2>&1 # 2262418
rm ${ARTbed}.gtf 

#### combine
cat ${SNPbed}.sorted.gff > $tmpDir/gn_SNPbed & 
cat ${ARTbed}.sorted.gff > $tmpDir/gn_ARTbed 
cat $tmpDir/gn_ARTbed $tmpDir/gn_SNPbed | sort -k1,1 -k4,4n | uniq | bgzip -c > $SNPDir/SNP_germlineSNV-seqArtifact.sorted.gff.gz # 34840978
tabix  $SNPDir/SNP_germlineSNV-seqArtifact.sorted.gff.gz

#tx
#### SNP
GFFtoTabix.py -b 300000 -i ${SNPbed}.tx.gtf > $tmpDir/log1 2>&1 # 18825510
rm ${SNPbed}.tx.gtf

#### seq_artifacts+germline_SNV
GFFtoTabix.py -b 300000 -i ${ARTbed}.tx.gtf > $tmpDir/log2 2>&1 # 1243549
rm ${ARTbed}.tx.gtf

#### combine
cat ${SNPbed}.tx.sorted.gff > $tmpDir/tx_SNPbed &
cat ${ARTbed}.tx.sorted.gff > $tmpDir/tx_ARTbed
cat $tmpDir/tx_ARTbed $tmpDir/tx_SNPbed | sort -k1,1 -k4,4n | uniq | bgzip -c > $SNPDir/SNP_germlineSNV-seqArtifact.tx.sorted.gff.gz # 19502715
tabix  $SNPDir/SNP_germlineSNV-seqArtifact.tx.sorted.gff.gz 
```




------------------------------------
dsRNAfinder (transcript mode)
------------------------------------
Notes:
	• EM-reassign bowtie2 bam embeded in cfPeak limited to 10-500nt insertion (PE merged into SE as insertion fragment before mapping), further adaptation of PE input and longer insertion fragment might be added 
	• dsRFinder do not need EM if we only focus on finding potential dsRNA regions, quantification of dsRNA could be considered as indepedent process, including re-mapping bam using STAR parameter that suit multi-mapped reads and quantify using tools like Complete-seq or TEtranscript


```bash
sample_id="REDItoolDnaRna_filterSNPinREDItools"
inBam="/BioII/lulab_b/baopengfei/shared_reference/REDItools/testREDItools/merge19_sort_subset01.bam" # 0.01
#inBam="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/tbam/NCpool/bam-EM/merge19_sort/merged.sorted.bam"
outDir="output/GSE71008_19_subset01/${sample_id}"
SNPgff="${SNPbed}.tx.sorted.gff.gz"
ARTgff="${ARTbed}.tx.sorted.gff.gz"
SNPARTgff="${SNPDir}/SNP_germlineSNV-seqArtifact.tx.sorted.gff.gz"
binSize=50
chrSize="${pre}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_sort_uniq_newTxID"
#cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID  | grep -v "^chr" | grep -v "_miR_" | grep -v "^hsa____">  ${chrSize}
binBed="${pre}/exSeek-dev/genome/hg38/tbed/transcriptome_sort_uniq_newTxID.bins.${binSize}.bed6"
gnFa="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa"
EERNum=3 # cutoff of of EER num in $binSize bin
extSize=25 # extended length each side of merged EER bins, usually half of binSize
minLen=50 # min length of extended merged EER bins
maxLen=2000 # max length of extended merged EER bins
mkdir -p $outDir
```

## run EM on bowtie2 bam
run within cfPeak


## run REDItools
```bash
source /BioII/lulab_b/baopengfei/mambaforge/bin/activate REDItools

REDItoolDnaRna.py -t ${cores} \
	-K ${SNPARTgff} \
	-i ${inBam} \
	-f ${gnFa} \
	-o ${outDir} \
	> ${outDir}/log 2>&1 &
#-K $SNPDir/SNP_germlineSNV-seqArtifact.tx.sorted.gff.gz \
#tx: hub 12cores?, 5h, 121354 records
#tx: cnode 6cores?, >24h, 120742 records
#gn: hub 12cores?, 2.5h, 40285769 records
#gn: cnode 6cores, 2.5h, 40070905 records 

#tx:  0.1  subset, cnode 6cores, 2.5h, 48074 records
#tx:  0.01subset, cnode 6cores, 0.5h, 8742 records
#gn: 0.01subset, cnode 6cores, 2min, 96503/97227 records


pid=`head -n3 ${outDir}/log | grep "Analysis ID" | sed s/"Analysis ID: "/""/g`
rediTab="${outDir}/DnaRna_${pid}/outTable_${pid}"
```

## get EER cluster
```bash
bash scripts/getEERcluster.sh \
	$rediTab \
	$chrSize \
	$gnFa \
	$binBed \
	$EERnum \
	$binSize \
	$extSize \
	$minLen \
	$maxLen \
	$tmpDir
```

## get intracellular stable EER cluster
```bash
source /BioII/lulab_b/baopengfei/mambaforge/bin/activate py37
${pre}/exSeek-dev/scripts/rnafold_dinushuffle_parallel.py  \
	${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa 50 1234 ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa.csv ${cores} \
	> ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa.log 2>&1
rm ${rediTab}_filter_sort.bed6.count.merge.ext.filterLen.fa_perm
#tx:  0.01subset, cnode 6cores, 0.1h, 922 records

## convert tx to gn (tx-only)
txBed=${rediTab}_filter_sort.bed6.count.merge.ext.filterLen
gnBed=${rediTab}_filter_sort.bed6.count.merge.ext.filterLen_gn.bed
{{
  grep -v '^chr' ${txBed} | $pre/exSeek-dev/bin/tbed2gbed <(cat  $pre/exSeek-dev/genome/hg38/bed/{long_DNA,long_RNA,tRNA,pri_miRNA,piRNA,rRNA}_newTxID.bed) /dev/stdin /dev/stdout
  awk 'BEGIN{{OFS="\t";FS="\t"}}/^chr/{{print $1,$2,$3,$4,$5,$6,0,0,0,1,$3-$2,0}}' ${txBed}
}} > ${gnBed}
```

## todo: inter-molecular
```bash
source /BioII/lulab_b/baopengfei/mambaforge/bin/activate IntaRNA

IntaRNA --threads 4 -q query.fa -t target.fa > outpre


```




------------------------------------
dsRNAfinder (genome mode)
------------------------------------
ssh hub

```bash
sample_id="REDItoolDnaRna_filterSNPinREDItools"
inBam="/BioII/lulab_b/baopengfei/shared_reference/REDItools/testREDItools/hg38_v38_sortbyCoord_subset01.bam" # 0.01
#inBam="/BioII/lulab_b/baopengfei/tmp/20230523-S048-1_BC11/hg38_v38_sortbyCoord.bam"
outDir="output/SLE_subset01/${sample_id}"
SNPgff="${SNPbed}.sorted.gff.gz"
ARTgff="${ARTbed}.sorted.gff.gz"
SNPARTgff="${SNPDir}/SNP_germlineSNV-seqArtifact.sorted.gff.gz"
binSize=50
chrSize="${pre}/../exOmics/DNA-seq/genome/chrom.size"
binBed="${pre}/../exOmics/DNA-seq/genome/hg38.bins.${binSize}.bed"
gnFa="/BioII/lulab_b/baopengfei/shared_reference/hg38/genome.fa"
EERNum=3 # cutoff of of EER num in $binSize bin
extSize=25 # extended length each side of merged EER bins
minLen=50 # min length of extended merged EER bins
maxLen=2000 # max length of extended merged EER bins
```

## run EM on STAR bam (holding)
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
same module of tx


## get EER cluster
same module of tx

## get intracellular stable EER cluster
same module of tx


## todo: annotation (gn mode only)
```bash
/BioII/lulab_b/baopengfei/projects/motif-RBP-EDA/scripts/lulab/sequential.assign.long.sh
```
