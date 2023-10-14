------------------------------------
# prepare environment  & reference
------------------------------------
## prepare&activate conda environment 
ssh baopengfei@101.6.121.24 # hub-remote
```bash
cd /BioII/lulab_b/baopengfei/shared_reference/REDItools/testREDItools

### REDItools env
mamba env create --name REDItools --file /BioII/lulab_b/baopengfei/shared_reference/REDItools/mamba_env.yml

### RNAfold env
mamba create -n py37 python=3.7 seaborn multiprocess
```

## prepare filtering reference
### download reference
#### SNP
All_20180418.vcf.gz: https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/
#### seq_artifacts+germline_SNV
somatic-hg38_1000g_pon.hg38.vcf.gz: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38

### SNP & seq_artifacts+germline_SNV vcf2bed (only needed for 1st run)
```bash
SNPDir="/BioII/lulab_b/baopengfei/projects/exOmics/DNA-seq/ref/SNV_ref"
alias bcftools="/BioII/lulab_b/baopengfei/anaconda3/bin/bcftools"
alias vcf2bed="/BioII/lulab_b/baopengfei/gitsoft/bedops_v2.4.40/vcf2bed"

#### SNP
#keep snp only
SNPpre="All_20180418"
SNPbed="${SNPDir}/${SNPpre}_snp_sort_addStrandfilterChr.hg38.bed"
bcftools view -v snps ${SNPDir}/${SNPpre}.vcf.gz | perl -lane 'BEGIN {srand(1984)} if (/^#/) { print } elsif (length($F[3]) == 1) { if (rand(1) > 0.5) {print} }' | bgzip > ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz
#vcf2bed 
gzip -dc ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz > ${SNPDir}/${SNPpre}_snp.hg38.vcf
vcf2bed < ${SNPDir}/${SNPpre}_snp.hg38.vcf > ${SNPDir}/${SNPpre}_snp_sort.hg38.bed
rm ${SNPDir}/${SNPpre}_snp.hg38.vcf
#add strand, rm non-auto-chr
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' ${SNPDir}/${SNPpre}_snp_sort.hg38.bed | grep "^chr" | grep -vE "alt|random|decoy|chrUn"  > ${SNPbed} 

#### seq_artifacts+germline_SNV
#keep snp only
SNPpre="somatic-hg38_1000g_pon"
ARTbed="${SNPDir}/${SNPpre}_snp_sort_addStrandfilterChr.hg38.bed"
bcftools view -v snps ${SNPDir}/${SNPpre}.hg38.vcf.gz | perl -lane 'BEGIN {srand(1984)} if (/^#/) { print } elsif (length($F[3]) == 1) { if (rand(1) > 0.5) {print} }' | bgzip > ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz
#vcf2bed 
gzip -dc ${SNPDir}/${SNPpre}_snp.hg38.vcf.gz > ${SNPDir}/${SNPpre}_snp.hg38.vcf
vcf2bed < ${SNPDir}/${SNPpre}_snp.hg38.vcf > ${SNPDir}/${SNPpre}_snp_sort.hg38.bed
rm ${SNPDir}/${SNPpre}_snp.hg38.vcf
#add strand, rm non-auto-chr
awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" "+"}' ${SNPDir}/${SNPpre}_snp_sort.hg38.bed | grep "^chr" | grep -vE "alt|random|decoy|chrUn"  > ${ARTbed} 
```


### convert bed to gff/gtf (only needed for 1st run, gn-mode-only)
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

### convert gtf to gff.gz (only needed for 1st run, gn-mode-only)
```bash
#### SNP
GFFtoTabix.py -i ${SNPbed}.gtf -u

#### seq_artifacts+germline_SNV
GFFtoTabix.py -i ${ARTbed}.gtf -u
```


### convert gn coordinate to tx coordinate (only needed for 1st run, tx-mode-only)
```bash
#### SNP
/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores ${cores} \
	-i ${SNPbed} \
	-o ${SNPbed}.tx
#1/50 converted, need add 4DNA gn2tx

#### seq_artifacts+germline_SNV
/usr/bin/Rscript \
/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/scripts/gbed2tbed.R \
	--cores ${cores} \
	-i ${ARTbed} \
	-o ${ARTbed}.tx
```

### convert bed to gff/gtf (only needed for 1st run, tx-mode-only)
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

### convert gtf to gff.gz (only needed for 1st run, tx-mode-only)
```bash
#### SNP
GFFtoTabix.py -i ${SNPbed}.tx.gtf -u

#### seq_artifacts+germline_SNV
GFFtoTabix.py -i ${ARTbed}.tx.gtf -u
```




------------------------------------
dsRNAfinder (genome mode)
------------------------------------

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
### make tx bin window (only needed for 1st run)
binSize=50
EERNum=3
pre="/BioII/lulab_b/baopengfei/projects/WCHSU-FTC"
chrSize="${pre}/../exOmics/DNA-seq/genome/chrom.size"
binBed="${pre}/hg38.bins.${binSize}.bed"
#make bins bed
bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}.tmp
#add strand
sed s/"+"/"-"/g ${binBed}.tmp >  ${binBed}.tmp2
cat ${binBed}.tmp ${binBed}.tmp2 | sort -k1,1 -k2,2n > ${binBed}

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
#
```

## todo: inter-molecular
#mamba create -n IntaRNA -c conda-forge -c bioconda intarna




------------------------------------
dsRNAfinder (transcript mode)
------------------------------------

## run REDItools
```bash
source activate REDItools

cores=10
tmpDir='tmp'
sample_id='test'
#inBam='/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/output/GSE71008_NCpool/tbam/NCpool/bam-EM/merge19_sort/merged.sorted.bam'
inBam='./testREDItools/tx_19.bam'
gnFa='/BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/fasta_newTxID/combine19.fa'
#outDir='output/GSE71008_19/REDItoolDnaRna'
outDir='output/GSE71008_19/test'

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
### make tx bin window (only needed for 1st run)
binSize=50
EERNum=3
pre='/BioII/lulab_b/baopengfei/projects/WCHSU-FTC'
chrSize='${pre}/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_sort_uniq_newTxID'
binBed='${pre}/exSeek-dev/genome/hg38/tbed/transcriptome_sort_uniq_newTxID.bins.${binSize}.bed6'

#make bins bed
cat /BioII/lulab_b/baopengfei/projects/WCHSU-FTC/exSeek-dev/genome/hg38/chrom_sizes/transcriptome_genome_sort_uniq_newTxID  | grep -v "^chr" | grep -v "_miR_" | grep -v "^hsa____">  ${chrSize}
bedtools makewindows -g ${chrSize} -w ${binSize} | awk '{print $1 "\t" $2 "\t" $3 "\t" "." "\t" "." "\t" "+"}'  > ${binBed}

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


bed12-->gtf: https://github.com/alejandrogzi/bed2gtf
tx <--> gn: https://github.com/mt1022/gppy (allow spliced regions)

(## optional: filter tRNA,miRNA?)
```

## todo: inter-molecular
#mamba create -n IntaRNA -c conda-forge -c bioconda intarna

