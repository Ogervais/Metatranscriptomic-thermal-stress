# Metatranscriptomic- red coral thermal stress 

Pipeline used for the analysis in Functional stability of Spirochaetota symbionts in the precious octocoral Corallium rubrum under heat stress paper (link)

## Table of content
* Description
* Software used
* Workflow diagram
* 1 Preprocess
* 2 Taxonomic analysis
* 3 Metatranscriptome assembly
* 4 Fonctional annotation

## Description

The input for this analysis was a paired-end reads of 150bp sequenced on Novaseq6000 from red coral RNA. 
Analyses were made on conda environment.

## Software used

* fastqc/v.0.11.9
* multiqc/v.1.13
* trimmomatic/v.0.39
* sortmerna/v.4.3.6
* bwa-mem2/v.2.2.1
* samtools/v1.6
* RNAspades/v.3.15.5
* cdhhits/v.4.8.1
* amos/v.3.1.0
* QUAST/v.5.0.2
* eggNOG-mapper/v2.1.11
* prodigal/v.2.6.3
* mmseqs2/v.13.45111
* diamond/v.2.1.9
* bowtie2/v.2.2.5
* subread/v.2.0.1

## Workflow diagram


## 1 Preprocess

### 1.1 Quality control

In a first step the quality of our reads were checked using `fastqc` (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 
To run fastQC in parallel on each fastq file `-t` was set to the same number of fastq file.

```
fastqc ./fastq/*.fastq \
    -o ./fastqc/output \
    -t 42
```

To inspect more easily our fastQC ouput we have run `multiqc` (https://multiqc.info/) to collate all the FastQC reports into one report.

```
multiqc \
    -o ./multiqc \
    ./fastqc/*_fastqc.zip
```

### 1.2 Trimming
For this step, various action were realized using `Trimmomatic` (https://github.com/usadellab/Trimmomatic):
- sequence adaptor were removed using TrueSeq3-PE 
- low quality reads were trimmed 
- reads below a minimum length of 36bp were removed
- duplicate reads were removed

Our script was made in order to have all our sample process without having to write the name of each sample.

```
OUTPUT=./trimmed

for f in *_R1.fastq.gz

do
n=${f%%_R1.fastq.gz}
trimmomatic PE -threads 40 \
${n}_R1.fastq.gz ${n}_R2.fastq.gz \
$OUTPUT/${n}_R1_trimmed.fastq.gz $OUTPUT/${n}_R1_unpaired_trimmed.fastq.gz \
$OUTPUT/${n}_R2_trimmed.fastq.gz $OUTPUT/${n}_R2_unpaired_trimmed.fastq.gz \
ILLUMINACLIP:./TrueSeq3-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

done
```

### 1.3 Removal of host RNA contamination

The metatranscriptomic data were extracted from a host (coral tissue) so in order to remove any RNA belonging to the host and food we use the transcriptome of the red coral (https://doi.org/10.21203/rs.3.rs-4582739/v1), mussel (GCA_025215535.1) and artemia (GCA_019857095.1) with `bwa` (https://github.com/bwa-mem2/bwa-mem2). Alignment result were then tranform in bam file using `samtools` (https://www.htslib.org/). We only kept target reads that was not aligned with the host transcriptome. 

```
#prepare database index
bwa-mem2 index -p <prefix> <in.fasta>
```

```
INPUT=./sortmerna
INPUT2=./bwa
DATABASE=./database/
OUTPUT=./bwa

#mapping using bwa-mem2

for f in *_R1_trimmed.fastq.gz

do 

n=${f%%_R1_trimmed.fastq.gz}

bwa-mem2 mem \
    -t 100 \
    ${DATABASE}/<db_prefix.fasta> \
    ${INPUT}/{n}_R1_trimmed.fastq.gz \
    ${INPUT}/{n}_R2_trimmed.fastq.gz > ${OUTPUT}/${n}_genome.sam

done
#tansform sam to bam file

for f in *_genome.sam

do

n=${f%%_genome.sam}

samtools view -@ 20 -bS ${INPUT2}/${n}_genome.sam > ${OUTPUT}/${n}_genome.bam
done

#extract unmapped reads

for f in *_genome.bam

do

n=${f%%_genome.bam}

samtools view -@ 20 -b -f 12 -F 256 \
    ${n}_genome.bam > ${n}_unmapped.bam
done

#sorted hits

for f in *_unmapped.bam

do

n=${f%%_unmapped.bam}

samtools sort -n -m 5G -@ 20 ${n}_unmapped.bam \
        -o ${n}_unmapped_sorted.bam
done

#split paired-end reads in fastq file
for f in *_unmapped_sorted.bam

do

n=${f%%_unmapped_sorted.bam}

samtools fastq -@ 20 ${n}_unmapped_sorted.bam \
        -1 ${n}_host_removed_R1.fastq.gz \
        -2 ${n}_host_removed_R2.fastq.gz \
        -0 -s -n
done

```


### 1.4 Removal of rRNA sequences

rRNA genes were partially remove during the wetlab protocol but because they are highly expressed it's still important to run a bioinformatic rRNA removal. SortmeRNA was used to remove rRNA (https://github.com/biocore/sortmerna) using various database you can download here: https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz

```
OUTPUT=./sortmerna
DATABASE=./database/rRNA/

for f in *_host_removed_R1.fastq.gz

do
n=${f%%_host_removed_R1.fastq.gz}
sortmerna \
    -ref ${DATABASE}/smr_v4.3_default_db.fasta \
    -ref ${DATABASE}/smr_v4.3_fast_db.fasta \
    -ref ${DATABASE}/smr_v4.3_sensitive_db.fasta \
    -ref ${DATABASE}/smr_v4.3_sensitive_db_rfam_seeds.fasta \
    -reads ${n}_host_removed_R1.fastq.gz\ 
    -reads ${n}_host_removed_R2.fastq.gz \ 
    -fastx \ 
    -paired_in \ 
    -out2 \ 
    -aligned ${OUTPUT}/${n}_rRNA_reads \ 
    -other ${OUTPUT}/${n}_non_rRNA_reads \ 
    -workdir ${HOME}/ \
    -threads 60 
done
```

## 2 Taxonomy analysis  

The taxonomic classification was made using `kaiju` (https://bioinformatics-centre.github.io/kaiju/) and nr database (https://bioinformatics-centre.github.io/kaiju/downloads.html). 

```
DATABASE=./database/kaiju/nr
INPUT=./sortmerna/result
OUTPUT=./taxonomy/kaiju
THREADS=22

##create database index

kaiju-makedb -t 50 -s nr

##run kaiju

for f in RNA-*_non_rRNA_reads_fwd.fq.gz #for each sample

do

n=${f%%_non_rRNA_reads.fwd.fq.gz} # strip part of file name

kaiju-multi \
    -z ${THREADS} \
    -t ${DATABASE}/nodes.dmp \
    -f ${DATABASE}/kaiju_db_nr.fmi \
    -i ${INPUT}/${n}_non_rRNA_reads_fwd.fq.gz \
    -j ${INPUT}/${n}_non_rRNA_reads_rev.fq.gz \
    -v \ #verbose
    -e 5 \
    -E 0.05 \
    -o ${OUTPUT}/${n}_kaiju.out

done 
```
In order to analyse the taxonomy of our sample a table was created by transforming result of kaiju. 

```
DATABASE=./database/kaiju/nr
INPUT=./taxonomy/kaiju
OUTPUT=./taxonomy/kaiju

#create table of taxonomy

for f in RNA-*_kaiju.out

do

n=${f%%_kaiju.out}

./kaiju/bin/kaiju2table \
        -t ${DATABASE}/nodes.dmp \
        -n ${DATABASE}/names.dmp \
        -r phylum \ #can be change by order, familly, genus, ...
        -p \
        -u \
        -o ${OUTPUT}/${n}_kaiju_phylum_summary.tsv \
        ${INPUT}/${n}_kaiju.out

./kaiju/bin/kaiju-addTaxonNames \
        -t ${DATABASE}/nodes.dmp \
        -n ${DATABASE}/names.dmp \
        -p \
        -i ${INPUT}/${n}_kaiju.out \
        -o ${OUTPUT}/${n}_kaiju_summary.tsv \

done
```

The result were converted to be compatible with Krona and be able to visualize the result.

```
DATABASE=./database/kaiju/nr
INPUT=./taxonomy/kaiju
OUTPUT=./taxonomy/kaiju
THREADS=22

##transform result into a readable file for krona

for f in RNA-*_kaiju.out

do

n=${f%%_kaiju.out}

kaiju2krona \
    -t ${DATABASE}/nodes.dmp \
    -n ${DATABASE}/names.dmp \
    -i ${OUTPUT}/${n}_kaiju.out \
    -o ${OUTPUT}/${n}_kaiju.out.krona

#create html file of krona

for f in RNA-*_kaiju.out.krona

do

n=${f%%_kaiju.out.krona}

ktImportText -o ${n}_kaiju.out.html ${n}_kaiju.out.krona

done
```


## 3 Metatranscriptome assembly

### 3.1 Assemble target reads

We have now only mRNA reads that will be assemble in contigs. All our sample were assembled using `RNAspades` (https://cab.spbu.ru/software/rnaspades/) and verified using `quast`.

```
INPUT=./sortmerna/result
OUTPUT=./assembly/spades
THREADS=40

spades.py --rna\
    --pe1-1 ${INPUT}/RNA-*_fwd.fq.gz \
    -pe1-2 ${INPUT}/RNA-*_rev.fq.gz \
    -t ${THREADS} \
    -o ${OUTPUT}/RNA-*_assembly
```

A mega assembly was created by merging all the assembly, detect and merge contig present more than one time
using `minimus2` () to assemble overlapping contigs. -> can be use to analyse individual transcript against a "database"
`cd-hit` () was used to remove contig redundancy before to merged contigs.

```
INPUT=./assembly/spades/RNA-*
MERGED1=${INPUT}/merged.fasta
MERGED2=${INPUT}/merged_cdhit.fasta
MERGED3=${INPUT}/merged_amos.afg
AFG=${INPUT}/merged_amos
THREADS=40

#merge the assembly in one file

cat ${INPUT}/transcripts.fasta > ${MERGED1}

#use cd-hits to remove redundant contigs 

cd-hit-est \
    -i ${MERGED1} \
    -o ${MERGED2} \
    -T ${THREADS} \
    -c 0.99 \ #identity threshold
    -d 100 \
    -aS 0.9 \ #alignment coverage for the shorter sequence if set to 0.9, the alignment must covers 90% of the sequence 
    -M 1400 #or more

#transform fasta file to afg file for AMOS

toAmos -o ${MERGED3} -s ${MERGED2}

#minimus2: use to assemble overlapping contigs (AMOS)

minimus2 ${AFG} \
    -D OVERLAP=100 \ #minimap overlap
    -D MINID=95 \ #minimum overlap %id for align
    -D THREADS=${THREADS}
```
### 3.2 Evaluate the quality of the metatranscriptome

Each metatranscriptome was then evaluates using `quast` (https://quast.sourceforge.net/) against our mega assembly.

```
INPUT=./assembly/spades
OUTPUT=./assembly/quast

quast ${INPUT}/RNA-*_assembly/transcripts.fasta ${INPUT}/RNA-*_assembly/transcripts.fasta \ #one for each assembly
    -r ${INPUT}/merged_amos.fasta \
    -o ${OUTPUT}/report
```

## 4. Fonctional annotation
### 4.1 Identification of open reading frame (ORFs) and annotation

`EggNOG-mapper` (https://github.com/eggnogdb/eggnog-mapper) was used for functional annotation of sequences. The parameter used were `prodigal` for gene prediction and `mmseqs2` for protein sequence searching. This annotation allow us to get COG, KEGG, GO, EC and the closest taxon using orthology.

```
#install database 

download_eggnog_data.py \ 
    --data_dir ./database/eggnog \
    -F
    -P
```

```
#create mmseqs database only for bacteria

create_dbs.py \
    -m mmseqs \
    --dbname bacteria \
    --taxa Bacteria
```

```
INPUT=./assembly/
OUTPUT=./assembly/eggnog
DATABASE=./database/eggnog

#gene prediction and functional annotation

emapper.py \
    --cpu 30 \
    -m  mmseqs \
    --itype metagenome \
    -i ${INPUT}/merged_amos.fasta \
    --data_dir ${DATABASE} \
    --mmseqs_db ${DATABASE}/bacteria.mmseqs/bacteria.mmseqs \
    -o ${OUTPUT}/eggnog_mmseqs \
    --genepred prodigal \
    --decorate_gff yes \
    --pfam_realign realign \
    --override

```

### 4.2 Annoation using virulence database

`diamond` (https://github.com/bbuchfink/diamond) was used to performed the annotation using two database VFDB (http://www.mgc.ac.cn/VFs/main.htm) and AMR (https://card.mcmaster.ca/download) on the fasta file from EggNOG-mapper providing ORFs. The same script was used for both database just by changing the database path and name.

```
DATABASE=./database/VFDB
INPUT=./assembly/eggnog
OUTPUT=./annotation/VFDB
THREADS=50

#create database index
diamond makedb \
    --in ${DATABASE}/VFBD_setB_pro.fas.gz \
    -d ${DATABASE}/VFDB_setB \
    --threads ${THREADS}

#annotate sample 

diamond blasp \
    -q ${INPUT}/eggnog_mmseqs.emapper.genepred.fasta \
    -d ${DATABASE}/VFDB_setB.dmnd \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore full_qseq full_sseq \
    -o ${OUTPUT}/VFG.tsv \
    --ultra-sensitive \
    -k 1 \
    --evalue 0.00001 \
    --log \
    --threads {THREADS}
```


### 4.3 Align target reads to assemblies

Each assembly were then align against the mega assembly using `bowtie2` (https://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

```
DATABASE=./assembly/spades
INPUT=./sortmerna/result
OUTPUT=./gene_expression

#build database index

bowtie2-build ${DATABASE}/merged_amos.fasta ${DATABASE}/mega_assembly

#align target to assembly

bowtie2 --very-sensitive -x ${DATABASE}/mega_assembly --threads 32 \
    -1 ${INPUT}/RNA-*_non_rRNA_reads_fwd.fq.gz \
    -2 ${INPUT}/RNA-*_non_rRNA_reads_rev.fq.gz \
    -S ${OUTPUT}/${n}.sam

#transform data in bam file

for f in *.sam 

do

n=${f%%.sam}

samtools view -@ 20 -bS ${n}.sam > ${n}.bam

done

#sorted hits

for f in *.bam 

do

n=${f%%.bam}

samtools sort -n -m 5G -@ 20 ${n}.bam \
    -o ${n}_sorted.bam

done 
```

`featurecounts` (https://subread.sourceforge.net/featureCounts.html) was then used to summarized reads count in order to analyse RNA expression.

```
GFF=./eggnog
INPUT=./gene_expression
OUTPUT=./gene_expression/featurecounts

featureCounts \
    -p \
    -O \
    -T 10 \
    -t 'CDS' \
    -a ${GFF}/eggnog_mmseqs.emapper.genepred.gff \
    -o ${OUTPUT}/featurecounts_eggnog.txt \
    -F 'GFF3' \
    -g ID \
    ${INPUT}/*_sorted.bam

```

## 5 Differential expression analysis



