#! /usr/bin/env bash

# Tiana Stastny, MOLB 7621, HW2

# 1. Identify size of largest overlap between CTCF and H3K4me3 locations.

datasets=$HOME/MOLB7621/data-sets

tfbs_bed=$datasets/bed/encode.tfbs.chr22.bed.gz
histone_bed=$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz

zcat $tfbs_bed | awk '$4 == "CTCF"' > ctcf-peaks.bed

answer_1=$(bedtools intersect -a ctcf-peaks.bed -b $histone_bed -wo \
	| awk '{print $NF}' \
	| sort -nr | head -n1)

echo "answer_1: $answer_1" > answers.yml


# 2. Calculate the GC content of nucleotides 19,000,000 to 19,000,500 on chr22 of
# hg19 genome build.  Report as fraction.

hg19chr22=$datasets/fasta/hg19.chr22.fa

echo -e "chr22\t19000000\t19000500" > tmp.bed

answer_2=$(bedtools nuc -fi $hg19chr22 -bed tmp.bed \
	| awk '(NR>1) {print $5}')

echo "answer_2: $answer_2" >> answers.yml

# 3. Find length of CTCF ChIP-seq peak that has largest mean signal in ctcf.hela.chr22.bg.gz

ctcfsignal=$datasets/bedtools/ctcf.hela.chr22.bg.gz

answer_3=$(bedtools map -o mean -a ctcf-peaks.bed -b $ctcfsignal -c 4 \
	| awk '{print $5, $3-$2}' \
	| sort -n \
	| tail -n1 \
	| awk '{print $2}')

echo "answer_3: $answer_3" >> answers.yml

# 4. Find gene promoter with highest median signal in ctcf.hela.chr22.bg.gz. Report gene name.

tsschr22=$datasets/bed/tss.hg19.chr22.bed.gz
genome=$datasets/bedtools/hg19.genome

bedtools flank -i $tsschr22 -g $genome -r 0 -l 1000 -s > tmp2.bed
sort -k2,2n tmp2.bed > tmpsort.bed

answer_4=$(bedtools map -a tmpsort.bed -b $ctcfsignal -c 4 -o median  \
	| awk '{print $4, $5}' \
	| sort -k2n \
	| tail -n1 \
	| awk '{print $1}')

echo "answer_4: $answer_4" >> answers.yml


#5. Find longest interval on chr22 that is not covered by genes.hg19.bed.gz.

geneint=$datasets/bed/genes.hg19.bed.gz

answer_5=$(bedtools complement -i $geneint -g $genome \
	| grep chr22 \
	| awk '{print $1, $2, $3, $3-$2}' \
	| sort -k4n \
	| tail -n1 \
	| awk '{print $1":"$2"-"$3}')

echo "answer_5: $answer_5" >> answers.yml










