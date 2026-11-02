#!/usr/bin/env bash
# Snippy example: E. faecalis isolates from host 22M21
# Author: Zhu Zhihao  Date: 2026.10.30
# Necessary files: metadata.txt (Prepare your own metadata. The first, fourth, and tenth columns are mandatory); ref/refseq_reference_genomes.tsv; script/plot_heatmap.R and script/plot_tree.R
# FASTQ: download, or place seq/{id}_1.fq.gz and seq/{id}_2.fq.gz yourself

mkdir -p snippy && cd snippy
# Move necessary files to ~/snippy
mkdir -p seq temp result
# The content of the directory should be: metadata.txt 

# 0. Environment and software
# conda create -y -n snippy -c bioconda -c conda-forge snippy=4.6.0 fastp=1.4.0 snp-dists=0.8.2 fasttree=2.2.0 minimap2=2.28 bedtools=2.31.1 samtools=1.20 phastaf=0.1.0 isescan=1.7.3
# gubbins needs python<=3.10, use a separate env: conda create -y -n gubbins -c bioconda -c conda-forge gubbins=3.4.3 python=3.10
# conda activate snippy

# Download GSA FASTQ (Run=CRR* in metadata.txt; rename to seq/{id}_1/2.fq.gz)
# conda install -y -c bioconda -c conda-forge iseq
# cut -f1,10 metadata.txt | tail -n +2 | while read id run; do
#   if [ -s seq/${id}_1.fq.gz ] && [ -s seq/${id}_2.fq.gz ]; then continue; fi
#   mkdir -p seq/tmp_$run
#   (cd seq/tmp_$run && iseq -i "$run" -g -a)
#   mv seq/tmp_$run/*_f1.fastq.gz seq/${id}_1.fq.gz
#   mv seq/tmp_$run/*_f2.fastq.gz seq/${id}_2.fq.gz
#   rm -rf seq/tmp_$run
# done

# 1. QC
cut -f1 metadata.txt | tail -n +2 | while read id; do
  mkdir -p temp/fastp/$id
  fastp -w 64 \
    -i seq/${id}_1.fq.gz -I seq/${id}_2.fq.gz \
    -o temp/fastp/$id/${id}_1.fq.gz -O temp/fastp/$id/${id}_2.fq.gz \
    -j temp/fastp/$id/${id}.json -h temp/fastp/$id/${id}.html
done

# 2. reference genome
# Use Species from metadata.txt to find the RefSeq URL in ref/refseq_reference_genomes.tsv, then download as ref/ref.fna
awk -F'\t' 'NR==FNR{if(NR==2)sp=$4;next} NR>1 && $6~"^"sp"($| )"{print $12;exit}' \
  metadata.txt ref/refseq_reference_genomes.tsv > ref/ref.url
wget -c -i ref/ref.url -O ref/ref.fna.gz
gunzip -c ref/ref.fna.gz > ref/ref.fna
samtools faidx ref/ref.fna

# 3. static mask (repeats + phage + IS)
mkdir -p temp/mask
# repeats
minimap2 -D -P -w19 -m200 ref/ref.fna ref/ref.fna \
  | cut -f1,3,4 | sort -k1,1 -k2,2n \
  | bedtools slop -i stdin -g ref/ref.fna.fai -r 1 -l 0 \
  | bedtools merge -i stdin > temp/mask/repeats.bed
# phage
phastaf --outdir temp/mask/phastaf_out ref/ref.fna
cp temp/mask/phastaf_out/phage.bed temp/mask/phage.bed
# is
isescan.py --seqfile ref/ref.fna --output temp/mask/isescan_out --nthread 8
awk 'BEGIN{OFS="\t"} /^[^#]/ && NF>=5 {s=$4-1; if(s<0)s=0; print $1,s,$5}' \
  temp/mask/isescan_out/*/*.gff \
  | sort -k1,1 -k2,2n | bedtools merge > temp/mask/is.bed
# merge
cat temp/mask/repeats.bed temp/mask/phage.bed temp/mask/is.bed \
  | awk 'NF>=3 && $1!~/^#/ {print $1"\t"$2"\t"$3}' \
  | sort -k1,1 -k2,2n | bedtools merge \
  > temp/mask/mask_static.bed

# 4. Snippy
# Prefer QC-ed paired-end FASTQ (--R1/--R2). --ctgs accepts assemblies but chops them into ~250 bp pseudo-reads and drops depth/quality.
# snippy --cpus 64 --outdir result/snippy/$id --ref ref/ref.fna --ctgs seq/${id}.fasta --mincov 10 --minfrac 0.9 --force
cut -f1 metadata.txt | tail -n +2 | while read id; do
  snippy --cpus 64 --outdir result/snippy/$id \
    --ref ref/ref.fna \
    --R1 temp/fastp/$id/${id}_1.fq.gz --R2 temp/fastp/$id/${id}_2.fq.gz \
    --mincov 10 --minfrac 0.9 --force
done

# 5. filter samples by ALIGNED/LENGTH (>=0.8); then link pass dirs
# Low ALIGNED/LENGTH shrinks the core alignment and can collapse pairwise SNP distances.
# Adjust 0.8, or edit result/core.highAlign.txt and re-run the snippy-core --mask step.
mkdir -p temp/core_all result/core
snippy-core --ref ref/ref.fna --prefix temp/core_all/core_all result/snippy/*
awk 'BEGIN{OFS=FS="\t"} $1=="ID"{print;next} $1!="Reference" && $2+0>0 && ($3+0)/($2+0)>=0.8{print}' \
  temp/core_all/core_all.txt | sort -k3,3n > result/core.highAlign.txt
mkdir -p result/snippy_pass
tail -n +2 result/core.highAlign.txt | cut -f1 | while read id; do
  ln -sfn "$(pwd)/result/snippy/$id" result/snippy_pass/$id
done

# 6. core SNPs with static mask (high-align samples only)
# Mask BED must use the same contig names as the reference.
# awk 'NR==FNR{r[$1];next} !($1 in r){print; exit 1}' ref/ref.fna.fai temp/mask/mask_static.bed
# Species-specific (repeats/phage/IS)
snippy-core --ref ref/ref.fna \
  --mask temp/mask/mask_static.bed \
  --prefix result/core/core result/snippy_pass/*

# 7. tree and pairwise distances
FastTree -gtr -nt result/core/core.aln > result/core/core.tree
snp-dists result/core/core.aln > result/core/core.snpdists.tsv
# heatmap (needs R: ggplot2, scales)
Rscript ../script/plot_heatmap.R result/core/core.snpdists.tsv figure
# tree (needs R: ggplot2, ggtree, treeio, ape)
Rscript ../script/plot_tree.R result/core/core.tree result/core/core.snpdists.tsv figure

# 8. optional: Gubbins then tree/distances
# mkdir -p result/gubbins
# snippy-clean_full_aln result/core/core.full.aln > result/gubbins/core.clean.full.aln
# conda activate gubbins
# run_gubbins.py -p result/gubbins/gubbins result/gubbins/core.clean.full.aln
# conda activate snippy
# snp-sites -c result/gubbins/gubbins.filtered_polymorphic_sites.fasta > result/gubbins/core_gubbins.aln
# FastTree -gtr -nt result/gubbins/core_gubbins.aln > result/gubbins/core_gubbins.tree
# snp-dists result/gubbins/core_gubbins.aln > result/gubbins/core_gubbins.snpdists.tsv
# Rscript ../script/plot_heatmap.R result/gubbins/core_gubbins.snpdists.tsv figure/gubbins
# Rscript ../script/plot_tree.R result/gubbins/core_gubbins.tree result/gubbins/core_gubbins.snpdists.tsv figure/gubbins
