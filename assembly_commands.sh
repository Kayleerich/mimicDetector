#!/bin/bash
strain=$1
fastq1='data/raw_reads/'${strain}'_L001_R1_001.fastq.gz'
fastq2='data/raw_reads/'${strain}'_L001_R2_001.fastq.gz'
trim1='data/raw_reads/'${strain}'_R1.trimmed.fastq.gz'
trim2='data/raw_reads/'${strain}'_R2.trimmed.fastq.gz'
dir_a='results/'${strain}'/assembly'

## conda environment: qc.yml
fastp -i ${fastq1} -o ${trim1} -I ${fastq2} -O ${trim1} -j results/{strain}/qc/${strain}_fastp.json -h results/${strain}/qc/{strain}_fastp.html
unicycler -1 ${trim1} -2 ${trim2} -o ${dir_a} -t 8
quast ${dir_a}/assembly.fasta -o ${dir_a}/quast -1 ${trim1} -2 ${trim2} -t 8
## conda environment: busco.yml
busco -f -i ${dir_a}/assembly.fasta -o ${strain}_unicycler_busco --out_path results/${strain}/busco -c 8 -m geno -l bacillales_odb10
## conda environment: bakta.yml
bakta --force --db bakta_downloads/db -o results/${strain}/bakta -p ${strain}_bakta -t 16 --genus Staphylococcus --species xylosus --gram + ${dir_a}/assembly.fasta

## to run:
# for file in data/raw_reads/*_L001_R1_001.fastq.gz; do
#     strain=$( echo $file | cut -f 3 -d '/' | cut -f 1 -d '.' )
#     sh assembly_commands.sh ${strain}
# done