#!/bin/bash

# Script variables: 
# ${DB} = Database prefix 
# ${QUERY} = Query file prefix 
# ${s} = Subsampling seed 
# ${m} = Motif file identifier, Glam2scan only 
# ${NAME} = Program and file identifier for output files  
# ${KSIZE} = k-mer size (range from 5 to 14) 
# ${MATRIX} = PAM30 or BLOSUM62 
# ${WS} = Wordsize, BLASTP only (2 or 3) 

# BLASTP: 
function blastp_benchmark { 
/usr/bin/time -o $3_$4mers_ws$6_$5_times.txt -a -f ‘real %E\\tuser %U\\tsys %S\\tmem %M’ blastp -query subsamples/$2.fasta -db $1 -outfmt 6 -out $3_$4mers_$2_ws$6_$5.out -evalue 10 -max_hsps 200 -word_size $6 -ungapped -matrix $5 -comp-based-stats 0 -seg no; 
} 
blastp_benchmark ${DB} ${QUERY}${s} ${NAME} ${KSIZE} ${MATRIX} ${WS} 

# DIAMOND*:  
function diamond_benchmark { 
/usr/bin/time -o $3_$4mers_$5_times.txt -a -f ‘real %E\\tuser %U\\tsys %S\\tmem %M’ diamond blastp -q subsamples/$2.fasta -d $1 -p 1 –more-sensitive -o $3_$4mers_$2_$5.out –comp-based-stats 0 –matrix $5 –masking 0 –evalue 10 –max-hsps 200; 
} 
diamond_benchmark ${DB} ${QUERY}${s} ${NAME} ${KSIZE} ${MATRIX} 
# *We replaced the --–more-sensitive option with --–ultra-sensitive for Diamond2 

# TOPAZ: 
function topaz_benchmark { 
/usr/bin/time -o $3_$4mers_ws3_$5_times.txt -a -f ‘real %E\\tuser %U\\tsys %S\\tmem %M’ /home/topaz/src/topaz search -T 1 -p $1 -E 10 -l 3 -B -o 100 -M -m $5 -H 200 -f subsamples/$2.fasta > $3_$4mers_$2_ws3_$5.out; 
} 
cat ${DB}.* > /dev/null; topaz_benchmark ${DB} ${QUERY}${s} ${NAME} ${KSIZE} ${MATRIX}  

# PHMMER: 
function phmmer_benchmark { 
/usr/bin/time -o $3_$4mers_$5_times.txt -a -f 'real %E\\tuser %U\\tsys %S\\tmem %M' phmmer --domtblout $3_$4mers_$2_$5.out --mx $5 --cpu 1 subsamples/$2.fasta $1 | sed -n '/>>/,/Inter/p' | grep -v 'Internal pipeline statistics summary:' > $3_$4mers_$2_$5_alns.out; 
} 
phmmer_benchmark ${DB} ${QUERY}${s} ${NAME} ${KSIZE} ${MATRIX} 

# GLAM2SCAN: 
function glam2_benchmark { 
glam2scan -o $3_$4mers_$2.out p $3_$4mers_$2.glam2 $1; 
} 
/usr/bin/time -o ${NAME}_${KSIZE}mers_times.txt -a -f 'real %E\\tuser %U\\tsys %S\\tmem %M' glam2_benchmark ${DB} ${QUERY}${m} ${NAME} ${KSIZE}  
