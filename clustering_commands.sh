#!/bin/bash
id=$1 # 90
len=$2 # 80
clustid='id'${id}'len'${len}

## conda environment file: ortho.yml 
## rename sequences to include isolate ID at beginning 
for file in results/*/bakta/*_bakta.faa; do str=$( echo $file | cut -f 2 -d '/' ); sed 's/^>/>'"$str"' /g' $file | sed 's/ /_/1' > tmp_faa/${str}.faa; done

## run orthofinder on all isolates 
orthofinder -f tmp_faa/ -t 24 

## cluster sequences for phylophlan
## create single fasta file 
cat tmp_faa/*.faa > all.faa
## cluster sequences at 90% identity over 80% of length, use -d 0 to make sure output includes isolate and protein names for each sequence
cd-hit -i all.faa -o clust.${clustid} -c 0.${id} -s 0.${len} -n 5 -d 0 
## count number of proteins per isolate in each cluster
for ((i=0;i<=$( grep '^>' clust.${clustid}.clstr | tail -n 1 | cut -f 2 -d ' ' );i++)); do sed -n '/^>Cluster '"$i"'$/,/^>Cluster /p' clust.${clustid}.clstr | grep -v '^>' | cut -f 2 -d '>' | cut -f 1 -d '_' | sort | uniq -c | sort -n > ${clustid}/counts/clust_${clustid}.${i}.prots.txt; done
## get number of isolates represented in each cluster 
wc -l ${clustid}/counts/clust_${clustid}.* | grep -v total | sed 's/^\s\+//' | sort -n | sed 's,${clustid}/,\t,g' | sed 's,.prots.txt,,g' > ${clustid}/counts_per_clust_${clustid}.txt

## identify clusters in the pangenome for different genome types: core99-100iso, soft95-99iso, shell15-95iso, cloud0-15iso, 95iso, unique, accessory
## core genome: 95-100% of isolates - iso95
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" '$1>=v*0.95' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/clust.${clustid}.95iso.clstr
## hard core genome: 99-100% of isolates - core99-100iso
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" '$1>=v*0.99' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.core99-100iso.clstr
## soft core genome: 95-99% of isolates - soft95-99iso
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" 'v*0.99>$1 && $1>=v*0.95' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.soft95-99iso.clstr
## shell genome: 15-95% of isolates - shell15-95iso
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" 'v*0.95>$1 && $1>=v*0.15' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.shell15-95iso.clstr
## cloud genome: 0-15% of isolates - cloud0-15iso
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" '$1<v*0.15' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.cloud0-15iso.clstr
## unique genes (in 1 strain) - unique
awk '$1==1' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.unique.clstr
## accessory genes (in >2 strains but not all) - accessory
awk -v v="$( ls ../results/*/bakta/*_bakta.faa | wc -l )" 'v>$1>2' ${clustid}/counts_per_clust_${clustid}.txt | cut -f 2 -d '.' | while IFS= read -r line; do rmln=$(($line+1)); sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | sed 's/^>Cluster '"$rmln"'$//g'; done > ${clustid}/pangenome/clust.${clustid}.accessory.clstr

## extract sequences for clusters in pangenome for different genome types
for file in ${clustid}/clust.${clustid}.*.clstr; do
    gen=$( echo $file | cut -f 3 -d '.' )
    ## get cluster number ID for cluster in genome
    grep '^>Cluster' ${clustid}/clust.${clustid}.${gen}.clstr | cut -f 2 -d ' ' > ${clustid}/clustID_${gen}_clust_${clustid}.txt
    ## get all protein sequences for each cluster in genome
    while IFS= read -r line; do 
        rmln=$(($line+1))
        sed -n '/^>Cluster '"$line"'$/,/^>Cluster '"$rmln"'/p' clust.${clustid}.clstr | grep -v '^>' | cut -f 2 -d '>' | cut -f 1 -d '.' | seqtk subseq all.faa - > ${clustid}/${gen}/clust_${clustid}.${line}.faa
    done < ${clustid}/clustID_${gen}_clust_${clustid}.txt
    ## get representative sequence for each cluster in genome
    grep '^>' ${clustid}/${gen}/clust_${clustid}.*.faa | cut -f 2 -d '>' | seqtk subseq clust.${clustid} - > ${clustid}/clust.${clustid}.${gen}_rep.faa
done

## conda environment: phylo.yml
## create custom phylophlan database from representative sequences of core genome (95iso) 
phylophlan_setup_database -i core_clust/${clustid}/clust.${clustid}.95iso_rep.faa -o phylophlan_databases/${clustid}.95iso_rep/ -d ${clustid}.95iso_rep -e .faa -t a
diamond makedb --in phylophlan_databases/${clustid}.95iso_rep/${clustid}.95iso_rep.faa --db phylophlan_databases/${clustid}.95iso_rep/${clustid}.95iso_rep
## create phylophlan config file
phylophlan_write_config_file --overwrite \
    -o phylophlan_data/custom_config_coreaa.cfg \
    -d a \
    --db_aa diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml 
## run phylophlan
phylophlan -i orthofinder_out/tmp_faa/ \
    --proteome_extension "faa" \
    -o phylophlan_out \
    -d ${clustid}.95iso_rep \
    --databases_folder phylophlan_databases \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    --min_num_entries 85 \
    -t a \
    -f phylophlan_data/custom_config_coreaa.cfg \
    --diversity low \
    --nproc 12 \
    --verbose 2>&1 | tee phylophlan_out/logs/phylophlan__output_isolates.log

## to run:
# sh clustering_commands.sh 90 80