#!/bin/bash 

# List of files to process
declare -a arr=("S15_Normal" "S15_SSL"
                "A2_Adenoma" "A2_Normal")

for i in "${arr[@]}"
do
  
  echo "i"

  #Step 1 Trim
  mkdir -p ./trimmed
  if [[ -f ./trimmed/"$i"_R1_val_1.fq.gz && -f ./trimmed/"$i"_R2_val_2.fq.gz ]]; then
    echo "Files already trimmed"
  else
    ~/TrimGalore-0.6.5/trim_galore --fastqc --clip_R2 5 --paired "$i"_R1.fastq.gz "$i"_R2.fastq.gz -o ./trimmed
  fi

  #Step 2 Align
  mkdir -p ./aligned
  if [[ -f ./aligned/"$i".bam ]]; then
    echo "File already aligned"
  else
    ~/Bismark-0.23.0/bismark --genome_folder '/references/hs/GRCh37/hg19' -1 ./trimmed/"$i"_R1_val_1.fq.gz -2 ./trimmed/"$i"_R2_val_2.fq.gz -B ./aligned/"$i" -p 4
  fi 
  
  #Step 3 deduplicate
  mkdir -p ./deduplicated
  if [[ -f ./deduplicated/"$i".deduplicated.bam ]]; then
    echo "File already deduplicated"
  else
    ~/Bismark-0.23.0/deduplicate_bismark -p --bam ./aligned/"$i"_pe.bam --output_dir ./deduplicated
  fi
  
  #Step 4 methylation calling
  mkdir -p ./meth_calls
  ~/Bismark-0.23.0/bismark_methylation_extractor -p --cytosine_report --comprehensive --genome_folder '/references/hs/GRCh37/hg19' --no_overlap --multicore 8 --buffer_size 2G ./deduplicated/"$i"_pe.deduplicated.bam -o ./meth_calls

done