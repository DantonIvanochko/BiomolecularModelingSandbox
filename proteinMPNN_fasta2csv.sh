#!/bin/bash

# convert proteinmpnn output to csv format

output_dir="./ProteinMPNN_outputs/"

for output_fa_dir in ${output_dir}/temp*Output/seqs/;
do
for output_fa_file in ${output_fa_dir}*.fa;
do
    output_fa_file_basename=$(basename $output_fa_file .fa);
    echo "Converting ${output_fa_file} to csv.";
    sed 's/$/,/' ${output_fa_file} |
    paste -d " "  - - | # merge alternating lines
    sed -e 's/\s\+//g' | # remove whitespace
    sed 's/>//g' | # remove > symbol from fasta format
    sed 's/T=//g' | # remove > symbol from fasta format
    sed 's/sample=//g' | # remove > symbol from fasta format
    sed 's/score=//g' | # remove > symbol from fasta format
    sed 's/global_//g' | # remove > symbol from fasta format
    sed 's/seq_recovery=//g' | # remove > symbol from fasta format
    sed "s/^/${output_fa_file_basename},/" | # add input file to first position of each line
    tail -n +2 >> ${output_fa_dir}/${output_fa_file_basename}_output.csv # sort on the 4th column (score) and write to csv
    sort -t, -k4,4 -n ${output_fa_dir}/${output_fa_file_basename}_output.csv -o ${output_fa_dir}/${output_fa_file_basename}_output.csv; # sort on the 4th column (score)
    sed -i '1s/^/input,T,sample,score,global_score,seq_recovery,seq,\n/' ${output_fa_dir}/${output_fa_file_basename}_output.csv # add column names and write to csv
done
done

