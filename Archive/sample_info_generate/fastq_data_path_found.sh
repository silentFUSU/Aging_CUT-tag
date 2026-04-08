#!/bin/bash

# Output file where the results will be written
data=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_dir.csv
cleaned_file=$(mktemp)
cat "$data" | tr -d '\r' | awk '{gsub(/[\x00-\x1F\x7F]+/, ""); print}' > "$cleaned_file"  

output_file="/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/all/HiC_datapath.txt"

# Clear the output file if it exists
> "$output_file"

while IFS=$',' read -r sampleID data_path; do
  # Skip header line
  if [[ $sampleID == "SampleID" ]]; then
    continue
  fi
#   echo $sampleID
  # Execute find to search recursively for files starting with SampleID and ending with .gz
  files=($(find "$data_path" -type f -name "${sampleID}*.gz"))

  # Check for exactly two files found
  if [[ ${#files[@]} -eq 2 ]]; then
    echo -e "$sampleID\t${files[0]}\t${files[1]}" >> "$output_file"
  else
    echo -e "$sampleID\tnot Found\tnot Found" >> "$output_file"
  fi

done < $cleaned_file