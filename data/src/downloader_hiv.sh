#!/bin/bash



while read -r line
do
	id=$(echo "$line" | awk -F',' '{printf "%s", $1}' | tr -d '"')
	name=$(echo "$line" | awk -F',' '{printf "%s", $27}' | tr -d '"')
	type=$(echo "$line" | awk -F',' '{printf "%s", $30}' | tr -d '"')

	good_type="VIRAL RNA"
	if [ "$type" == "$good_type" ]; then
		day=$(echo "$name" | cut -d'_' -f 2 | tr -dc '0-9')
		pat=$(echo "$name" | cut -d'_' -f 1 | tr -dc '0-9')
		out_name="pat"$pat"_day"$day"_hiv"
		echo "$out_name"
		fastq-dump --split-files -O "../fastq/hiv/" $id
		mv ../fastq/hiv/"$id"_1.fastq ../fastq/hiv/"$out_name"_1.fastq
		mv ../fastq/hiv/"$id"_2.fastq ../fastq/hiv/"$out_name"_2.fastq
	fi

done < ../download_metadata.csv
