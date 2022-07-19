#!/bin/bash



while read -r line
do
	id=$(echo "$line" | awk -F',' '{printf "%s", $1}' | tr -d '"')
	name=$(echo "$line" | awk -F',' '{printf "%s", $27}' | tr -d '"')
	type=$(echo "$line" | awk -F',' '{printf "%s", $30}' | tr -d '"')

	good_type="TRANSCRIPTOMIC"
	if [ "$type" == "$good_type" ]; then
		day=$(echo "$name" | cut -d'_' -f 2 | tr -dc '0-9')
		pat=$(echo "$name" | cut -d'_' -f 1 | tr -dc '0-9')
		out_name="pat"$pat"_day"$day"_abr"
		echo fastq/abr/"$out_name"_1.fastq.gz

		if [ -f ../fastq/abr/"$out_name"_1.fastq.gz ]
		then
			echo "file already present, skipped"
		else
			echo "downloading"
			fastq-dump --split-files --gzip -O "../fastq/abr/" $id
			mv ../fastq/abr/"$id"_1.fastq.gz ../fastq/abr/"$out_name"_1.fastq.gz
			mv ../fastq/abr/"$id"_2.fastq.gz ../fastq/abr/"$out_name"_2.fastq.gz
		fi

	fi

done < ../download_metadata.csv
