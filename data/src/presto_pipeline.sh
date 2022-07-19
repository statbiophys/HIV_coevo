#!/bin/bash




fastq_dir="../fastq/abr/"
clones_dir="../abr_clones/"

sing_file_pipeline() {

	file_name=${1}

	echo "Unzipping..."
	zcat ${fastq_dir}${file_name}_1.fastq.gz > ${fastq_dir}${file_name}_1.fastq
	zcat ${fastq_dir}${file_name}_2.fastq.gz > ${fastq_dir}${file_name}_2.fastq


	FilterSeq.py quality -s ${fastq_dir}${file_name}_1.fastq -q 30 --outname ${file_name}_1_FS --log FS1.log
	FilterSeq.py quality -s ${fastq_dir}${file_name}_2.fastq -q 30 --outname ${file_name}_2_FS --log FS2.log
	rm -f *log
	rm ${fastq_dir}${file_name}_1.fastq
	rm ${fastq_dir}${file_name}_2.fastq


	PairSeq.py -1 ${fastq_dir}${file_name}_1_FS_quality-pass.fastq -2 ${fastq_dir}${file_name}_2_FS_quality-pass.fastq
	rm -f ${fastq_dir}${file_name}_1_FS_quality-pass.fastq
	rm -f ${fastq_dir}${file_name}_2_FS_quality-pass.fastq


	AssemblePairs.py align -1 ${fastq_dir}${file_name}_2_FS_quality-pass_pair-pass.fastq -2 ${fastq_dir}${file_name}_1_FS_quality-pass_pair-pass.fastq \
		--coord presto --rc tail --outname ${file_name} --log AP.log
	rm -f ${fastq_dir}${file_name}_1_FS_quality-pass_pair-pass.fastq
	rm -f ${fastq_dir}${file_name}_2_FS_quality-pass_pair-pass.fastq
	rm -f *log


	CollapseSeq.py -s ${fastq_dir}${file_name}_assemble-pass.fastq -n 0 --inner --uf PRCONS --cf CONSCOUNT --act sum --outname ${file_name}
	rm -f ${fastq_dir}${file_name}_assemble-pass.fastq


	#ParseHeaders.py table -s ${fastq_dir}${file_name}_collapse-unique.fastq -f ID DUPCOUNT


	seqtk="/home/andrea/Documents/Immunology/software/seqtk/seqtk"

	$seqtk seq -a ${fastq_dir}${file_name}_collapse-unique.fastq > ${fastq_dir}${file_name}_assemble.fasta
	rm -f ${fastq_dir}${file_name}_collapse-unique.fastq

	AssignGenes.py igblast -s ${fastq_dir}${file_name}_assemble.fasta -b /home/andrea/.local/share/igblast --organism human --loci ig --format blast

	MakeDb.py igblast -i ${fastq_dir}${file_name}_assemble_igblast.fmt7 -s ${fastq_dir}${file_name}_assemble.fasta -r /home/andrea/.local/share/germlines/imgt/human/vdj
	rm -f ${fastq_dir}${file_name}_assemble_igblast.fmt7
	rm -f ${fastq_dir}${file_name}_assemble.fasta


	python3 final_parser.py ${fastq_dir}${file_name}_assemble_igblast_db-pass.tsv ${clones_dir}${file_name}_clones.tsv
	rm -f ${file_name}_assemble_igblast_db-pass.tsv
}


for entry in ${fastq_dir}*
do
	if [[ "$entry" == *abr_1.fastq.gz ]]
	then
		name=$(echo $entry | cut -d/ -f4 | cut -d. -f1 |  cut -d_ -f1-3)
		echo $name
		if [ -f "$clones_dir""$name"_clones.tsv ]
		then
			echo "file already present, skipped"
		else
			sing_file_pipeline $name
		fi
	fi
done
