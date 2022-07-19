# Instructions for dataset downloading processing


## Downloading fastq files


The open access database is available under the Sequence Read Archive, SRA, under the BioProject ID PRJNA543982.

https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA543982

The list of fastq files to download is **download\_metadata.csv**, which contains all the SSR identifiers.
The scripts that run through that table and download all the fastq files are **src/downloader_abr.sh** for the antibodies data and **src/downloader_hiv.sh** for the HIV sequences.
These scripts use *fastq-dump*, therefore the SRA toolkits needs to be present on the system (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).

Note that in this repository most of the files generated in the pipeline are not present because their large dimensions.
Only the final tables used in the analysis are saved: for the HIV they are contained in the **hiv_msa** directory, for the immune system in the **abr_fams** directory.


## Assembling the HIV sequences

The HIV sequence pairs are assembled from the fastq files using the python notebook script **src/hiv_assembler.ipynb**.

The required libraries for the script are *numpy*, *pandas*, *biopython* (https://biopython.org/).

The procedure for assembling the pairs is the following:
- The sequence 1 and the reverse complement of sequence 2 are locally aligned (using `pairwise2.align.localms` of biopython) with a large gap penalty.
-  Using the information of the alignment, the assembled read is composed as follows: the first part of the sequence 1 until the overlapping region, the overlapping region in which each position contains the letter with larger quality between the two sequences, the final non-overlapping part of the reversed complement of sequence 2.
- The assembled sequences with average sequencing quality smaller than 32 are discarded.
- The assembled sequences with length smaller than 400bp and larger that 425bp are discarded.
- The assembled sequences that count only one read (singletons) are discarded.

The obtained sequences are saved into the directory **hiv_seqs/** as fasta files.
The identifiers of the fasta for each unique sequence are: `sequenceIndex_count`.
The script writes the log file **hiv_seqs_info.tsv** which counts all the failed quality checks and the final count for each sample.
Note that for two samples, pat6_day2808 and pat2_day4656, the procedure fails to find any good sequence.



## Framing  and multiple aligning the HIV sequences


In order to track Single Nucleotide Polymorphism along the different time samples of a patient, the sequences need to be multiple aligned.
Moreover, to distinguish between synonymous and non-synonymous mutations the sequences need to be framed, i.e. identifying the proper codons in the sequence of nucleotides.
These two operations are performed in the python notebook **src/hiv_frame.ipynb**.
The framing operation is performed first.
For each sequence, the number of stop codons is counted for the three possible choices of the frame.
The frame that minimises the number of stop codons is then chosen.
The number of stop codons in such a frame is most of the time zero, as expected from productive sequences.
However, a small fraction of sequences with few stop codons remains.
The ones with more that two counts of stop codons are discarded.
The framed sequences are saved in **hiv_framed_seqs/**.
The file **hiv_seqs_info.tsv** is updated.

The Multiple Alignment Sequence among different time samples of a given patient is performed with Clustal Omega (http://www.clustal.org/omega/).
To preserve the codon structure in the multiple alignment, this operation is performed in the amino-acid space.
Therefore, each sequence is converted from nucleotides to amino-acids, then all the sequences of all the time points are aligned, and finally the gap-structure of the alignment is mapped back to the nucleotide sequences.
The resulting fasta files for each patient are written in the folder **hiv_msa/**.
The identifiers of those fasta files are `SequenceIndex_count_timeIndex`.



## Assembling and blasting the antibody clonotypes


The assemble has been done with PRESTO (https://presto.readthedocs.io/en/stable/) and IgBlast (https://changeo.readthedocs.io/en/stable/examples/igblast.html). The script is **src/presto_pipeline.sh**.
The operations performed in the script are:
- A quality filter, which discard reads with average quality less than 30.
- Assembling each read with its reverse complement.
- Collapse equal reads defining a counter of multiplicity.
- Blasting the sequence against the database of human antibodies and identifying the proper V gene, J gene, and CDR3.

The final files (saved in the directory **abr_clones** contain the list of all the clonotypes with the sequences and all the V and J gene information.


## Building the IgH lineages

For each patient, the IgH clonotypes have been cluster into lineages/families.
The procedure is performed over all the clonotypes at any time point of the patient.
They are first grouped together in different classes defined by a specific V gene, J gene (annotated with IgBlast at the previous step) and CDR3 length.
We call this groups VJL classes.
We discarded all the classes having length smaller than 30bp and also all the clonotypes having count equal to 1 (singletons).

Within each VJL class lineages are generated by single-linkage clustering with a threshold of 0.1 on the hamming distance rescaled by the CDR3 length.
This procedure and threshold is similar to the one suggested in \cite{gupta2017hierarchical} and used in \cite{nourmohammad2019fierce}.
The clustering is performed with a very fast prefix tree algorithm (\url{https://github.com/statbiophys/ATrieGC})
The resulting clusters define the B-cell receptor families/lineages.
We perform this operation in the notebook **src/abr_fam_builder.ipynb**, by merging together all the clonotypes at any time point of a given patient.
