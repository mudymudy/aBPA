#!/usr/bin/env nextflow

params.data = ""
reads = Channel.of(params.data)

params.output = ""
resultsDir = Channel.of(params.output)

params.gcompleteness = 50
geneCompleteness = Channel.of(params.gcompleteness)

params.coverage = 0.5
normalizedCoverage = Channel.of(params.coverage)

params.threads = 10
threadsGlobal = Channel.of(params.threads)

params.tax_id = 0
taxID = Channel.of(params.tax_id)

params.genomes = 100
downloadGenomes = Channel.of(params.genomes)

params.clustering = 0.95
cdHitCluster = Channel.of(params.clustering)

params.core = 0.01
pangenomeThreshold = Channel.of(params.core)

params.clean = "strict"
pangenomeMode = Channel.of(params.clean)

params.config = ""
configFile = Channel.of(params.config)

params.help = false


// Enable DSL2
nextflow.enable.dsl=2

def print_help() {
	println "\n\033[1;31mSYNOPSIS\033[0m"

	println "\n\033[1;33mUSAGE\033[0m"
	println "\nnextflow run aBPA.nf --data <PATH> --output <PATH> --tax_id <INT> --config <PATH/FILE> [..OPTIONS..]"
	
	println "\n\033[1;33mMANDATORY\033[0m"
	println "  --data <PATH>		Set data file PATH"
	println "  --output <PATH>		Set output directory PATH"
	println "  --tax_id <INT>		Set taxonomical ID value <INT>"
	println "  --config <PATH>		Set config file PATH"
	
	println "\n\033[1;33mOPTIONS\033[0m"
	println "  --threads <INT>		Set number of threads (default: 10)"
	println "  --gcompleteness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
	println "  --coverage <INT/FLOAT>	Set mean depth of coverage threshold (default: 0.5)"
	println "  --genomes <INT>		Set number of genomes to download (default: 100)"
	println "  --clustering <INT/FLOAT>	Set clustering threshold (default 0.95)"
	println "  --core-threshold <FLOAT>	Set core genome threshold (default: 0.01)"
	println "  --clean-mode <STRING>	Set pangenome mode (default: strict)"
	println "  --help			Print help page and exit"
	
	println "\n\033[1;31mDESCRIPTION\033[0m"
	println "\n\033[1;33m--data <PATH>\033[0m"
	println "Please specify the full PATH of your data. Example: /home/user/mydata/data"
	
	println "\n\033[1;33m--output <PATH>\033[0m"
	println "Please specify the full PATH of your output folder. You need to make the folder first before running the program."

	println "\n\033[1;33m--tax_id <INT>\033[0m"
	println "Please specify the taxonomical ID for your bacteria. It should be a discrete and unique number."
	
	println "\n\033[1;33m--config <PATH>\033[0m"
	println "\nPlease set file PATH of your config.tab file. Example: /home/user/me/aBPA/config/config.tab"
	println "config.tab file should contain 3 fields separated by tab. First field should have the sample name, second field softclipping value <INT> and third field group ID."
	println "\nExample: \n	SAMPLE1	5	NONUDG\n	SAMPLE2	2	UDG"

	println "\n\033[1;33m--threads <INT>\033[0m"
	println "Set amount of threads to be used globally."

	println "\n\033[1;33m--gcompleteness <INT>\033[0m"
	println "Set gene breadth of coverage threshold as percentage. Genes that have a value less than <INT> will be considered absent."

	println "\n\033[1;33m--coverage <INT/FLOAT>\033[0m"
	println "Set gene normalized coverage threshold. Currently aBPA is using the simplest statistic for normalization: (Gene mean depth/Global mean depth)."

        println "\n\033[1;33m--genomes <INT/FLOAT>\033[0m"
	println "Set amount of FASTA/GENBANK files to be downloaded. Bear in mind disk space."

        println "\n\033[1;33m--clustering <INT/FLOAT>\033[0m"
	println "Set clustering threshold <INT/FLOAT> for FASTA database.\nA value of 0.9 means any group of sequences with identity values equal or bigger than 0.9 will be clustered together and a consensus representative sequence will be produced."

        println "\n\033[1;33m--core-threshold <INT/FLOAT>\033[0m"
	println "Set threshold for core genome building. Similarly as clustering flag but during pangenome step."

        println "\n\033[1;33m--clean-mode <INT/FLOAT>\033[0m"
	println "Set behaviour of pangenome building. Visit Panaroo documentation to know more about this.\n\n"

    exit 0
}


if (params.help) {
    print_help()
}


/*
 * As the name suggest, it will just generate the whole pipeline directory structure.
 */

process dirStructure {

	input:
	path makeDir

	script:
	"""
	#!/bin/bash
	mkdir -p "${makeDir}/NCBI/GFF"
	mkdir -p "${makeDir}/NCBI/FASTA"
	mkdir -p "${makeDir}/CLUSTERING"
	mkdir -p "${makeDir}/PROKKA/GFF"
	mkdir -p "${makeDir}/PANGENOME"
	mkdir -p "${makeDir}/ALIGNMENTS"
	mkdir -p "${makeDir}/NORMALIZATION"
	mkdir -p "${makeDir}/MATRIX"
	mkdir -p "${makeDir}/PLOTS"
	mkdir -p "${makeDir}/HETEROPLASMY/intermediate_files"
	mkdir -p "${makeDir}/HETEROPLASMY/distributions"
	"""
}

/*
 * entrez() will download FASTA and GenBank files from NCBI as long as a taxonomical ID was provided (mandatory value)
 */

process entrez {
	conda "${projectDir}/envs/entrez.yaml"
        publishDir "${makeDir}/results/NCBI", mode: 'copy', overwrite: true

	input:
	val gE
	val txid
	path makeDir

	output:

	path '*.gz'


	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${txid}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | while read url; do
        
	if [ "\$counter" -ge "${gE}" ]; then
		break
	fi

	if [ -z "\$url" ]; then
		continue
	fi
        
	fname="\$(basename "\$url")"
	wget "\$url/\${fname}_genomic.gbff.gz"
        
	wget "\$url/\${fname}_genomic.fna.gz"

	counter="\$((counter + 1))"	
	done
	"""
}



process unzipFiles {
	input:
	path gz_files

	output:
	path "*.fna" , emit: fastaFiles
	path "*.gbff" , emit: gffFiles

	script:
	"""
	gunzip -f ${gz_files}
	
	"""
}

process fastaDatabase {
	conda "${projectDir}/envs/biopython.yaml"
	
	input:
	path gffFiles, stageAs: 'gff/*'

	output:
	path "clustered_sequences.fasta" , emit: theFastaDatabase

	script:

	"""
	parsing_and_contatenating.py gff

	"""
}


process clustering {

	conda "${projectDir}/envs/cdhit.yaml"

	input:
	path fastaDB
	val clustering
	val threadsGlobal

	output:
	path "clustered_non_redundant_genes.fasta", emit: clusteredDatabase
	
	script:
	"""
	#!/bin/bash
	cd-hit-est -i $fastaDB -o clustered_non_redundant_genes.fasta -c $clustering -n 10 -T $threadsGlobal
	"""
}

process prokkaMakeAnnotations {
	conda "${projectDir}/envs/prokka.yaml"

	input:
	path clusteredSeqsDB
	val threadsGlobal
	path gffFiles , stageAs: 'gff/*'
	path fastaFiles, stageAs: 'fasta/*'

	output:
	path '*_genomic' , emit: prokkaOut
	path 'filteredGFF/*gff', emit: prokkaGFF
	script:
	"""
	ls -l gff | awk 'NR==2{print \$NF}' > first.txt
	name=\$(cat first.txt | awk -F'/' '{print \$NF}')
	echo -e "\$name"
	species=\$(head -n 20 gff/"\$name" | grep "ORGANISM" | awk '{print \$2, \$3}' | sed -e 's/ /_/g')
	echo -e "\$species"
	for i in fasta/*; do
		name=\$(basename "\$i")
		prokka --outdir "\${name%.fna}" --addgenes  --addmrna --species "\$species" --proteins clusteredSeqsDB --force --cpus $threadsGlobal "\$i"
	done
	
	mkdir -p filteredGFF

	for sample in *_genomic; do
		name=\$(basename "\${sample%_genomic}")
		mv "\$sample"/*gff filteredGFF/"\${name}.gff"
	done

	"""
}

process makePangenome {
	conda "${projectDir}/envs/panaroo.yaml"
	
	input:
	path prokkaGFF, stageAs: 'filteredGFF/*'
	val pangenomeMode
	val pangenomeThreshold
	val threadsGlobal

	output:
	path 'pan_genome_reference.fa' , emit: panSequence
	path 'gene_presence_absence.Rtab' , emit: initialMatrix
	path 'aligned_gene_sequences/*' , emit: alignedGenesSeqs

	script:
	"""
	panaroo -i filteredGFF/*.gff -o ./ --clean-mode $pangenomeMode -a core --core_threshold $pangenomeThreshold -t $threadsGlobal

	"""
}


process  formattingPangenome {
	conda "${projectDir}/envs/seqtk.yaml"

	input:
	path panGenomeReference, stageAs: 'pan_genome_reference.fa'

	output:
	path 'panGenomeReference.fasta', emit: panGenomeReference

	script:
	"""
	seqtk seq $panGenomeReference > panGenomeReference.fasta
	"""
}


process alignment {
	conda "${projectDir}/envs/alignment.yaml"

	input:
	path reads
	path panRef, stageAs: 'panGenomeReference.fasta'
	val threadsGlobal
	

	output:
	path '*_DMC_P.bam', emit: postAlignedBams
	path '*_final.fastq', emit: postAlignedReads
	

	script:
	"""
	for sample in reads/*; do
		bwa index $panRef
		name=\$(basename "\$sample")
		bwa aln -l 16500 -n 0.01 -o 2 -t $threadsGlobal $pan_genome_reference.fasta "\$sample" > "\${name%.fastq*}.sai"
		bwa samse $pan_genome_reference.fasta "\${name%.fastq*}.sai" "\$sample" > "\${name%.fastq*}.sam"
		samtools view -bS "$output"/ALIGNMENTS/"${name%.fastq*}.sam" > "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
		samtools quickcheck "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
		samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
		samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam"
		samtools view -b -@ 10 -F 4 "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"
		samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"
		/miniforge3/envs/alignment/bin/bam trimBam "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam" -L "$softclipping" -R "$softclipping" --clip
		samtools view -q 25 -o "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam"
		samtools view -e 'length(seq)>34' -O BAM -o "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam"
		samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam"
		samtools coverage "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name}"_genomicsMetrics.txt
		samtools fastq -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_final.fastq"
	done
	"""
}


workflow {
	dirStructure(resultsDir)
	dwnld = entrez(downloadGenomes, taxID, resultsDir)
	unzipFiles(dwnld)
	concatenatedSeqs = fastaDatabase(unzipFiles.out.gffFiles)
	clustering(fastaDatabase.out.theFastaDatabase, cdHitCluster, threadsGlobal)
	prokkaMakeAnnotations(clustering.out.clusteredDatabase, threadsGlobal, unzipFiles.out.gffFiles, unzipFiles.out.fastaFiles)
	makePangenome(prokkaMakeAnnotations.out.prokkaGFF, pangenomeMode, pangenomeThreshold, threadsGlobal)
	formattingPangenome(makePangenome.out.panSequence)
	alignment(reads, formattingPangenome.out.panGenomeReference, threadsGlobal,  )
}
