#!/usr/bin/env nextflow

params.data = ""
reads = Channel.of(params.data)

params.output = ""
resultsDir = Channel.of(params.output)

params.gcompleteness = 50
geneCompleteness = Channel.of(params.gcompleteness)

params.coverageDown = 0.5
normalizedCoverageDown = Channel.of(params.coverageDown)

params.coverageUp = 4
normalizedCoverageUp = Channel.of(params.coverageUp)

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
	path configFile

	output:
	path '*_DMC_P.bam', emit: postAlignedBams
	path '*_final.fastq', emit: postAlignedReads
	

	script:
	"""
	for sample in $reads/*; do
		bwa index $panRef
		name=\$(basename "\$sample")
		softClip=\$(grep "\$name" $configFile | awk '{print \$2}')
		bwa aln -l 16500 -n 0.01 -o 2 -t $threadsGlobal $panRef "\$sample" > "\${name%.fastq*}.sai"
		bwa samse $panRef "\${name%.fastq*}.sai" "\$sample" > "\${name%.fastq*}.sam"
		samtools view -bS "\${name%.fastq*}.sam" > "\${name%.fastq*}.bam"
		samtools quickcheck "\${name%.fastq*}.bam"
		samtools sort -o "\${name%.fastq*}_sorted.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}.bam"
		samtools index "\${name%.fastq*}_sorted.bam"
		rm "\${name%.fastq*}.bam"
		samtools view -b -@ 10 -F 4 "\${name%.fastq*}_sorted.bam" > "\${name%.fastq*}_sorted_mappedreads.bam"
		samtools index "\${name%.fastq*}_sorted_mappedreads.bam"
		bam trimBam "\${name%.fastq*}_sorted_mappedreads.bam" "\${name%.fastq*}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
		samtools view -q 25 -o "\${name%.fastq*}_qc.bam" "\${name%.fastq*}_softclipped.bam"
		samtools view -e 'length(seq)>34' -O BAM -o "\${name%.fastq*}_lg.bam" "\${name%.fastq*}_qc.bam"
		samtools sort -o "\${name%.fastq*}_DMC_P.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}_lg.bam"
		samtools coverage "\${name%.fastq*}_DMC_P.bam" > "\${name}"_genomicsMetrics.txt
		samtools fastq -@ $threadsGlobal "\${name%.fastq*}_DMC_P.bam" > "\${name%.fastq*}_final.fastq"
	done

	rm *sam *sai *_lg.bam *_qc.bam *_sorted_mappedreads.bam*
	"""
}


process alignmentSummary {
	conda "${projectDir}/envs/alignment.yaml"
	
	input:
	path configFile
	path bamfiles, stageAs: 'bam/*'
	

	output:
        path 'postPangenomeAlignment*bam' , emit: postAlignmentFiles
	path 'completenessSummary.tab', emit: completenessSummary
	path '*_refLength.txt', emit: refLenght
	path '*_rawCoverage.txt' , emit: rawCoverage


	script:
	"""
	#!/bin/bash
	awk '{print \$NF}' $configFile | uniq > groups.txt

	while read -r groupID; do
		groupName=\$(echo "\$groupID")
		grep -w "\$groupID" $configFile | awk '{print \$1}' > "\$groupName"ID
	done < groups.txt

        for groupFile in *ID; do
                IDs=\$(basename "\${groupFile%ID}")

		echo "Group ID: \$IDs"

                bamFiles=()

                while read -r sampleName; do
                        bamFiles+=(./bam/"\${sampleName%.fastq*}_DMC_P.bam")
                        echo "Adding BAM file: ./bam/\${sampleName%.fastq*}_DMC_P.bam"
		done < "\$groupFile"

                if [ \${#bamFiles[@]} -eq 1 ]; then
                        cp "\${bamFiles[0]}" postPangenomeAlignment_"\${IDs}".bam
                elif [ \${#bamFiles[@]} -gt 1 ]; then
                        samtools merge postPangenomeAlignment_mergedGroup"\${IDs}".bam "\${bamFiles[@]}"
                fi
        done

        for i in postPangenomeAlignment*bam; do
                samplename=\$(basename ./bam/"\${i%.bam}")
                samtools index "\$i"
                samtools depth -a "\$i" > "\${samplename}_rawCoverage.txt"
                samtools idxstats "\$i" | awk '{sum += \$2} END {print sum}' > "\${samplename}_refLength.txt"
                samtools coverage "\$i" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> completenessSummary.tab
        done
	"""
}

process normalizationFunction {

	input:
	path refLength, stageAs: 'refLength/*'
	path rawCoverage, stageAs: 'rawCoverage/*'


	output:
	path 'geneNormalizedSummary.txt', emit: geneNormalizedSummary
	path 'globalMeanCoverage.txt' , emit: globalMeanCoverage

	script:
	"""
	#!/bin/bash

	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > geneNormalizedSummary.txt
	echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > globalMeanCoverage.txt

	for i in rawCoverage/*_rawCoverage.txt; do
		name=\$(basename "\${i%_rawCoverage.txt}")

		#Compute global mean coverage
		globalMean=\$(awk -v name="\$name" '{sum += \$3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check log file"}' "\$i")
		finalCount=\$(awk -v name="\$name" '{fcount++} END {print fcount}' "\$i")
		refCount=\$(cat refLength/"\${name}_refLength.txt")
		echo -e "\$name\t\$finalCount\t\$refCount\t\$globalMean" >> globalMeanCoverage.txt

		#Normalize coverage per gene
		awk -v globalMean="\$globalMean" -v name="\$name" -v sampleCoverage="\$finalCount" -v refCount="\$refCount" '
			{
			if (\$2 > geneLength[\$1]) {
				geneLength[\$1] = \$2
			}
			sumgene[\$1] += \$3
			countgene[\$1]++
			}
		END {
			for (gene in geneLength) {
				geneMean = sumgene[gene] / countgene[gene]
				normalizedGeneScaled = (geneMean / globalMean) * geneLength[gene]
				normalizedGeneSimple = (geneMean / globalMean)
				normalizedGenomeSimple = (geneMean / globalMean) * (geneLength[gene] / sampleCoverage)
				normalizedGenomeScaled = (geneMean / globalMean) * (geneLength[gene] / refCount)
				print name"\t"gene"\t"normalizedGeneSimple"\t"normalizedGeneScaled"\t"normalizedGenomeSimple"\t"normalizedGenomeScaled
			}
		}
		' "\$i" >> geneNormalizedSummary.txt
	done
	"""
}

process updateNormalization {

	input:
	path normalized, stageAs: 'normalized/*'
	path completeness, stageAs: 'completeness/*'


	output:
	path 'geneNormalizedUpdated.tab', emit: geneNormalizedUpdated


	script:
	"""
	#!/bin/bash
	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > geneNormalizedUpdated.tab

	sed -i -e 's/~/_/g' normalized/geneNormalizedSummary.txt

	awk 'NR>1{print \$1"XYZ"\$2, \$3, \$4, \$5, \$6}' normalized/geneNormalizedSummary.txt > TMP1

	awk '{print \$1"XYZ"\$2, \$3}' completeness/completenessSummary.tab > TMP2

	while read -r ID completeness;do

		if grep -wq "\${ID}" TMP1; then
			oldLine=\$(grep -w "\${ID}" TMP1)
			specificCompleteness=\$(grep -w "\${ID}" TMP2 | awk '{print \$NF}')
			echo -e "\${oldLine}\t\${specificCompleteness}" >> geneNormalizedUpdated.tab
		fi

	done < TMP2

	sed -i -e 's/XYZ/\t/g' geneNormalizedUpdated.tab
	sed -i -e 's/ /\t/g' geneNormalizedUpdated.tab

	rm TMP1 TMP2 
	"""
}

process applyCoverageBounds {
	
	input:
	path geneNormalizedUpdated, stageAs: 'geneNormalizedUpdated.tab'
	val normalizedCoverageDown
	val normalizedCoverageUp

	output:
	path 'geneNormalizedUpdatedFiltered.tab', emit: geneNormalizedUpdatedFiltered

	script:
	"""
	awk -v UpBound="\$normalizedCoverageUp" '\$3 < \$UpBound {print \$0}' $geneNormalizedUpdated > TMP1
	awk -v DownBound="\$normalizedCoverageDown" '\$3 > \$DownBound {print \$0}' TMP1 > geneNormalizedUpdatedFiltered.tab
	"""
}

process plotCoveragevsCompleteness {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdated, stageAs: 'gNS/'
	val gcompleteness
	val coverage

	output:
	path 'plotCoverage_vs_Completeness.png', emit: plotCoverage_vs_Completeness
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdated.tab $gcompleteness $coverage
	"""
}





process plotCoveragevsCompletenessOnFiltered {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdated, stageAs: 'gNS/'
	val gcompleteness
	val coverage

	output:
	path 'plotCoverage_vs_Completeness.png', emit: plotCoverage_vs_Completeness
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdated.tab $gcompleteness $coverage
	"""
}







process makeMatrix {
	conda "${projectDir}/envs/plot.yaml"

	input:
	path pangenomeRtab, stageAs: 'pangenome/*'
	path gMC, stageAs: 'gMC/*'
	path normalized, stageAs: 'normalized/*'

	output:
	path 'matrix.tab', emit: matrix
	path '*_final.csv', emit: finalCsv
	path 'sample_names', emit: sampleNames
	path 'INDEX', emit: INDEX

	script:
	"""
	awk 'NR==1{print \$0}' pangenome/gene_presence_absence.Rtab > matrix.tab
	awk 'NR>1 {print \$0}' pangenome/gene_presence_absence.Rtab | sort -k 1 -t \$'\t' >> matrix.tab
	awk 'NR>1 {print \$1}' pangenome/gene_presence_absence.Rtab | sort -k 1 -t \$'\t' > INDEX

	awk 'NR>1 {print \$1}' gMC/globalMeanCoverage.txt > sample_names

	while read -r name; do

		echo -e "Gene\tnormalizedCoverage\tcompleteness" > "\${name}"_index.tmp
		grep -w "\$name" normalized/geneNormalizedUpdated.tab | awk '{print \$2, \$3, \$NF}' >> "\${name}"_index.tmp

	done < sample_names


	for i in *_index.tmp; do

		sed -i -e 's/ /\t/g' "\$i"
		lambda.py "\$i"

	done
	"""
}


process buildHeatmap {
	conda "${projectDir}/envs/heatmap.yaml"

	input:
	path fCSV, stageAs: 'fCSV/*'
	path INDEX, stageAs: 'INDEX/*'
	path matrix, stageAs: 'matrix/*'
	path names, stageAs: 'names/*'

	output:
	path 'final_matrix.tab', emit: finalMatrix
	path 'presenceAbsence*.png', emit: presenceAbsence
	path 'maskedMatrixGenesOnlyAncient.txt', emit: maskedMatrixGenesOnlyAncient
	path 'maskedMatrixGenesUbiquitous.txt', emit: maskedMatrixGenesUbiquitous
	path 'maskedMatrixGenesNoUbiquitous.txt', emit: maskedMatrixGenesNoUbiquitous
	path 'sampleOrdernoUbiquitous.txt', emit: sampleOrdernoUbiquitous 
	path 'sampleOrderonlyAncient.txt', emit: sampleOrderonlyAncient

	
	script:
	"""
	for i in fCSV/*_final.csv; do

		name=\$(basename "\$i")
		sed  -e 's/,/\t/g' "\$i" | awk 'NR>1{print \$0}' > "\${name%_index.tmp_final.csv}"_INDEX.Z

	done


	for i in *_INDEX.Z; do
		name=\$(basename "\$i")
		#Create the FINAL_INDEX file
		echo "\${name%_INDEX.Z}" > "\${name}"_FINAL_INDEX
    
		#Process the INDEX file
		while read -r gene; do
			toprint=\$(echo "\$gene 0")
			if grep -wq "\$gene" "\$i"; then
				grep -w "\$gene" "\$i" >> "\${name}"_FINAL_INDEX
			else
				echo "\$toprint" >> "\${name}"_FINAL_INDEX
			fi
		done < INDEX/INDEX
    
		#Extract the last column
		awk '{print \$NF}' "\${name}"_FINAL_INDEX > "\${name}"_FINALCOLUMN

	done


	paste matrix/matrix.tab *_FINALCOLUMN > final_matrix.tab
	tr '\n' ' ' < names/sample_names > names_heatmap

	heatmap.py final_matrix.tab names_heatmap
	"""
}


process makeConsensus {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'
	path bamFiles, stageAs: 'BAM/*'

	output:
	path 'extractedSequences*.fq', emit: extractedSequencesReads 
	path 'extractedSequences*.fasta', emit: extractedSequencesFasta

	script:
	"""
	for b in BAM/*; do
		basename=\$(basename "\$b")
		bcftools mpileup -f $panGenomeRef "\$b" | bcftools call -c | vcfutils.pl vcf2fq > extractedSequences"\${basename%.bam}".fq
		seqtk seq -a extractedSequences"\${basename%.bam}".fq > extractedSequences"\${basename%.bam}".fasta
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
	alignment(reads, formattingPangenome.out.panGenomeReference, threadsGlobal, configFile)
	alignmentSummary(configFile, alignment.out.postAlignedBams)
	normalizationFunction(alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage)
	updateNormalization(normalizationFunction.out.geneNormalizedSummary, alignmentSummary.out.completenessSummary)
	plotCoveragevsCompleteness(updateNormalization.out.geneNormalizedUpdated, geneCompleteness, normalizedCoverageDown)
	makeMatrix(makePangenome.out.initialMatrix , normalizationFunction.out.globalMeanCoverage, updateNormalization.out.geneNormalizedUpdated)
	buildHeatmap(makeMatrix.out.finalCsv, makeMatrix.out.INDEX ,makeMatrix.out.matrix, makeMatrix.out.sampleNames)
	makeConsensus(formattingPangenome.out.panGenomeReference, alignmentSummary.out.postAlignmentFiles)
	applyCoverageBounds(updateNormalization.out.geneNormalizedUpdated, normalizedCoverageDown, normalizedCoverageUp)
}
