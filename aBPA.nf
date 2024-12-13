#!/usr/bin/env nextflow

params.data = ""
reads = Channel.of(params.data)

params.output = ""
resultsDir = Channel.of(params.output)

params.completeness = 50
geneCompleteness = Channel.of(params.completeness)

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

params.core = 0.95
pangenomeThreshold = Channel.of(params.core)

params.clean = "strict"
pangenomeMode = Channel.of(params.clean)

params.config = ""
configFile = Channel.of(params.config)

params.outgroup = ""
outTax = Channel.of(params.outgroup)

params.trustedGenomes = false
def trustedDataChannel = params.trustedGenomes ? Channel.fromPath(params.trustedGenomes) : Channel.empty()

params.missing = 0.01
missingProb = Channel.of(params.missing)

params.gaps = 2
gapFraction = Channel.of(params.gaps)

params.seed = 16500
seedAlignment = Channel.of(params.seed)

params.mapq = 30
mappingQuality = Channel.of(params.mapq)

params.minlength = 34
minReadLength = Channel.of(params.minlength)

params.maxlength = 300
maxReadLength = Channel.of(params.maxlength)

params.genotyper = "gatk"
genotypeMethod = Channel.of(params.genotyper)


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
	println "  --completeness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
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

	println "\n\033[1;33m--completeness <INT>\033[0m"
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

	path '*fna', emit: fastaFiles
	path '*gbff', emit: gffFiles

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
	downloadMeAndCheck() {
		wget "\$url/\${fname}_genomic.gbff.gz"
	        wget "\$url/\${fname}_genomic.fna.gz"

		gunzip -f "\${fname}_genomic.gbff.gz"
		gunzip -f "\${fname}_genomic.fna.gz"


		if [ -f "\${fname}_genomic.gbff.gz" ] || [ -f "\${fname}_genomic.fna.gz" ]; then	
			rm -f "\${fname}_genomic.gbff.gz" "\${fname}_genomic.fna.gz"
			return 1
		fi

		return 0
	}

	while ! downloadMeAndCheck; do
		echo "Files were corrupted. Retrying"
		sleep 3
	done
	
	# IF AFTER DOING gunzip -f WE STILL FIND THE EXTENSION *gz, then remove both files with the same {fname} (even if is just "\${fname}_genomic.gbff.gz" or "\${fname}_genomic.fna.gz" or both)
	# If we remove any file, then:
	#	1. We don't add +1 to the counter (it doesn't make sense to add +1 if the files were corrupted)
	#	2. We download both *gbff.gz and *fna.gz again and we proceed to test the integrity again by checking if there is any *gz extension after gunzip.

	counter="\$((counter + 1))"	
	fna_count=\$(ls -1 *.fna 2>/dev/null | wc -l)
	gbff_count=\$(ls -1 *.gbff 2>/dev/null | wc -l)

	# Continue downloading until the count matches gE
	if [ "\$fna_count" -lt "\$gE" ] || [ "\$gbff_count" -lt "\$gE" ]; then
		echo "Still need more files. Current count: \$fna_count .fna files, \$gbff_count .gbff files."
		continue
	fi


	done
	cat .command.out >> entrez.log
	"""
}


process trustedData {

	input:
	path data

	output:
	path '*fna', emit: trustedFasta
	path '*gbff', emit: trustedGenBank

	script:
	"""
	cp $data/* ./

	"""
}

process fastaDatabase {
	conda "${projectDir}/envs/biopython.yaml"
	
	input:
	path gffFiles, stageAs: 'gff/*'
	path fastaFiles, stageAs: 'fasta/*'
		
	output:
	path 'clustered_sequences.fasta' , emit: theFastaDatabase
	path 'cleanedFasta/*fna', emit: validFasta
	path 'cleanedGff/*gbff', emit: validGff
	path 'fastaDatabase.log', emit: fastaDatabaseLogFile

	script:
	"""
	parseTest.py gff > parseTest.txt
	grep "is not a valid GenBank file" parseTest.txt | awk '{print \$1}' | sed -e 's/gff\\///g' > blackListed.txt
	
	while read -r removeMe; do
	rm fasta/"\${removeMe%gbff}fna" gff/*"\$removeMe"	
	done < blackListed.txt

	parsing_and_contatenating.py gff
	
	mv fasta cleanedFasta
	mv gff cleanedGff
	
	for genome in cleanedFasta/*; do
		name=\$(basename "\$genome")
		awk '!/^>/ { gsub(/[^ACGTN]/, "N") }1' "\$genome" > tmp && mv tmp "\$genome"
	done

	cat .command.out >> fastaDatabase.log
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
	path "clustering.log", emit: clusteringLog
	script:
	"""
	#!/bin/bash
	cd-hit-est -i $fastaDB -o clustered_non_redundant_genes.fasta -c $clustering -T $threadsGlobal -d 0 -g 1 -M 0
	cat .command.out >> clustering.log
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
	path '*_fromProkka' , emit: prokkaOut
	path 'filteredGFF/*gff', emit: prokkaGFF
	path 'prokkaMakeAnnotations.log', emit: prokkaLogfile

	script:
	"""
	ls -l gff | awk 'NR==2{print \$NF}' > first.txt
	name=\$(cat first.txt | awk -F'/' '{print \$NF}')
	echo -e "\$name"
	species=\$(head -n 20 gff/"\$name" | grep "ORGANISM" | awk '{print \$2, \$3}' | sed -e 's/ /_/g')
	echo -e "\$species"
	for i in fasta/*; do
		name=\$(basename "\$i")
		prokka --outdir "\${name%.fna}_fromProkka" --species "\$species" --proteins clusteredSeqsDB --rawproduct --cpus $threadsGlobal "\$i"
	done
	
	mkdir -p filteredGFF

	for sample in *_fromProkka; do
		name=\$(basename "\${sample%_fromProkka}")
		mv "\$sample"/*gff filteredGFF/"\${name}.gff"
	done

	cat .command.out >> prokkaMakeAnnotations.log
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
	path 'makePangenome.log', emit: panarooLog

	script:
	"""
	panaroo -i filteredGFF/*.gff -o ./ --clean-mode $pangenomeMode -a core --core_threshold $pangenomeThreshold -t $threadsGlobal

	cat .command.out >> makePangenome.log
	"""
}


process  formattingPangenome {
	conda "${projectDir}/envs/formattingPangenome.yaml"

	input:
	path panGenomeReference, stageAs: 'pan_genome_reference.fa'

	output:
	path 'panGenomeReference.fasta', emit: panGenomeReference
	path 'panGenomeReference.dict', emit: panGenomeReferenceDictionary
	path 'panGenomeReference.fasta.fai', emit: panGenomeReferenceIndex

	script:
	"""
	seqtk seq $panGenomeReference > panGenomeReference.fasta
	picard CreateSequenceDictionary -R panGenomeReference.fasta
	samtools faidx panGenomeReference.fasta
	"""
}


process alignment {
	conda "${projectDir}/envs/alignment.yaml"

	input:
	path reads
	path panRef, stageAs: 'panGenomeReference.fasta'
	val threadsGlobal
	path configFile
	val missingProb
	val seedAlignment
	val gapFraction
	val minReadLength
	val maxReadLength

	output:
	path '*_DMC_P.bam', emit: postAlignedBams
	path '*_final.fastq', emit: postAlignedReads
	

	script:
	"""
	for sample in $reads/*; do
    		bwa index $panRef
    		name=\$(basename "\$sample")
    		softClip=\$(grep "\$name" $configFile | awk '{print \$2}')
    
    		# Making read groups
    		rg_id="\${name%.fastq*}"  # sample name as id
    		rg_sm="\${name%.fastq*}" # sample name again
    		rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
    		rg_lb="lib1"            # group id
    		rg_pu="unit1"           # not sure what Ill put here

    		bwa aln -l $seedAlignment -n $missingProb -o $gapFraction -t $threadsGlobal $panRef "\$sample" > "\${name%.fastq*}.sai"
    		bwa samse -r "@RG\\tID:\$rg_id\\tSM:\$rg_sm\\tPL:\$rg_pl\\tLB:\$rg_lb\\tPU:\$rg_pu" \
		$panRef "\${name%.fastq*}.sai" "\$sample" > "\${name%.fastq*}.sam"
    
    		samtools view -bS "\${name%.fastq*}.sam" > "\${name%.fastq*}.bam"
    		samtools quickcheck "\${name%.fastq*}.bam"
    		samtools sort -o "\${name%.fastq*}_sorted.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}.bam"
    		samtools index "\${name%.fastq*}_sorted.bam"
    		rm "\${name%.fastq*}.bam"
    		samtools view -b -@ $threadsGlobal -F 4 "\${name%.fastq*}_sorted.bam" > "\${name%.fastq*}_sorted_mappedreads.bam"
    		samtools index "\${name%.fastq*}_sorted_mappedreads.bam"
    		bam trimBam "\${name%.fastq*}_sorted_mappedreads.bam" "\${name%.fastq*}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
    		samtools view -q 30 -o "\${name%.fastq*}_qc.bam" "\${name%.fastq*}_softclipped.bam"
    		samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name%.fastq*}_lg.bam" "\${name%.fastq*}_qc.bam"
    		samtools sort -o "\${name%.fastq*}_DMC_P.bam" -O bam -@ $threadsGlobal "\${name%.fastq*}_lg.bam"
    		samtools coverage "\${name%.fastq*}_DMC_P.bam" > "\${name}_genomicsMetrics.txt"
    		samtools fastq -@ $threadsGlobal "\${name%.fastq*}_DMC_P.bam" > "\${name%.fastq*}_final.fastq"
	done

	rm *sam *sai *_lg.bam *_qc.bam *_sorted_mappedreads.bam*
	cat .command.out >> alignment.log
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

	#Getting some stats AND Re group the data so genotyping is simpler
        for i in postPangenomeAlignment*bam; do
                samplename=\$(basename ./bam/"\${i%.bam}")
                samtools index "\$i"
                samtools depth -a "\$i" > "\${samplename}_rawCoverage.txt"
                samtools idxstats "\$i" | awk '{sum += \$2} END {print sum}' > "\${samplename}_refLength.txt"
                samtools coverage "\$i" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> completenessSummary.tab
		mv "\$i" ./"\${i%.bam}_TMP.bam"
		mv "\$i".bai ./"\${i%.bam.bai}_TMP.bam.bai"
		picard AddOrReplaceReadGroups I="\${i%.bam}_TMP.bam" O="\${samplename}.bam" RGLB="\${samplename}" RGSM="\${samplename}" RGPU=Illumina RGPL=ILLUMINA RGID="\${samplename}" RGDS="\${samplename}"
		samtools index "\${samplename}.bam"
	done

	rm *TMP.bam*

	cat .command.out >> alignmentSummary.log
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
	val completenessBound

	output:
	path 'geneNormalizedUpdatedFiltered.tab', emit: geneNormalizedUpdatedFiltered

	script:
	"""
	awk 'NR==1{print \$0}' $geneNormalizedUpdated > header
	awk -v UpBound=$normalizedCoverageUp '\$3 < UpBound {print \$0}' $geneNormalizedUpdated > TMP1
	awk -v DownBound=$normalizedCoverageDown '\$3 > DownBound {print \$0}' TMP1 > TMP2
	awk -v completenessBound=$completenessBound '\$NF> completenessBound {print \$0}' TMP2 > TMP3
	cat header TMP3 > geneNormalizedUpdatedFiltered.tab
	"""
}

process plotCoveragevsCompleteness {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdated, stageAs: 'gNS/'
	val completeness
	val coverage

	output:
	path 'plotCoverage_vs_Completeness.png', emit: plotCoverage_vs_Completeness
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdated.tab $completeness $coverage
	"""
}





process plotCoveragevsCompletenessOnFiltered {
	conda "${projectDir}/envs/plot.yaml"
	
	input:
	path geneNormalizedUpdatedFiltered, stageAs: 'gNS/*'
	val completeness
	val coverage

	output:
	path 'plotCoverageVsCompletenessFiltered.png', emit: plotCoverageVsCompletenessFiltered
	
	script:
	"""
	plot_cvg_vs_completeness.py gNS/geneNormalizedUpdatedFiltered.tab $completeness $coverage
	mv plotCoverage_vs_Completeness.png ./plotCoverageVsCompletenessFiltered.png
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
		grep -w "\$name" normalized/geneNormalizedUpdatedFiltered.tab | awk '{print \$2, \$3, \$NF}' >> "\${name}"_index.tmp

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
	path 'genesAbovePercentSeries.txt', emit: genesAbovePercentSeries
	path 'blackListedQualityChecked.txt', emit: blackListed

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
        # I need to add a quality control step right here. User samples can be false positives sometimes or just super low quality and have 0 genes after filtering
        # Then, black list unwanted samples and exclude them from the final_matrix.tab document.
        for checkSample in *_FINALCOLUMN; do
                sampleName=\$(basename "\${checkSample%_INDEX.Z_FINALCOLUMN}")
                counts=\$(grep -c "0" "\$checkSample")
                totalLines=\$(wc -l "\$checkSample" | awk '{print \$1 - 1}')
                proportion=\$(awk -v absence="\$counts" -v record="\$totalLines" 'BEGIN { print (absence / record ) }')         

                if  (( \$(awk -v p="\$proportion" 'BEGIN { print (p > 0.95) }' ) )); then
                        echo "\$sampleName" >> blackListedQualityChecked.txt
                fi
        done
	
	# Do this only if blacklisted
	if [[ -s blackListedQualityChecked.txt ]]; then
		while read -r removeMe; do
			mv "\${removeMe}_INDEX.Z_FINALCOLUMN" "\${removeMe}LowQualitySample"    
			grep -v "\${removeMe}" names/sample_names > names/sample_names.tmp
			mv names/sample_names.tmp names/sample_names
		done < blackListedQualityChecked.txt    
	else
		touch blackListedQualityChecked.txt
	
	fi

	paste matrix/matrix.tab *_FINALCOLUMN > final_matrix.tab
	tr '\n' ' ' < names/sample_names > names_heatmap

	heatmap.py final_matrix.tab names_heatmap
	
	"""
}


process bcftoolsConsensus {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'
	path bamFiles, stageAs: 'BAM/*'

	output:
	path 'extractedSequences*.fasta', emit: consensusSequences

	script:
	"""
	for b in BAM/*; do
		basename=\$(basename "\$b")
		bcftools mpileup -f $panGenomeRef "\$b" | bcftools call -c | vcfutils.pl vcf2fq > extractedSequences"\${basename%.bam}".fq
		seqtk seq -a extractedSequences"\${basename%.bam}".fq > extractedSequences"\${basename%.bam}".fasta
	done
	"""
}


process gatkConsensus {
	conda "${projectDir}/envs/gatk.yaml"

	input:
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'
	path bamFiles, stageAs: 'BAM/*'
	path panGenomeRefDictionary, stageAs: 'panGenomeRef.dict'
	path panGenomeReferenceIndex, stageAs: 'panGenomeRef.fasta.fai'
	output:
	path '*GenotypedNormalizedConsensusSeq.fasta', emit: gatkConsensusSequences
	path '*GenotypedNormalized.vcf.gz', emit: gatkGenotypes

	script:
	"""

	for b in BAM/*; do
		basename=\$(basename "\$b")
		gatk3 -T UnifiedGenotyper --min_base_quality_score 30 \
			--genotype_likelihoods_model BOTH --annotateNDA \
			--genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES \
			-I "\$b" -R panGenomeRef.fasta \
			-o "\${basename%.bam}"Genotyped.vcf 

		bcftools norm -f panGenomeRef.fasta "\${basename%.bam}"Genotyped.vcf > "\${basename%.bam}"GenotypedNormalized.vcf
		bgzip -i "\${basename%.bam}"GenotypedNormalized.vcf
		bcftools index "\${basename%.bam}"GenotypedNormalized.vcf.gz
		bcftools consensus -a N -M N -f panGenomeRef.fasta "\${basename%.bam}"GenotypedNormalized.vcf.gz -o "\${basename%.bam}"GenotypedNormalizedConsensus.fasta
		seqtk seq "\${basename%.bam}"GenotypedNormalizedConsensus.fasta > "\${basename%.bam}"GenotypedNormalizedConsensusSeq.fasta
	done
	"""
}







/*
 * First seqtk seq every gene.aln file, and then fix FASTA header. Add samples gene sequences to the right alignments. check sequence lenght and fill with Ns if lenght smaller than modern strains.
 * Finally check number of FASTA headers and if not == to $genomes then add the missing and fill sequence with Ns.
 */


process filterGeneAlignments {
	conda "${projectDir}/envs/seqtk.yaml"

	input:
	path genesAln, stageAs: 'genes/*'
        path extractedSequencesFasta, stageAs: 'sampleGenes/*'
	path fFiles, stageAs: 'FNA/*'
	val genomes
	path outgroupSeq, stageAs: 'outgroup'
	path blackListed, stageAs: 'blackListed.txt'

	output:
	path 'filteredGenes/*_Filtered.fasta', emit: genesAlnSeq
	path 'sampleNames.txt', emit: sampleNames
	path 'specialCases/*fasta', emit: specialCases

	script:
	"""
	#!/bin/bash

	mkdir -p AlnSeq/

	echo -e "Fixing FASTA headers and extension of sequences with seqtk in existing gene alignments"
	for file in genes/*.aln.fas; do
		name=\$(basename "\${file%.aln.fas}")
		seqtk seq "\${file}" | awk '/^>/ {sub(/;.*/, "", \$0)} {print}' > AlnSeq/"\${name}_AlnSeq.fasta"
	done
	echo -e "Done\\n"

	echo -e  "Fixing FASTA headers and extension of sequences with seqtk in existing gene alignments, but for alignments ending with *fasta"
	for file in genes/*.fasta; do
		name=\$(basename "\${file%.fasta}")
		seqtk seq "\${file}" | awk '/^>/ {sub(/;.*/, "", \$0)} {print}' > AlnSeq/"\${name}_AlnSeq.fasta"
	done
	echo -e "Done\\n"


	echo -e "Replace ~ characters with _ in gene MSA filenames"
	for i in AlnSeq/*_AlnSeq.fasta; do
    		# Create a temporary name with '_TMP' to avoid overwriting
    		tmpname="\${i}_TMP"
    		mv "\$i" "\$tmpname"
		    
	    	# Fix the name by replacing ~ with _
	    	fixedname=\$(basename "\${tmpname%_TMP}" | sed -e 's/~/_/g')
	    	
    		# Rename the file with the fixed name
    		mv "\$tmpname" "AlnSeq/\$fixedname"
	done
	echo -e "Done\\n"
	

	echo -e "BlackList low quality user samples and then make a file with user sample names"
	if [[ -s blackListed.txt ]]; then
		while read -r removeMe; do
			mv sampleGenes/extractedSequences"\${removeMe}.fasta" sampleGenes/"\${removeMe}blackListed.txt"
			echo ""\${removeMe}" has been removed from analysis due to low quality."
		done < blackListed.txt
	else
		echo -e "Every sample passed quality checks."
	fi
	echo -e "Done\\n"

	echo -e "Collect sample names after quality checks"
	for sample in sampleGenes/*fasta; do
		sampleName=\$(basename "\${sample%.fasta}")
		echo "\${sampleName}" >> userSampleNames.txt
		echo -e ""\${sampleName}" has been added to userSampleNames.txt"
	done
	echo -e "Done\\n"


	echo -e "Adding user sample genes sequences to each particular gene MSA and replace gene name with sample name"
        for i in AlnSeq/*_AlnSeq.fasta; do
                name=\$(basename "\${i%_AlnSeq.fasta}")
                for sample in sampleGenes/*fasta; do
			sed -i -e 's/~/_/g' "\$sample"
			sampleName=\$(basename "\${sample%.fasta}")
			echo -e "Adding "\${sampleName}" gene sequences into "\${name}" MSA file"
                        grep -w -A 1 "\$name" "\$sample" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${i}"
                done
        done
	echo -e "Done\\n"

	echo -e "Adding outgroup gene sequences into each gene MSA file"
	sed -i -e 's/~/_/g' outgroup
	for i in AlnSeq/*_AlnSeq.fasta; do
		name=\$(basename "\${i%_AlnSeq.fasta}")
		echo -e "Adding outgroup sequences into "\${name}" MSA"
		grep -w -A 1 "\$name" outgroup | awk -v outgroup="outgroup" '/^>/ {sub(/^>.*/, ">" outgroup, \$0)} {print}' >> "\${i}"
	done
	echo -e "Done\\n"

	echo -e "Taking the total lenght of the gene sequence and then fill incomplete user samples gene sequences with n"
	for i in AlnSeq/*_AlnSeq.fasta; do
		geneName=\$(basename "\$i")
		numberOfColumns=\$(awk 'NR==2 {print length}' "\$i")
		awk -v numCols="\$numberOfColumns" '{
			if (\$0 ~ /^>/) {
				print
			} else {
			while ( length(\$0) < numCols) {
				\$0 = \$0 "n"
			}
			print
			}
		}' "\$i" > tmp && mv tmp "\${i}"
	done
	echo -e "Done\\n"


	# Make a file with downloaded modern genomes names
	for i in FNA/*.fna; do
		fnames=\$(basename "\${i%.fna}")
		echo "\${fnames}" >> modernSampleNames.txt
	done
	echo outgroup >> modernSampleNames.txt
	
	# Make a file with every sample combined
	cat modernSampleNames.txt userSampleNames.txt > sampleNames.txt


	echo -e "Checking if there are Panaroo headers artifacts"
	for file in AlnSeq/*_AlnSeq.fasta ; do
		geneName=\$(basename "\${file}")
		echo -e "Reading \$geneName MSA"
		
		while read -r sampleName; do
			newVariableName=">\$sampleName"
			matches=\$(grep "\$sampleName" "\$file")

			if [[ -n "\$matches" ]]; then
	
				while IFS= read -r matchedLine; do
	
					if [[ "\$matchedLine" != "\$newVariableName" ]]; then

						sed -i -e "s/\${matchedLine}/\${newVariableName}/g" "\$file"
						echo -e "\$matchedLine name was found but it should have been "\$newVariableName" instead. Fixed"
					fi
				done <<< "\$matches"
			else
				echo -e "It seems everything was okay for \$newVariableName in \$file."

			fi

		done < modernSampleNames.txt
	done
	echo -e "Done\\n"


	echo -e "If there are missing modern strain, append the sample and fill it with -"
	#(gaps;absence of gene; because we trust modern genomes assemblies?)
	for file in AlnSeq/*_AlnSeq.fasta ; do

		sampleValue=\$(awk '/^>/ {print \$0}' "\$file" | wc -l)
                numberOfColumns=\$(awk 'NR==2 {print length}' "\$file")
		totalSamples=\$(wc -l < sampleNames.txt)
		
		if (( sampleValue < totalSamples)); then
			
			while read -r strain; do
				if ! grep -wq "\$strain" "\$file"; then
					echo ">\$strain" >> "\$file"
					fakeSeq=\$(printf '%*s' "\$((numberOfColumns))" | tr ' ' '-')
					echo "\$fakeSeq" >> "\$file"
				fi 
			done < modernSampleNames.txt
		fi
	done
	echo -e "Done \\n"


	echo -e "If missing user sample == true, then append it and fill it with n's" 
	#(can't treat them as gaps because there is uncertainty)
        for file in AlnSeq/*_AlnSeq.fasta ; do

                sampleValue=\$(awk '/^>/ {print \$0}' "\$file" | wc -l)
                numberOfColumns=\$(awk 'NR==2 {print length}' "\$file")
                totalSamples=\$(wc -l < sampleNames.txt)
                
                if (( sampleValue < totalSamples)); then

                        while read -r strain; do
                                if ! grep -wq "\$strain" "\$file"; then
                                        echo ">\$strain" >> "\$file"
                                        fakeSeq=\$(printf '%*s' "\$((numberOfColumns))" | tr ' ' 'n')
                                        echo "\$fakeSeq" >> "\$file"
                                fi
                        done < sampleNames.txt
                fi
        done


	# Replacing ~ with _ in list of genes files
	for txtFile in *txt; do
		sed -i -e 's/~/_/g' "\$txtFile"
	done


	# Make directory to store special cases
	mkdir -p specialCases
	# Make directory where the final gene MSA will be stored
	mkdir -p filteredGenes
	
	echo -e "Check if there is a missing sample that should match perfectly in genes MSA, if not found, send MSA to special cases"
	for i in AlnSeq/*_AlnSeq.fasta; do
		geneName=\$(basename "\$i")
		samplesPresent=true

		# First I need to check if every sample is present or not.
		while read -r sample; do
			if ! grep -wq "\$sample" "\$i"; then
				samplesPresent=false
				break
			fi
		done < sampleNames.txt

		# If samplesPresent=True, then means that every sample was present so we can send the file to Filtered.
		if [ "\$samplesPresent" = true ]; then
			while read -r sample; do
				grep -w -A 1 "\$sample" "\$i" >> filteredGenes/"\${geneName}"	
			done < sampleNames.txt
		

		else
			# If missing sample, then move it to special cases
			mv "\$i" specialCases/
			echo -e ""\$i" file was found to have problems. Sent to specialCases"
		fi
	
	done
	echo -e "Done\\n"

	# After everything, check again if the number of samples in genes MSAs is different than the real number of samples, if true send those MSA to special cases
	echo -e "Looking for more special cases"
	for file in filteredGenes/*.fasta; do
		value=\$(awk '/^>/ {print \$0}' "\$file" | wc -l)
		sampleN=\$(wc -l < sampleNames.txt)
		if [ \$value -ne \$sampleN ]; then
			mv "\$file" specialCases/
			echo -e ""\$file" file was found to have problems. Sent to specialCases"
		fi
	done
	echo -e "Done\\n"


	# Dealing with fragmented Panaroo gene alignments multi-entries (nothing else to do here since is a panaroo problem, but I'll try to save as much as possible)
	for file in specialCases/*.fasta; do
		name=\$(basename "\${file%_AlnSeq.fasta}")
	# Identifying repeated headers and save them to .dpd 
		awk '/^>/ {count[\$0]++} END {for (header in count) if (count[header] > 1) print substr(header, 2)}' "\$file" > specialCases/"\${name}.dpd"

	# Extracting sequences for each repeated entry into temporary files
		while read -r entry; do
			awk -v sampleName="\$entry" '
			\$0 ~ "^>" sampleName {print_header=1; next}
			/^>/ {print_header=0}
			print_header {print}' "\$file" > specialCases/"\${name}Seqs_\${entry}"
		done < specialCases/"\${name}.dpd"

	# Remove repeated entries and their sequences from the original FASTA file
		awk -v dpdFile=specialCases/"\${name}.dpd" '
		BEGIN {
		while (getline < dpdFile) {
			repeated[\$0] = 1
			}
		}
		/^>/ { header = substr(\$0, 2)
		if (repeated[header]) {
			skip = 1
		} else {
			skip = 0}
		}
		!skip' "\$file" > specialCases/"\${name}_cleaned.fasta"

		for indexSeqs in specialCases/"\${name}Seqs_"*; do
			geneName=\$(basename "\${indexSeqs%Seqs_*}")
			sampleName=\$(basename "\${indexSeqs##*Seqs_}")

			# Read the longest line based on letter count
			sequence=\$(awk '
			{gsub(/[^a-zA-Z]/, "", \$0); len=length(\$0)}
			len > max_length {max_length=len; longest=\$0}
			END {print longest}' "\$indexSeqs")

        # Finally add the selected sequence back to the cleaned original FASTA file with the header as well
			echo ">\${sampleName}" >> specialCases/"\${name}_cleaned.fasta"
			echo "\$sequence" >> specialCases/"\${name}_cleaned.fasta"
		done

	# Cleaning temporary files
	#	rm specialCases/"\${name}Seqs_"* specialCases/"\${name}.dpd"

	done


        # If missing user sample == true, then append it and fill it with n's (can't treat them as gaps because there is uncertainty)
        for file in specialCases/*_cleaned.fasta ; do

                sampleValue=\$(awk '/^>/ {print \$0}' "\$file" | wc -l)
                numberOfColumns=\$(awk 'NR==2 {print length}' "\$file")
                totalSamples=\$(wc -l < sampleNames.txt)

                if (( sampleValue < totalSamples)); then

                        while read -r strain; do
                                if ! grep -wq "\$strain" "\$file"; then
                                        echo ">\$strain" >> "\$file"
                                        fakeSeq=\$(printf '%*s' "\$((numberOfColumns))" | tr ' ' 'n')
                                        echo "\$fakeSeq" >> "\$file"
                                fi
                        done < sampleNames.txt
                fi
        done

	# Turns out these broken entries were also incomplete
        for i in specialCases/*_cleaned.fasta; do

                numberOfColumns=\$(awk 'NR==2 {print length}' "\$i")
                echo "\$numberOfColumns"

                awk -v numCols="\$numberOfColumns" '{
                        if (\$0 ~ /^>/) {
                                print
                        } else {
                        while ( length(\$0) < numCols) {
                                \$0 = \$0 "n"
                        }
                        print
                        }
                }' "\$i" > specialCases/tmp && mv specialCases/tmp "\${i}"
        done


	cp specialCases/*_cleaned.fasta filteredGenes/

	for file in filteredGenes/*_cleaned.fasta; do
		name=\$(basename "\${file%_cleaned.fasta}_AlnSeq.fasta")
		mv "\$file" filteredGenes/"\$name"
	done


	# Make INDEX file so we can sort every file in the same order before building MSA based on first file.
	echo -e "Creating INDEX file"
	firstFile=\$(ls -1 filteredGenes/ | awk 'NR==1 {print \$0}') && awk '/^>/ {print \$0}' filteredGenes/"\$firstFile" > filteredGenes/INDEX
	echo -e "Done\\n"

	for fasta in filteredGenes/*.fasta; do
    		echo "Sorting \$fasta file"
		name=\$(basename "\${fasta%_AlnSeq.fasta}")

		while read -r header; do
        		awk -v headerName="\$header" '
            		\$0 ~ headerName {
                		print \$0     # Print the matched header
                		getline      # Get the sequence line after the header and store it in \$0
                		print \$0     # Print the sequene
            		}
        		' "\$fasta" >> filteredGenes/"\${name}_sorted"
    		done < filteredGenes/INDEX

    		echo "Finished sorting \$fasta to \${name}_sorted"
	done

	# Check that this worked
	for i in filteredGenes/*_sorted; do
    		awk '/^>/ {print \$0}' "\$i" > filteredGenes/currentHeaders

    		# Compare with INDEX file. It should be the same but if not, point that out

    		if ! diff -q filteredGenes/INDEX filteredGenes/currentHeaders >/dev/null; then
        		echo "Headers in \$i differ from INDEX." >> filteredGenes/notSorted
        		mv "\$i" specialCases/
    		else
        		echo "Headers in \$i match INDEX." >> filteredGenes/sortedSuccesfully
    		fi
	done

	if [ -f filteredGenes/notSorted ]; then
    		echo "There is one or more files with sorting problems. Check notSorted file and specialCases directory for more details"
	else
    		echo "It seems every file is sorted. Moving on"
	fi
	
	# Clean up temporary file
	#rm filteredGenes/currentHeaders


	echo -e "Check that every MSA contains the same amount of nucleotides"
	# (outgroup or user samples can generate insertions sometimes?)
	for gene in filteredGenes/*_sorted; do
		name=\$(basename "\$gene")
		echo "Reading \$gene file"
		awk '!/^>/ {print length}' "\$gene" >> filteredGenes/"\${name%_sorted}"_testingIntegrity.txt
		echo -e "Done"
	done

	for file in filteredGenes/*testingIntegrity.txt; do
		name=\$(basename "\$file")
	    	if [[ -f "\$file" ]]; then
	        	# Get unique values in the first column
        		uniqueValues=\$(cut -d ' ' -f 1 "\$file" | sort | uniq)

        		# Count the number of unique values
        		uniqueCount=\$(echo "\$uniqueValues" | wc -l)

        		# Delete the file if there is only one unique value
        		if [[ \$uniqueCount -eq 1 ]]; then
	            		rm "\$file"
            			echo "Deleted \$file (only 1 unique value)"
        		else
	            		echo "Kept \$file (more than 1 unique value, means problem)"
				mv "\$file" filteredGenes/"\${name%_testingIntegrity.txt}_problematicFile.txt"
        		fi
    		fi
	done

	for problem in filteredGenes/*_problematicFile.txt; do
		name=\$(basename "\${problem%_problematicFile.txt}")
		mv filteredGenes/"\${name}"_sorted specialCases/
	done

	for file in filteredGenes/*_sorted; do
		name=\$(basename "\$file")
		mv "\$file" filteredGenes/"\${name%_sorted}_Filtered.fasta"
	done
	
	cat .command.out >> filterGeneAlignments.log
	"""
}


/*   Add individual gene sequences to these particular gene.aln files
 *   To do this we need: If there is one or more samples missing in any particular gene.aln file, then we get the gene lenght and we add Ns.
 *   At the end , we will have every gene.aln file filled with every sample. So we can now concatenate every gene.aln file based on INDEX.
 */



process makeMSA {

	input:
	path genesAlnSeq, stageAs: 'genes/*'
	path maskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path maskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path maskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path genesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path sampleNames, stageAs: 'sampleNames.txt'

	output:
	path 'genesAbovePercentMSA.fasta', emit: genesAbovePercentMSA
	path 'maskedMatrixGenesNoUbiquitousMSA.fasta', emit: maskedMatrixGenesNoUbiquitousMSA
	path 'maskedMatrixGenesOnlyAncientMSA.fasta', emit: maskedMatrixGenesOnlyAncientMSA
	path 'maskedMatrixGenesUbiquitousMSA.fasta', emit: maskedMatrixGenesUbiquitousMSA
	

	script:
	"""
	#!/bin/bash

	touch genesAbovePercentMSA.fasta
	touch maskedMatrixGenesUbiquitousMSA.fasta
	touch maskedMatrixGenesOnlyAncientMSA.fasta
	touch maskedMatrixGenesNoUbiquitousMSA.fasta

	sed -i -e 's/~/_/g' genesAbovePercentSeries.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesNoUbiquitous.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesOnlyAncient.txt
	sed -i -e 's/~/_/g' maskedMatrixGenesUbiquitous.txt



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste genesAbovePercentMSA.fasta "\$file" > TMP; mv TMP genesAbovePercentMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < genesAbovePercentSeries.txt
	sed -i -e 's/\t//g' genesAbovePercentMSA.fasta
	awk -F'>' '/^>/ {print ">" \$2} !/^>/' genesAbovePercentMSA.fasta > TMPg; mv TMPg genesAbovePercentMSA.fasta



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesNoUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesNoUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesNoUbiquitous.txt
        sed -i -e 's/\t//g' maskedMatrixGenesNoUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesNoUbiquitousMSA.fasta > TMPnU; mv TMPnU maskedMatrixGenesNoUbiquitousMSA.fasta



	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesOnlyAncientMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesOnlyAncientMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesOnlyAncient.txt
        sed -i -e 's/\t//g' maskedMatrixGenesOnlyAncientMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesOnlyAncientMSA.fasta > TMPo2; mv TMPo2 maskedMatrixGenesOnlyAncientMSA.fasta

	while read -r gene; do
		file="genes/\${gene}_Filtered.fasta"
		if [[ -f "\$file" ]] ; then
			paste maskedMatrixGenesUbiquitousMSA.fasta "\$file" > TMP; mv TMP maskedMatrixGenesUbiquitousMSA.fasta
		else
			echo "Warning: File \$file not found, skipping..."
		fi
	done < maskedMatrixGenesUbiquitous.txt

        sed -i -e 's/\t//g' maskedMatrixGenesUbiquitousMSA.fasta
        awk -F'>' '/^>/ {print ">" \$2} !/^>/' maskedMatrixGenesUbiquitousMSA.fasta > TMPU2; mv TMPU2 maskedMatrixGenesUbiquitousMSA.fasta

	"""
}


process pMauve {
	conda "${projectDir}/envs/pMauve.yaml"

	input:
	path gffFiles, stageAs: '*'

	output:
	path 'pMauveAlignment.xmfa', emit: pMauveAlignment
	path 'pMauveAlignment.bbcols', emit: pMauveBbcols
	path 'pMauveAlignment.backbone', emit: pMauveBackbone
	path 'pMauveAlignmentCoreGenome', emit: pMauveCoreGenome

	script:
	"""
	progressiveMauve  --output=pMauveAlignment *
	mv pMauveAlignment ./pMauveAlignment.xmfa
	stripSubsetLCBs pMauveAlignment.xmfa pMauveAlignment.bbcols pMauveAlignmentCoreGenome 500
	"""
}

process xmfaToFasta {
	conda "${projectDir}/envs/biopython.yaml"

	input:
	path coreGenome, stageAs: 'pMauveAlignmentCoreGenome.xmfa'

	output:
	path 'pMauveFastaMSA.fasta', emit: pMauveFastaMSA

	script:
	"""
	convertXmfaToFasta.py pMauveAlignmentCoreGenome.xmfa
	"""
}

process filterMauveFasta {
	conda "${projectDir}/envs/seqtk.yaml"

	input:
	path mauveFastaMSA, stageAs: 'pMauveFastaMSA.fasta'

	output:
	path 'concatenatedSeqtkMauveFastaMSA.fasta', emit: concatenatedSeqtkMauveFastaMSA

	script:
	"""
	seqtk seq pMauveFastaMSA.fasta > seqtkMauveFastaMSA.fasta
	sed -i -e '/^>/ s/\\/.*//' seqtkMauveFastaMSA.fasta 

	awk '
	{
	    if (/^>/) {
	        header = \$1  
	        gsub(/^>/, "", header)  
	    } else {
	        seq[header] = seq[header] \$0  
	    }
	}
	END {
	    for (id in seq) {
	        print ">" id  
	        print seq[id] 
	    }
	}' seqtkMauveFastaMSA.fasta > concatenatedSeqtkMauveFastaMSA.fasta
	
	"""
}

process startingTree {
	conda "${projectDir}/envs/iqtree.yaml"	

	input:
	path concatenatedSeqtkMauveFastaMSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'

	output:
	path 'startingTreeMauveFasta.treefile', emit: startingTreeMauveFasta
	path 'startingTreeMauveFasta.iqtree', emit: startingTreeMauveFastaLog
	path 'kappaValue', emit: kappa

	script:
	"""
	iqtree -s concatenatedSeqtkMauveFastaMSA.fasta --prefix startingTreeMauveFasta -T 10 -B 1000 -m MFP
	
	awk -F':' '
		/A-C:/ {ACtransversion=\$2 + 0}
		/A-G:/ {AGtransition=\$2 + 0}
		/A-T:/ {ATtransversion=\$2 + 0}
		/C-T:/ {CTtransition=\$2 + 0}
		/C-G:/ {CGtransversion=\$2 + 0}
		/G-T:/ {GTtransversion=\$2 + 0}
		END {
		transitionRate= ((AGtransition + CTtransition)/2)
		transversionRate= (( ACtransversion + ATtransversion + CGtransversion + GTtransversion ) / 4)
		kappa = (transitionRate / transversionRate  )		

		print kappa }' startingTreeMauveFasta.iqtree > kappaValue

	"""
}


process findRecombinationSpots {
	conda "${projectDir}/envs/clonalframe.yaml" 

	input:
	path MSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'
	path startingTree, stageAs: 'startingTreeMauveFasta.treefile'
	path kappa, stageAs: 'kappaValue'
	

	output:
	path 'recombinantOutputs.importation_status.txt', emit: recombinationMap
	
	script:
	"""
	kappa=\$(cat kappaValue)
	ClonalFrameML startingTreeMauveFasta.treefile concatenatedSeqtkMauveFastaMSA.fasta recombinantOutputs -kappa "\${kappa}"

	"""
}

process mapRecombinantsToGenes {
	conda "${projectDir}/envs/blast.yaml"

	input:
	path recombinationMap, stageAs: 'recombinantOutputs.importation_status.txt'
	path MSA, stageAs: 'concatenatedSeqtkMauveFastaMSA.fasta'
	path db, stageAs: 'database/*'
	path GFF, stageAs: 'gff/*'

	output:
	stdout

	script:
	"""
	#!/bin/bash

	awk '!/^NODE/ && NR>1 {print \$0}' recombinantOutputs.importation_status.txt > filteredRecombinationMap

	while read -r name beg end; do
	    awk -v name="\$name" -v start="\$beg" -v end="\$end" '
	    /^>/ {
	        # Process headers
	        headerSequence = substr(\$0, 2)  # remove the ">" thing to get the header name
	        isHeader = (headerSequence == name)  # check if the header matches the target name
	
	        if (isHeader) {
	            fullSeq = ""  # Reset fullSeq for the new header
	            currentHeader = name "_seq"  # make unique header for this sequence
	        }
	    }
	    !/^>/ && isHeader {
	        fullSeq = fullSeq \$0   
	    }
	    /^>/ && fullSeq != "" {
	        # store fullSeq for each unique header when a new header starts
	        # remove sequences shorter than 30 bp? To avoid mapping uncertainty
	        seqTesting = substr(fullSeq, start, end - start + 1)
	        if (length(seqTesting) >= 30) {
	            seqArray[currentHeader] = seqTesting
        	}
	        fullSeq = ""
	    }
	    END {
	        # this is for the last entry
	        if (fullSeq != "") {
	            seqTesting = substr(fullSeq, start, end - start + 1)
	            if (length(seqTesting) >= 30) {
	                seqArray[currentHeader] = seqTesting
	            }
	        }
	        
	        # print the sequences
	        for (header in seqArray) {
	            printf(">%s\\n%s\\n", header, seqArray[header])
	        }
	    }
	    ' concatenatedSeqtkMauveFastaMSA.fasta >> "\${name%.fna}"_TMP.fasta
	done < filteredRecombinationMap
	
	for i in *_TMP.fasta; do
		name=\$(basename "\${i%_TMP.fasta}")
		awk '
		    BEGIN { count = 0 }  
		    /^>/ { 
		        count++  
		        \$0 = \$0 "_" count  # this just add a counter ID for each header to make them unique
		        print
		    } 
		    !/^>/ { 
		        print  # print the sequence line only
		    }
		' "\$i"  > "\$name"_newfileWithExtractedHeadersAndSequences.fasta
	done

	rm *TMP.fasta
	
	# Mapping sequences to PanGenomeReference

	for sample in *newfileWithExtractedHeadersAndSequences.fasta; do
		name=\$(basename "\${sample%_newfileWithExtractedHeadersAndSequences.fasta}")
		blastn -query "\$sample" -db database/panGenomeReferenceDB -out blastResults"\$name" -outfmt 6
	done	


	echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tseqLength" > blastSummaryResults.tab

	for sample in *newfileWithExtractedHeadersAndSequences.fasta; do
    		awk '/^>/ {name = \$0 ; getline ; seqLength = length(\$0); print name, seqLength}' "\$sample" >> recombinantsIDplusLengths.txt
	done

	sed -i -e 's/>//g' recombinantsIDplusLengths.txt
	
	cat blastResults* >> TMP1.tab
	sed -i -e 's/~/_/g' TMP1.tab

	while read -r sample seqLength; do
    		grep -w "\$sample" TMP1.tab | awk -v value="\$seqLength" '{print \$0, value}' >> blastSummaryResults.tab
	done < recombinantsIDplusLengths.txt

	sed -i -e 's/ /\\t/g' blastSummaryResults.tab

	rm TMP1.tab
	"""
}



process treeThreshold {
	conda "${projectDir}/envs/iqtree.yaml"
	
	input:
	path genesMSA, stageAs: 'genesAbovePercentMSA.fasta'

	output:
	path 'genesAbovePercentMSA.contree', emit: genesAbovePercentMSAContree
	path 'genesAbovePercentMSA.iqtree', emit: genesAbovePercentMSAIqtree
	path 'genesAbovePercentMSA.log', emit: genesAbovePercentMSALog
	path 'genesAbovePercentMSA.treefile', emit: genesAbovePercentMSATreefile

	script:
	"""
	iqtree -s genesAbovePercentMSA.fasta --prefix genesAbovePercentMSA -T 10 -B 1000 -m MFP 
	"""
}

process treeUbiquitous {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ubiquitousMSA, stageAs: 'maskedMatrixGenesUbiquitousMSA.fasta'

        output:
        path 'maskedMatrixGenesUbiquitousMSA.contree', emit: maskedMatrixGenesUbiquitousMSAContree
	path 'maskedMatrixGenesUbiquitousMSA.iqtree', emit: maskedMatrixGenesUbiquitousMSAIqtree
	path 'maskedMatrixGenesUbiquitousMSA.log', emit: maskedMatrixGenesUbiquitousMSALog
	path 'maskedMatrixGenesUbiquitousMSA.treefile', emit: maskedMatrixGenesUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesUbiquitousMSA.fasta --prefix maskedMatrixGenesUbiquitousMSA -T 10 -B 1000 -m MFP
        """
}

process treeNoUbiquitous {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path noUbiquitousMSA, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.fasta'

        output:
        path 'maskedMatrixGenesNoUbiquitousMSA.contree', emit: maskedMatrixGenesNoUbiquitousMSAContree
	path 'maskedMatrixGenesNoUbiquitousMSA.iqtree', emit: maskedMatrixGenesNoUbiquitousMSAIqtree
	path 'maskedMatrixGenesNoUbiquitousMSA.log', emit: maskedMatrixGenesNoUbiquitousMSALog
	path 'maskedMatrixGenesNoUbiquitousMSA.treefile', emit: maskedMatrixGenesNoUbiquitousMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesNoUbiquitousMSA.fasta --prefix maskedMatrixGenesNoUbiquitousMSA -T 10 -B 1000 -m MFP
        """
}


process treeAncient {
        conda "${projectDir}/envs/iqtree.yaml"
        
        input:
        path ancientMSA, stageAs: 'maskedMatrixGenesOnlyAncientMSA.fasta'

        output:
        path 'maskedMatrixGenesOnlyAncientMSA.contree', emit: maskedMatrixGenesOnlyAncientMSAConqtree
	path 'maskedMatrixGenesOnlyAncientMSA.iqtree', emit: maskedMatrixGenesOnlyAncientMSAIqtree
	path 'maskedMatrixGenesOnlyAncientMSA.log', emit: maskedMatrixGenesOnlyAncientMSALog
	path 'maskedMatrixGenesOnlyAncientMSA.treefile', emit: maskedMatrixGenesOnlyAncientMSATreefile

        script:
        """
        iqtree -s maskedMatrixGenesOnlyAncientMSA.fasta --prefix maskedMatrixGenesOnlyAncientMSA -T 10 -B 1000 -m MFP
        """
}

process outgroupEntrez {
	conda "${projectDir}/envs/entrez.yaml"

	input:
	val outID
	
	output:
	path '*fna', emit: outgroupFasta

	script:
	"""
	#!/bin/bash
	counter=0
	esearch -db assembly -query "txid${outID}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | while read url; do
        
	if [ "\$counter" -ge 1 ]; then
		break
	fi

	if [ -z "\$url" ]; then
		continue
	fi
        
	fname="\$(basename "\$url")"
	wget "\$url/\${fname}_genomic.fna.gz"

	counter="\$((counter + 1))"	
	
	gunzip -f "\${fname}_genomic.fna.gz"
	
	done	
	"""
}


process makeReads {
	conda "${projectDir}/envs/art.yaml"

	input:
	path outgroupFasta, stageAs: 'outgroupFasta.fasta'
	
	output:
	path 'outgroupReads.fq', emit: outgroupReads


	script:
	"""
	art_illumina -i outgroupFasta.fasta -l 150 -f 100 -o outgroupReads
	"""
}


process outgroupAlignmentFAndiltering {
	conda "${projectDir}/envs/alignment.yaml"
	
	input:
	path outgroupReads, stageAs: 'outgroupReads.fq'
	path panGenomeRef, stageAs: 'panGenomeReferenceSeq.fasta'
	val threadsGlobal

	output:
	path 'outgroupFastaPostAlignment.bam', emit: outgroupFastaPostAlignment

	script:
	"""
	bwa index panGenomeReferenceSeq.fasta
	bwa mem -B 1 -E 1 panGenomeReferenceSeq.fasta outgroupReads.fq -t $threadsGlobal > outgroupFasta.sam
	samtools view -bS outgroupFasta.sam > outgroupFasta.bam
	samtools quickcheck outgroupFasta.bam
	samtools sort -o outgroupFastaSorted.bam -O bam -@ $threadsGlobal outgroupFasta.bam
	samtools index outgroupFastaSorted.bam
	samtools view -b -@ 10 -F 4 outgroupFastaSorted.bam > outgroupFastaSortedMappedreads.bam
	samtools index outgroupFastaSortedMappedreads.bam
	samtools sort -o outgroupFastaPostAlignment.bam -O bam -@ $threadsGlobal outgroupFastaSortedMappedreads.bam
	"""
}


process makeOutgroupConsensus {
	conda "${projectDir}/envs/consensus.yaml"

	input:
	path outgroupFastaPostAlignment, stageAs: 'outgroupFastaPostAlignment.bam'
	path panGenomeRef, stageAs: 'panGenomeRef.fasta'

	output:
	path 'extractedSequencesOutgroup.fasta', emit: extractedSequencesOutgroupFasta
	path 'extractedSequencesOutgroup.fq', emit: extractedSequencesOutgroupFastq

	script:
	"""
	bcftools mpileup -f panGenomeRef.fasta outgroupFastaPostAlignment.bam | bcftools call -c | vcfutils.pl vcf2fq > extractedSequencesOutgroup.fq
	seqtk seq -a extractedSequencesOutgroup.fq > extractedSequencesOutgroup.fasta

	"""
}


process blastMe {
	conda "${projectDir}/envs/blast.yaml"

	input:
	path panSeq, stageAs: 'panGenomeReference.fasta'

	output:
	path 'panGenomeReferenceDB*', emit: panGenomeReferenceDB

	script:
	"""
	makeblastdb -in panGenomeReference.fasta -dbtype nucl -out panGenomeReferenceDB
	"""
}


process getResults {

	input:
	path makeDir
	path checkedFastas, stageAs: 'checkedFasta/*'
	path checkedGffs, stageAs: 'checkedGff/*'
	path fastaDatabaseLog, stageAs: 'fastaDatabase.log'
	path fastaDatabaseSeqs, stageAs: 'clusteredSequences.fasta'
	path clusteredSequences, stageAs: 'clusteredNonRedundantGenes.fasta'
	path clusteringLogFile, stageAs: 'clustering.log'
	path prokkaGff, stageAs: 'prokkaGff/*'
	path prokkaLog, stageAs: 'prokka.log'
	path panarooLog, stageAs: 'makePangenome.log'
	path genesMSA, stageAs: 'geneMSA/*'
	path panrefG, stageAs: 'pangenomeReferenceGenome.fasta'
	path geneNormalizedUpdated, stageAs: 'geneNormalizedUpdated.tab'
	path globalMeanCoverage, stageAs: 'globalMeanCoverage.tab'
	path PostAlignmentFiles, stageAs: '*'
	path RefLenghts, stageAs: '*'
	path RawCoverage, stageAs: '*'
	path CompletenessSummary, stageAs: 'completenessSummary.tab'
	path FinalMatrix, stageAs: 'finalMatrix.tab'
	path presenceAbsenceplots, stageAs: '*'
	path MaskedMatrixGenesOnlyAncient, stageAs: 'maskedMatrixGenesOnlyAncient.txt'
	path MaskedMatrixGenesUbiquitous, stageAs: 'maskedMatrixGenesUbiquitous.txt'
	path MaskedMatrixGenesNoUbiquitous, stageAs: 'maskedMatrixGenesNoUbiquitous.txt'
	path GenesAbovePercentSeries, stageAs: 'genesAbovePercentSeries.txt'
	path GenesAbovePercentMSAIqtree, stageAs: 'genesAbovePercentMSA.iqtree'
	path GenesAbovePercentMSALog, stageAs: 'genesAbovePercentMSA.log'
	path GenesAbovePercentMSATreefile, stageAs: 'genesAbovePercentMSA.treefile'
	path MaskedMatrixGenesUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesUbiquitousMSA.iqtree'
	path MaskedMatrixGenesUbiquitousMSALog, stageAs: 'maskedMatrixGenesUbiquitousMSA.log'
	path MaskedMatrixGenesUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesUbiquitousMSA.treefile'
	path MaskedMatrixGenesNoUbiquitousMSAIqtree, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.iqtree'
	path MaskedMatrixGenesNoUbiquitousMSALog, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.log'
	path MaskedMatrixGenesNoUbiquitousMSATreefile, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.treefile'
	path MaskedMatrixGenesOnlyAncientMSAIqtree, stageAs: 'maskedMatrixGenesOnlyAncientMSA.iqtree'
	path MaskedMatrixGenesOnlyAncientMSALog, stageAs: 'maskedMatrixGenesOnlyAncientMSA.log'
	path MaskedMatrixGenesOnlyAncientMSATreefile, stageAs: 'maskedMatrixGenesOnlyAncientMSA.treefile'

	output:
	stdout

	script:
	"""
	#!/bin/bash

	mkdir -p "${makeDir}/modernData/"

	mv checkedFasta/* "${makeDir}/modernData/"
	mv checkedGff/* "${makeDir}/modernData/"
	mv fastaDatabase.log "${makeDir}/modernData/"

	mkdir -p "${makeDir}/clusteredSequences"
	mv clusteredSequences.fasta "${makeDir}/clusteredSequences/"
	mv clusteredNonRedundantGenes.fasta "${makeDir}/clusteredSequences/"
	mv clustering.log "${makeDir}/clusteredSequences/"

	mkdir -p "${makeDir}/prokkaResults/"

	mv prokkaGff/* "${makeDir}/prokkaResults/"
	mv prokka.log "${makeDir}/prokkaResults/"

	mkdir -p "${makeDir}/pangenomeFiles"

	mv makePangenome.log "${makeDir}/pangenomeFiles/"
	mv geneMSA/* "${makeDir}/pangenomeFiles/"
	mv pangenomeReferenceGenome.fasta "${makeDir}/pangenomeFiles/"

	mkdir -p "${makeDir}/alignmentResults"

	mv geneNormalizedUpdated.tab "${makeDir}/alignmentResults/"
	mv globalMeanCoverage.tab "${makeDir}/alignmentResults/"
	mv *_refLength.txt "${makeDir}/alignmentResults/"
	mv *Coverage.txt "${makeDir}/alignmentResults/"
	mv *bam "${makeDir}/alignmentResults/"
	mv completenessSummary.tab "${makeDir}/alignmentResults/"

	mkdir -p "${makeDir}/matrixResults"
	mv finalMatrix.tab "${makeDir}/matrixResults"

	mkdir -p "${makeDir}/plotsResults"
	mv *png "${makeDir}/plotsResults"

	mkdir -p "${makeDir}/MSAs"
	mv maskedMatrixGenesUbiquitous.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesOnlyAncient.txt "${makeDir}/MSAs/"
	mv maskedMatrixGenesNoUbiquitous.txt "${makeDir}/MSAs/"
	mv genesAbovePercentSeries.txt "${makeDir}/MSAs/"
	mv *MSA* "${makeDir}/MSAs/"

	"""
}

workflow {
	if (!params.trustedGenomes) {
		entrez(downloadGenomes, taxID, resultsDir)
		fastaFiles = entrez.out.fastaFiles
		gffFiles = entrez.out.gffFiles

	} else {
		trustedData(trustedDataChannel)
		fastaFiles = trustedData.out.trustedFasta
		gffFiles = trustedData.out.trustedGenBank

	}


	fastaDatabase(gffFiles, fastaFiles)
	clustering(fastaDatabase.out.theFastaDatabase, cdHitCluster, threadsGlobal)
	prokkaMakeAnnotations(clustering.out.clusteredDatabase, threadsGlobal, fastaDatabase.out.validGff, fastaDatabase.out.validFasta)
	makePangenome(prokkaMakeAnnotations.out.prokkaGFF, pangenomeMode, pangenomeThreshold, threadsGlobal)
	formattingPangenome(makePangenome.out.panSequence)
	blastMe(formattingPangenome.out.panGenomeReference)
        outgroupEntrez(outTax)
        makeReads(outgroupEntrez.out.outgroupFasta)
        outgroupAlignmentFAndiltering(makeReads.out.outgroupReads, formattingPangenome.out.panGenomeReference, threadsGlobal)
        makeOutgroupConsensus(outgroupAlignmentFAndiltering.out.outgroupFastaPostAlignment, formattingPangenome.out.panGenomeReference)
	alignment(reads, formattingPangenome.out.panGenomeReference, threadsGlobal, configFile, missingProb, seedAlignment, gapFraction, minReadLength, maxReadLength)
	alignmentSummary(configFile, alignment.out.postAlignedBams)
	normalizationFunction(alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage)
	updateNormalization(normalizationFunction.out.geneNormalizedSummary, alignmentSummary.out.completenessSummary)


	if (params.genotyper == "gatk") {
		gatkConsensus(formattingPangenome.out.panGenomeReference, alignmentSummary.out.postAlignmentFiles, formattingPangenome.out.panGenomeReferenceDictionary, formattingPangenome.out.panGenomeReferenceIndex)
		extractedSequencesFasta = gatkConsensus.out.gatkConsensusSequences
		vcfFile = gatkConsensus.out.gatkGenotypes

	} else if (params.genotyper == "bcftools") {
		bcftoolsConsensus(formattingPangenome.out.panGenomeReference, alignmentSummary.out.postAlignmentFiles)
		extractedSequencesFasta = bcftoolsConsensus.out.consensusSequences

	} else {
		error "Invalid option for --genotyper. Please choose 'gatk' or 'bcftools'."
	}

	plotCoveragevsCompleteness(updateNormalization.out.geneNormalizedUpdated, geneCompleteness, normalizedCoverageDown)
        applyCoverageBounds(updateNormalization.out.geneNormalizedUpdated, normalizedCoverageDown, normalizedCoverageUp, geneCompleteness)
	makeMatrix(makePangenome.out.initialMatrix , normalizationFunction.out.globalMeanCoverage, applyCoverageBounds.out.geneNormalizedUpdatedFiltered)
	buildHeatmap(makeMatrix.out.finalCsv, makeMatrix.out.INDEX ,makeMatrix.out.matrix, makeMatrix.out.sampleNames)
	plotCoveragevsCompletenessOnFiltered(applyCoverageBounds.out.geneNormalizedUpdatedFiltered, geneCompleteness,normalizedCoverageDown)
	filterGeneAlignments(makePangenome.out.alignedGenesSeqs, extractedSequencesFasta, fastaDatabase.out.validFasta, downloadGenomes, makeOutgroupConsensus.out.extractedSequencesOutgroupFasta, buildHeatmap.out.blackListed)
	pMauve(fastaDatabase.out.validFasta)
	makeMSA(filterGeneAlignments.out.genesAlnSeq, buildHeatmap.out.maskedMatrixGenesNoUbiquitous, buildHeatmap.out.maskedMatrixGenesOnlyAncient, buildHeatmap.out.maskedMatrixGenesUbiquitous, buildHeatmap.out.genesAbovePercentSeries, filterGeneAlignments.out.sampleNames)
	treeThreshold(makeMSA.out.genesAbovePercentMSA)
	treeUbiquitous(makeMSA.out.maskedMatrixGenesUbiquitousMSA)
	treeNoUbiquitous(makeMSA.out.maskedMatrixGenesNoUbiquitousMSA)
	treeAncient(makeMSA.out.maskedMatrixGenesOnlyAncientMSA)
	xmfaToFasta(pMauve.out.pMauveCoreGenome)
	filterMauveFasta(xmfaToFasta.out.pMauveFastaMSA)
	startingTree(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA)
	findRecombinationSpots(filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, startingTree.out.startingTreeMauveFasta, startingTree.out.kappa)
	mapRecombinantsToGenes(findRecombinationSpots.out.recombinationMap, filterMauveFasta.out.concatenatedSeqtkMauveFastaMSA, blastMe.out.panGenomeReferenceDB, prokkaMakeAnnotations.out.prokkaGFF)
	getResults(
	resultsDir, fastaDatabase.out.validFasta , fastaDatabase.out.validGff , fastaDatabase.out.fastaDatabaseLogFile , fastaDatabase.out.theFastaDatabase, 
	clustering.out.clusteredDatabase, clustering.out.clusteringLog, prokkaMakeAnnotations.out.prokkaGFF, prokkaMakeAnnotations.out.prokkaLogfile,  makePangenome.out.panarooLog,
	filterGeneAlignments.out.genesAlnSeq, formattingPangenome.out.panGenomeReference, updateNormalization.out.geneNormalizedUpdated, normalizationFunction.out.globalMeanCoverage,
	alignmentSummary.out.postAlignmentFiles, alignmentSummary.out.refLenght, alignmentSummary.out.rawCoverage, alignmentSummary.out.completenessSummary, buildHeatmap.out.finalMatrix,
	buildHeatmap.out.presenceAbsence, buildHeatmap.out.maskedMatrixGenesOnlyAncient, buildHeatmap.out.maskedMatrixGenesUbiquitous, buildHeatmap.out.maskedMatrixGenesNoUbiquitous,
	buildHeatmap.out.genesAbovePercentSeries, treeThreshold.out.genesAbovePercentMSAIqtree ,treeThreshold.out.genesAbovePercentMSALog , treeThreshold.out.genesAbovePercentMSATreefile,
	treeUbiquitous.out.maskedMatrixGenesUbiquitousMSAIqtree, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSALog, treeUbiquitous.out.maskedMatrixGenesUbiquitousMSATreefile, 
	treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSAIqtree , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSALog , treeNoUbiquitous.out.maskedMatrixGenesNoUbiquitousMSATreefile,
	treeAncient.out.maskedMatrixGenesOnlyAncientMSAIqtree , treeAncient.out.maskedMatrixGenesOnlyAncientMSALog , treeAncient.out.maskedMatrixGenesOnlyAncientMSATreefile
	)
}
