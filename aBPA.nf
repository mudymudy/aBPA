#!/usr/bin/env nextflow

// Default parameter values
params.data = ""
params.output = ""
params.lcompleteness = 50
params.coverage = 0.5
params.softclipping = 5
params.threads = 10
params.tax_id = null
params.genomes = 100
params.clustering = 0.95
params.core = 0.01
params.clean = "strict"

// Enable DSL2
nextflow.enable.dsl=2

// Help function
process print_help {
    	when:
    	params.help
	
    	script:
    	"""
    	RED='\\033[1;31m'
    	YELLOW='\\033[1;33m'
    	NC='\\033[0m' # No Color
	
    	echo -e "\${RED}SYNOPSIS\${NC}"
    	echo -e "\\n\${YELLOW}USAGE\${NC} \\n\\nnextflow run aBPA.nf -params-file params.json"
    	echo -e "\\n\${YELLOW}OPTIONS\${NC}"
    	echo -e "\\nMandatory"
    	echo -e "  --data <PATH>            Set data file path"
    	echo -e "  --output <PATH>          Set output directory"
    	echo -e "  --taxid <INT>            Set taxonomic ID value"
	
    	echo -e "\\nOptional"
    	echo -e "  --threads <INT>          Set number of threads (default: 10)"
    	echo -e "  --lcompleteness <INT>     Set gene completeness threshold (default: 50)"
    	echo -e "  --coverage <FLOAT>       Set mean depth of coverage threshold (default: 0.5)"
    	echo -e "  --soft-clipping <INT>    Set soft-clipping value (default: 5)"
    	echo -e "  --genomes <INT>          Set number of genomes to download (default: 100)"
    	echo -e "  --clustering <FLOAT>     Set clustering threshold (default: 0.95)"
    	echo -e "  --core-threshold <FLOAT> Set core genome threshold (default: 0.01)"
    	echo -e "  --clean-mode <STRING>    Set pangenome clean mode (default: strict)"
    	echo -e "  --help                  Print this help message and exit"
	    
    	echo -e "\\n\${RED}DESCRIPTION\${NC}"
    	echo -e "\\n\${YELLOW}--data <PATH>\${NC}: Specify the full path of your data (e.g., /home/user/data)"
    	echo -e "\\n\${YELLOW}--output <PATH>\${NC}: Specify the full path of your output directory. Create it before running the program."
    	exit 0
    	"""
}

// Process for downloading GenBank and FASTA files based on taxonomic ID
process download_genbank_fasta {
	// Declare conda environment
	conda 'envs/entrez.yaml'
    
    	input:
    	val tax_id from params.tax_id
    	val genomes from params.genomes
    	val output_dir from params.output

    	script:
    	"""
    	echo -e "Setting up directory structure\n"

    	mkdir -p "${output_dir}/NCBI/FASTA"
    	mkdir -p "${output_dir}/NCBI/GFF"
    	mkdir -p "${output_dir}/CLUSTERING"
    	mkdir -p "${output_dir}/PROKKA/GFF"
    	mkdir -p "${output_dir}/PANGENOME"
    	mkdir -p "${output_dir}/ALIGNMENTS"
	mkdir -p "${output_dir}/NORMALIZATION"
    	mkdir -p "${output_dir}/MATRIX"
    	mkdir -p "${output_dir}/PLOTS"
    	mkdir -p "${output_dir}/HETEROPLASMY/intermediate_files"
    	mkdir -p "${output_dir}/HETEROPLASMY/distributions"

    	echo -e "Done\n"

    	echo -e "Downloading GenBank and FASTA files based on taxonomic ID\n"
    	counter=0

    	esearch -db assembly -query "txid${tax_id}[Organism]" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | while read -r url; do
      	if [ "\$counter" -ge "${genomes}" ]; then
	        break
      	fi
      
      	fname=\$(basename "\$url")
      
      	wget -P "${output_dir}/NCBI/GFF" "\${url}/\${fname}_genomic.gbff.gz"
      	wget -P "${output_dir}/NCBI/FASTA" "\${url}/\${fname}_genomic.fna.gz"
      
      	counter=\$((counter + 1))
      	echo -e "Strains downloaded: \$counter of ${genomes}"
    	done

    	echo -e "Done\n"
    	"""
}

process parse_and_build_fasta_db {
	// Declare conda environment
	conda 'envs/biopython.yaml'
    
	input:
    	val output_dir from params.output

	script:
	"""
	echo -e "Parsing and building FASTA database"

	gzip -d "${output_dir}/NCBI/GFF/*"
	gzip -d "${output_dir}/NCBI/FASTA/*"
	python scripts/parsing_and_contatenating.py "${output_dir}/NCBI/GFF"

	mv clustered_sequences.fasta "${output_dir}/CLUSTERING/"

	echo -e "Done\n"
	"""
}

process clustering_seqs {
	// Declare conda environment
	conda 'envs/cdhit.yaml'

	input:
	val output_dir from params.output
	val clustering from params.clustering
	val threads from params.threads

	script:
	"""
	echo -e "Clustering gene sequences"

	cd-hit-est -i "$output"/CLUSTERING/clustered_sequences.fasta -o "$output"/CLUSTERING/clustered_non_redundant_genes.fasta -c "$clustering" -n "$threads"

	echo -e "Done\n"
	"""
}



process prokka {
	// Declare conda environment
	conda 'envs/prokka.yaml'

	input:
	val output_dir from params.output
	val threads from params.threads
	
	script:
	"""
	echo -e "Annotating FASTA sequences\n"

	ls -l "$output"/NCBI/GFF/*gbff | awk 'NR==1{print $NF}' > "$output"/NCBI/GFF/first
	name=$(cat "$output"/NCBI/GFF/first | awk -F'/' '{print $NF}')
	echo -e "$name"
	species=$(head -n 20 "$output"/NCBI/GFF/"$name" | grep "ORGANISM" | awk '{print $2, $3}' | sed -e 's/ /_/g')
	echo -e "$species"
	for i in "$output"/NCBI/FASTA/*; do
	        name=$(basename "$i")
	        prokka --outdir "$output"/PROKKA/"${name%.fna}" --addgenes  --addmrna --species "$species" --proteins "$output"/CLUSTERING/clustered_non_redundant_genes.fasta --force --cpus "$threads" "$i"
	done


	for sample in "$output"/PROKKA/*; do
	        name=$(basename "$sample")
	        mv "$sample"/*gff "$output"/PROKKA/GFF/"${name}.gff"
	done

	echo -e "Done\n"
	"""
}


process panaroo {
	// Declare conda environment
	conda 'envs/panaroo.yaml'

	input:
	val output_dir from params.output
	val threads from params.threads
	val clean from params.clean
	val core from params.core

	script:
	"""
	echo -e "Building pangenome from annotated files\n"

	panaroo -i "$output"/PROKKA/GFF/*.gff -o "$output"/PANGENOME/ --clean-mode "$clean" -a core --core_threshold "$core" -t "$threads"

	echo -e "Done\n"
	"""

}


process alignment {
	// Declare conda environment
	conda 'envs/alignment.yaml'
	input:
	val output_dir from params.output
	val threads from params.threads
	val softclipping from params.softclipping

	script:
	"""
	echo -e "Formatting pangenome sequence"

	seqtk seq "$output"/PANGENOME/pan_genome_reference.fa > "$output"/PANGENOME/pan_genome_reference.fasta

	echo -e "Done\n"

	echo -e "Aligning samples against pangenome reference"

	for sample in "$data"/*; do
		echo -e "\nGenerating the reference genome index files . . ."
		bwa index "$output"/PANGENOME/pan_genome_reference.fasta
	        name=$(basename "$sample")

		echo -e "\nRunning alignment against reference . . ."
		bwa aln -l 16500 -n 0.01 -o 2 -t "$threads" "$output"/PANGENOME/pan_genome_reference.fasta "$sample" > "$output"/ALIGNMENTS/"${name%.fastq*}.sai"
	
		echo -e "\nConverting SAI to SAM . . ."
		bwa samse "$output"/PANGENOME/pan_genome_reference.fasta "$output"/ALIGNMENTS/"${name%.fastq*}.sai" "$sample" > "$output"/ALIGNMENTS/"${name%.fastq*}.sam"
	
		echo -e "\nConverting SAM to BAM and sorting . . ."
		samtools view -bS "$output"/ALIGNMENTS/"${name%.fastq*}.sam" > "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
	
		echo -e "\nChecking sanity of BAM file . . ."
		samtools quickcheck "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
	
		echo -e "\nSorting BAM . . ."
		samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
	
		echo -e "\nGenerating BAM index . . ."
		samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam"
	
		echo -e "\nGetting only mapped reads . . . "
		samtools view -b -@ 10 -F 4 "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"
		samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"
	
		echo -e "\nRemoving 5 bases at each end of every reads . . ."
		~/miniforge3/envs/alignment/bin/bam trimBam "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam" -L "$softclipping" -R "$softclipping" --clip
	
		echo -e "\nGetting only reads with 25 mapping quality or more . . ."
		samtools view -q 25 -o "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam"
	
		echo -e "\nRemoving reads smaller than 34bp . . ."
		samtools view -e 'length(seq)>34' -O BAM -o "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam"
	
		echo -e "\nSorting the output . . ."
		samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam"
	
		echo -e "\nComputing basic statistics . . ."
		samtools coverage "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name}"_genomicsMetrics.txt
	
		echo -e "\nConverting to FASTQ . . ."
		samtools fastq -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_final.fastq"

		rm "$output"/ALIGNMENTS/*_softclipped.bam "$output"/ALIGNMENTS/*_qc.bam "$output"/ALIGNMENTS/*_lg.bam "$output"/ALIGNMENTS/*sai "$output"/ALIGNMENTS/*sam "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
	
	done

	echo -e "Done\n"
	"""
}

process raw_extracting {
	// Declare conda environment
	conda 'envs/alignment.yaml'

	input:
	val output_dir from params.output


	script:
	"""
	echo -e "Extracting raw coverage per gene\n"
	
	
	for i in "$output"/ALIGNMENTS/*_DMC_P.bam; do
	        name=$(basename "$i")
		samtools index "$i"
		samtools depth -a "$i" > "$output"/NORMALIZATION/"${name%_DMC_P.bam}_rawCoverage.txt"
		samtools idxstats "$i" | awk '{sum += $2} END {print sum}' > "$output"/NORMALIZATION/"${name%_DMC_P.bam}_refLength.txt"
		samplename=$(basename "${name%_DMC_P.bam}")
		samtools coverage "$i" | awk -v samplename="$samplename" 'NR>1 {print samplename, $1, $6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t $'\t' >> "$output"/NORMALIZATION/completenessSummary.tab

	done

	echo -e "Done\n"
	"""
}


process normalize_array {
	// Declare conda environment
	conda 'envs/normalization.yaml'

	input:
	val output_dir from params.output


	script:
	"""
	echo -e "Normalizing coverage per gene\n"
	
	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > "$output"/NORMALIZATION/geneNormalizedSummary.txt
	echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > "$output"/NORMALIZATION/globalMeanCoverage.txt
	
	for i in "$output"/NORMALIZATION/*_rawCoverage.txt; do
    	name=$(basename "${i%_rawCoverage.txt}")
    	echo -e "Sample being processed: $name\n"
	
    	# Compute global mean coverage
    	globalMean=$(awk -v name="$name" '{sum += $3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check " name}' "$i")
    	finalCount=$(awk -v name="$name" '{fcount++} END {print fcount}' "$i")
    	refCount=$(cat "$output"/NORMALIZATION/"${name}_refLength.txt")
    	echo -e "$name\t$finalCount\t$refCount\t$globalMean" >> "$output"/NORMALIZATION/globalMeanCoverage.txt
	
    	# Normalize coverage per gene
    	awk -v globalMean="$globalMean" -v name="$name" -v sampleCoverage="$finalCount" -v refCount="$refCount" '
        	{
            	if ($2 > geneLength[$1]) {
	                geneLength[$1] = $2
            	}
            	sumgene[$1] += $3
            	countgene[$1]++
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
    	' "$i" >> "$output"/NORMALIZATION/geneNormalizedSummary.txt
	
	done
	
	
	echo -e "Done\n"
	
	echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > "$output"/NORMALIZATION/geneNormalizedUpdated.tab
	
	sed -i -e 's/~/_/g' "$output"/NORMALIZATION/geneNormalizedSummary.txt
	
	awk 'NR>1{print $1"XYZ"$2, $3, $4, $5, $6}' "$output"/NORMALIZATION/geneNormalizedSummary.txt > "$output"/NORMALIZATION/TMP1
	
	awk '{print $1"XYZ"$2, $3}' "$output"/NORMALIZATION/completenessSummary.tab > "$output"/NORMALIZATION/TMP2
	
	while read -r ID completeness;do
	
		if grep -wq "${ID}" "$output"/NORMALIZATION/TMP1; then
			oldLine=$(grep -w "${ID}" "$output"/NORMALIZATION/TMP1)
			specificCompleteness=$(grep -w "${ID}" "$output"/NORMALIZATION/TMP2 | awk '{print $NF}')
			echo -e "${oldLine}\t${specificCompleteness}" >> "$output"/NORMALIZATION/geneNormalizedUpdated.tab
		fi
	
	done < "$output"/NORMALIZATION/TMP2
	"""
}


process normalization_and_plots {
	// Declare conda environment
	conda 'envs/normalization.yaml'

	input:
	val output_dir from params.output
	val lcompleteness from params.lcompleteness
	val coverage from params.coverage


	script:
	"""
	sed -i -e 's/XYZ/\t/g' "$output"/NORMALIZATION/geneNormalizedUpdated.tab
	sed -i -e 's/ /\t/g' "$output"/NORMALIZATION/geneNormalizedUpdated.tab
	
	
	rm "$output"/NORMALIZATION/TMP1 "$output"/NORMALIZATION/TMP2 "$output"/NORMALIZATION/geneNormalizedSummary.txt "$output"/NORMALIZATION/completenessSummary.tab
	mv "$output"/NORMALIZATION/geneNormalizedUpdated.tab "$output"/NORMALIZATION/geneNormalizedSummary.tab

	python scripts/plot_cvg_vs_completeness.py "$output"/NORMALIZATION/geneNormalizedSummary.tab "$lcompleteness" "$coverage"


	mv plotCoverage_vs_Completeness.png "$output"/PLOTS/plotCoverage_vs_Completeness.png
	
	
	
	awk 'NR==1{print $0}' "$output"/PANGENOME/gene_presence_absence.Rtab > "$output"/MATRIX/matrix.tab
	awk 'NR>1 {print $0}' "$output"/PANGENOME/gene_presence_absence.Rtab | sort -k 1 -t $'\t' >> "$output"/MATRIX/matrix.tab
	awk 'NR>1 {print $1}' "$output"/PANGENOME/gene_presence_absence.Rtab | sort -k 1 -t $'\t' > "$output"/MATRIX/INDEX
	
	awk 'NR>1 {print $1}' "$output"/NORMALIZATION/globalMeanCoverage.txt > "$output"/MATRIX/sample_names
	
	while read -r name; do
	
		echo -e "Gene\tnormalizedCoverage\tcompleteness" > "$output"/MATRIX/"${name}"_index.tmp
		grep -w "$name" "$output"/NORMALIZATION/geneNormalizedSummary.tab | awk '{print $2, $3, $NF}' >> "$output"/MATRIX/"${name}"_index.tmp
	
	done < "$output"/MATRIX/sample_names
	
	
	for i in "$output"/MATRIX/*_index.tmp; do
	
		sed -i -e 's/ /\t/g' "$i"
		python scripts/lambda.py "$i"
	
	done
	
	
	mv *_final.csv "$output"/MATRIX/
	
	for i in "$output"/MATRIX/*_final.csv; do
	
		name=$(basename "$i")
		sed  -e 's/,/\t/g' "$i" | awk 'NR>1{print $0}' > "$output"/MATRIX/"${name%_index.tmp_final.csv}"_INDEX.Z
	
	done
	
	
	for i in "$output"/MATRIX/*_INDEX.Z; do
    	name=$(basename "$i")
    	# Create the FINAL_INDEX file
    	echo "${name%_INDEX.Z}" > "$output"/MATRIX/"${name}"_FINAL_INDEX
	    	
    	# Process the INDEX file
    	while read -r gene; do
        	toprint=$(echo "$gene 0")
        	if grep -wq "$gene" "$i"; then
            	grep -w "$gene" "$i" >> "$output"/MATRIX/"${name}"_FINAL_INDEX
        	else
            	echo "$toprint" >> "$output"/MATRIX/"${name}"_FINAL_INDEX
	        fi
    	done < "$output"/MATRIX/INDEX
	    
    	# Extract the last column
    	awk '{print $NF}' "$output"/MATRIX/"${name}"_FINAL_INDEX > "$output"/MATRIX/"${name}"_FINALCOLUMN
	
	done
	
	
	paste "$output"/MATRIX/matrix.tab "$output"/MATRIX/*_FINALCOLUMN > "$output"/MATRIX/final_matrix.tab
	
	
	rm "$output"/MATRIX/*_INDEX.Z "$output"/MATRIX/*_final.csv "$output"/MATRIX/*_FINAL_INDEX "$output"/MATRIX/*INDEX.Z_FINALCOLUMN "$output"/MATRIX/*_index.tmp "$output"/MATRIX/INDEX "$output"/MATRIX/matrix.tab
	
	mv "$output"/MATRIX/final_matrix.tab "$output"/MATRIX/matrix.tab
	
	tr '\n' ' ' < "$output"/MATRIX/sample_names > "$output"/MATRIX/names_heatmap
	
	python scripts/heatmap.py "$output"/MATRIX/matrix.tab "$output"/MATRIX/names_heatmap



	mv *png "$output"/PLOTS
	""
}


workflow {
    if (params.help) {
        print_help()
    } else {
        download_genbank_fasta()
        parse_and_build_fasta_db()
        clustering_seqs()
	prokka()
	panaroo()
	alignment()
	raw_extracting()
	normalize_array()
	normalization_and_plots()
    }
}













#THIS PROCESS IS STILL IN DEVELOPTMENT, PLEASE DO NOT INCLUDE

: ' THIS NEEDS FIXING
for sample in "$data"/*; do
	samtools mpileup
	awk '{print $2,$4,$5,$6,$7}' "$output"/HETEROPLASMY/"$file" > "$output"/HETEROPLASMY/"${file}_filtered_pileup"
	awk '$3 ~ /[acgtACTG]/' "$output"/HETEROPLASMY/"${file}_filtered_pileup" > "$output"/HETEROPLASMY/"${file}_het1"
	sed 's/\([atcg]\)/\U\1/g' "$output"/HETEROPLASMY/"${file}_het1" > "$output"/HETEROPLASMY/"${file}_het2"
	sed 's/,/./g' "$output"/HETEROPLASMY/"${file}_het2" > "$output"/HETEROPLASMY/"${file}_het3"
	sed 's/\^F//g' "$output"/HETEROPLASMY/"${file}_het3" > "$output"/HETEROPLASMY/"${file}_het4"
	sed 's/\$//g' "$output"/HETEROPLASMY/"${file}_het4" > "$output"/HETEROPLASMY/"${file}_het5"
	awk '$3 !~ /-[0-9]+/ && $3 !~ /\+[0-9]+/' "$output"/HETEROPLASMY/"${file}_het5" > "$output"/HETEROPLASMY/"${file}_het6"
	#TODO WHAT THE HECK ARE THESE LINES
	awk '$3 ~ /-[0-9]+/ || $3 ~ /\+[0-9]+/' "$output"/HETEROPLASMY/"${file}_het5" > "$output"/HETEROPLASMY/"${file}_weird_lines"
	awk '$2 != 1' "$output"/HETEROPLASMY/"${file}_het6" > "$output"/HETEROPLASMY/"${file}_het7"
	awk '{ dot_count = gsub(/\./, "", $3); $5 = dot_count - $2; print }' "$output"/HETEROPLASMY/"${file}_het7" > "$output"/HETEROPLASMY/"${file}_het8"
	awk '{print $5}' "${file}_het8" | paste "$output"/HETEROPLASMY/"${file}_het7" - > "$output"/HETEROPLASMY/"${file}_het9"
	awk '{ $5 = ($5 < 0) ? -$5 : $5; print }' "$output"/HETEROPLASMY/"${file}_het9" > "$output"/HETEROPLASMY/"${file}_het10"
	awk '$5 != $2' "$output"/HETEROPLASMY/"${file}_het10" > "$output"/HETEROPLASMY/"${file}_het11"
	awk '{new_column = ($5 / $2) * 100; print $0, new_column}' "$output"/HETEROPLASMY/"${file}_het11" > "$output"/HETEROPLASMY/heterozygosis.txt
	awk '$6<=10' "$output"/HETEROPLASMY/heterozygosis.txt | wc -l > "$output"/HETEROPLASMY/less_than_or_equal_to_10%
	awk '$6>=90' "$output"/HETEROPLASMY/heterozygosis.txt | wc -l > "$output"/HETEROPLASMY/bigger_than_or_equal_to_90%
#Removing sites where there is only 1 read from the main file (OPTIONAL)
#awk 2 != 1 filtered_pileup > mpile_filtered.txt
	mv "$output"/HETEROPLASMY/*_than_* "$output"/HETEROPLASMY/distributions/
	mv "$output"/HETEROPLASMY/*_het* "$output"/HETEROPLASMY/intermediate_files/
done

'
