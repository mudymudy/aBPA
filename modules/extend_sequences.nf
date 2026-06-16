/*
 * EXTEND_sEQUENCES{} process will extend gene sequences from pangenome reference.
 */


process EXTEND_SEQUENCES {

	conda "${projectDir}/envs/extend_sequences.yaml"

	input:
	path gff
	path fasta
	tuple path(final_graph), path(gene_data)
	val extend
	val parallel
	path gene_list

	output:
	path 'extended_pangenome_reference_sequence.fasta' , emit: extended_reference
	path 'pangenome_length.txt', emit: pangenome_length
	path 'unextended_pangenome_reference.fasta', emit: unextended_reference
	path 'final_list_genes.txt', emit: final_list_genes

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	#fetch data from graph
	awk '
		BEGIN {OFS= "\t"}

		/^[[:space:]]*node/ {      #if node, activate the script and re-set values
	
			activator = 1
			centroid = ""
			gene_name = ""
			next
	
		}
	
		activator {
	
			if (\$1 == "longCentroidID" && centroid == "" && \$2 ~ /^"[0-9]/ && \$2 !~ /refound/) {  #id has to start with number
				split(\$2, parts, ";")
				centroid = parts[1]

			}
	
				else if (\$1 == "name") {
				gene_name = \$2
				print centroid, gene_name
				activator = 0
	
			}
		}' $final_graph > seq_ids_per_node
	
	sed -i -e 's/"//g' -e 's/ /\t/g' seq_ids_per_node
	
	
	while read -r seqID gene; do
		awk -F',' -v cluster_id="\$seqID" -v gene_name="\$gene" '\$3 == cluster_id { print \$1, \$3, \$4, gene_name }' $gene_data >> updated_seq_ids_per_node
	done < seq_ids_per_node


	#Now we need to do some cleaning
	awk '{
		centroids[\$2] = centroids[\$2]"\t"\$4   #\$2 and \$3 should be the same for different lines if there are more copies of a centroid
		centroid_count[\$2]++
	}
	END { 
		for (id in centroids) {
			if ( centroid_count[id] > 1) {
				print id, centroids[id], centroid_count[id]
			}
		}
	}' updated_seq_ids_per_node > redundant_genes.txt

	awk '{print \$1}' redundant_genes.txt > exclude_centroids.txt

	awk 'FNR==NR { ex[\$1]=1; next }
	{
  	if (\$2 in ex) {
	    if (!seen[\$2]++) {
      	print > "updated_seq_ids_per_node.cleaned"   #final file
    	} else {
      	print > "redundant_centroids_lines.txt"      #exclude redundant lines
    	}
  	} else {
	    print > "updated_seq_ids_per_node.cleaned"
  	}
	}
	' exclude_centroids.txt updated_seq_ids_per_node
	

	mv updated_seq_ids_per_node.cleaned ./updated_seq_ids_per_node


	
	#now we read updated_seq_ids_per_node to generate new gff files

	while read -r sample seqID gene_tag gene_name; do

		grep -w "\$gene_tag" "\${sample}".gff >> "\${sample}"_selected_lines.gff

	done < updated_seq_ids_per_node


	#now get the bed files from _selected_lines.gff

	make_bed() {
	file=\$1
	name=\$(basename "\${file%_selected_lines.gff}")

		awk '\$3 == "CDS" || \$3 == "gene" {
    		sample = \$1
    		start = \$4 - 1  #GFF is 1-based, BED is 0-based
    		end = \$5
    		strand = \$7
    		attr = \$9

    		name = "."
    		if (attr ~ /ID=/) {
        		match(attr, /ID=[^;]+/)
        		name = substr(attr, RSTART+3, RLENGTH-3)
    		}

    		print sample, start, end, name, ".", strand
		}' OFS="\t" "\$file" > "\$name".bed

	}
	export -f make_bed
	find ./ -name "*_selected_lines.gff" | parallel -j $parallel make_bed


	#finally faidx fasta files and then do bedtools slop

	index_fasta() {
	file=\$1

		samtools faidx "\$file"
	}
	export -f index_fasta
	find ./ -name "*.fasta" | parallel -j $parallel index_fasta

	#now extend the sequences
	extend_sequences() {
	file=\$1
	name=\$(basename "\${file%.bed}")

		#Make extended individual sequences
		bedtools slop -i "\$file" -g "\$name".fasta.fai -b $extend > "\$name"_extended_sequences.bed
		bedtools getfasta -bed "\$name"_extended_sequences.bed -fi "\$name".fasta -name > "\$name"_extended_sequences.fasta

		#here I need to check: if sequence is in - strand, then get the complementary reverse.
		awk '\$NF == "-" {print \$4}' "\$name"_extended_sequences.bed > "\$name"_negative_strand
		awk '\$NF == "+" {print \$4}' "\$name"_extended_sequences.bed > "\$name"_positive_strand

		#get positive lines
		while read -r gene_tag; do
			grep -A1 -w "\$gene_tag" "\$name"_extended_sequences.fasta >> "\$name"_tmp_positive.fasta
		done < "\$name"_positive_strand

		#get negative lines and reverse them
		while read -r gene_tag; do
    			#get the header+sequence
    			grep -A1 -w "\$gene_tag" "\$name"_extended_sequences.fasta > "\$name"_tmp_seq.fasta

 				if [ -s "\${name}_tmp_seq.fasta" ]; then
        			seqtk seq -r "\${name}_tmp_seq.fasta" >> "\${name}_reverse_complement.fasta"
    			fi

		done < "\$name"_negative_strand

		cat "\$name"_tmp_positive.fasta "\$name"_reverse_complement.fasta 2>/dev/null > "\$name"_extended_sequences.fasta

	}
	export -f extend_sequences
	find ./ -name "*.bed" | parallel -j $parallel extend_sequences

	rm *_extended_sequences.bed


	#make unextended reference
	unextended() {
	file=\$1
	name=\$(basename "\${file%.bed}")


		bedtools getfasta -bed "\$file" -fi "\$name".fasta -name > "\$name"_unextended_sequences.fasta

		#if sequence is in - strand, then get the complementary reverse but for un-extended
		awk '\$NF == "-" {print \$4}' "\$file" > "\$name"_negative_strand_unextended
		awk '\$NF == "+" {print \$4}' "\$file" > "\$name"_positive_strand_unextended


		while read -r gene_tag; do
			grep -A1 -w "\$gene_tag" "\$name"_unextended_sequences.fasta >> "\$name"_tmp_positive_unextended.fasta
		done < "\$name"_positive_strand_unextended

		while read -r gene_tag; do
    			#get the header+sequence
    			grep -A1 -w "\$gene_tag" "\$name"_unextended_sequences.fasta > "\$name"_tmp_seq_unextended.fasta

    			#apply reverse complement
 				if [ -s "\${name}_tmp_seq_unextended.fasta" ]; then
	    			seqtk seq -r "\$name"_tmp_seq_unextended.fasta >> "\$name"_reverse_complement_unextended.fasta
				fi

		done < "\$name"_negative_strand_unextended

		cat "\$name"_tmp_positive_unextended.fasta "\$name"_reverse_complement_unextended.fasta 2>/dev/null > "\$name"_unextended_sequences.fasta

	}
	export -f unextended
	find ./ -name "*.bed" | parallel -j $parallel unextended



	#now we need the gene names back!


	rename_extended_sequences() {
	file=\$1
	name=\$(basename "\${file%_extended_sequences.fasta}")

		awk '

			FNR==NR {

				gene_tag = \$3
				gene_name = \$4
				pair[gene_tag] = gene_name
				next
			}
			
			/^>/ { 
		
				split(substr(\$0,2), subStrings, "::")
				gene_tag_header = subStrings[1]
	
				if (gene_tag_header in pair) {
	
					print ">" pair[gene_tag_header]
	
				} 
	
				next
	
			}
	
			{
	
				print
	
			}
			'	updated_seq_ids_per_node "\$file" > "\$name"_extended_reference.fasta
	
	}
	export -f rename_extended_sequences
	find ./ -name "*_extended_sequences.fasta" | parallel -j $parallel rename_extended_sequences
	



	rename_unextended_sequences() {
	file=\$1
	name=\$(basename "\${file%_unextended_sequences.fasta}")

		awk '

			FNR==NR {

				gene_tag = \$3
				gene_name = \$4
				pair[gene_tag] = gene_name
				next
			}
			
			/^>/ { 
		
				split(substr(\$0,2), subStrings, "::")
				gene_tag_header = subStrings[1]
	
				if (gene_tag_header in pair) {
	
					print ">" pair[gene_tag_header]
	
				} 
	
				next
	
			}
	
			{
	
				print
	
			}
			'	updated_seq_ids_per_node "\$file" > "\$name"_unextended_reference.fasta
	
	}
	export -f rename_unextended_sequences
	find ./ -name "*_unextended_sequences.fasta" | parallel -j $parallel rename_unextended_sequences



	#now we just concatenate them

	cat *_extended_reference.fasta > pre_extended_pangenome_reference_sequence.fasta
	cat *_unextended_reference.fasta > pre_unextended_pangenome_reference_sequence.fasta


	#QC for gene list
	while read -r gene; do
		grep -A 1 -w "\$gene" pre_extended_pangenome_reference_sequence.fasta >> extended_pangenome_reference_sequence.fasta
		grep -A 1 -w "\$gene" pre_unextended_pangenome_reference_sequence.fasta >> unextended_pangenome_reference_sequence.fasta
	done < $gene_list


	#get pangenome length
	seqtk seq unextended_pangenome_reference_sequence.fasta | awk '!/^>/ {line_length += length(\$0)} END {print line_length}' > pangenome_length.txt

	#Report
	grep -w "name" final_graph.gml | awk '{print \$2}' | sed -e 's/"//g' > genes_in_graph.txt

	graphGenes=\$(wc -l < genes_in_graph.txt)
	echo -e "Genes in final_graph.gml: \$graphGenes" >> Pangenome_report.txt

	msaGenes=\$(wc -l < gene_list.txt)
	echo -e "Genes that have MSA from Panaroo: \$msaGenes" >> Pangenome_report.txt

	while read -r geneID; do
	    if grep -qw "\$geneID" gene_list.txt; then
        	echo -e "Gene \$geneID from graph was found to have a multiple sequence alignment from Panaroo" >> Pangenome_report.txt
    	else
	        echo -e "WARNING: Gene \$geneID from graph was not found to have a multiple sequence alignment from Panaroo." >> Pangenome_report.txt
	    fi
	done < genes_in_graph.txt

	#final list of genes
	awk '/^>/ {print \$0}' unextended_pangenome_reference_sequence.fasta | sed -e 's/>//g' | sort > final_list_genes.txt
	
	#output
	mv unextended_pangenome_reference_sequence.fasta ./unextended_pangenome_reference.fasta
	cp unextended_pangenome_reference.fasta ${params.output}/PANGENOME
	cp Pangenome_report.txt ${params.output}/PANGENOME
	"""
}
