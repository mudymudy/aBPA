/*
 * MAKE_PANGENOME{} will build a gene pangenome based on prokka gff output.
 */
 
process MAKE_PANGENOME {
	conda "${projectDir}/envs/panaroo.yaml"
	
	input:
	path prokka_gff
	val pangenomeMode
	val pangenomeThreshold
	val threads
	val alignment

	output:
	path 'pan_genome_reference.fa' , emit: panSequence
	path 'gene_presence_absence.Rtab' , emit: initialMatrix
	path 'aligned_gene_sequences/*' , emit: alignedGenesSeqs
	tuple path('final_graph.gml'), path('gene_data.csv'), emit: pangenome_metadata
	path 'gene_list.txt', emit: gene_list

	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/PANGENOME

	#make pangenome
	panaroo -i *.gff -o ./ --clean-mode $pangenomeMode -a $alignment --merge_paralogs --core_threshold $pangenomeThreshold -t $threads

	#output gene list from aligned sequences
    gene_list() {
    file=\$1
                name=\$(basename "\$file")
                echo -e "\${name}" | sed -e 's/.aln.fas//g' -e 's/.fasta//g' -e 's/aligned_gene_sequences//g' -e '/^\$/d' -e 's/~/_/g' >> gene_list.txt
    }
    export -f gene_list
    find ./aligned_gene_sequences/ -name "*" | parallel -j $threads gene_list

	#replace ~ with _ on final graph
	sed -i -e 's/~/_/g' final_graph.gml

	cat .command.out >> makePangenome.log
	cp makePangenome.log ${params.output}/PANGENOME
	cp gene_data.csv ${params.output}/PANGENOME
	cp final_graph.gml ${params.output}/PANGENOME
	cp summary_statistics.txt ${params.output}/PANGENOME
	"""
}
