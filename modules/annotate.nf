/*
 * ANNOTATE{} process will annotate fasta files to generate gff files suitable for panaroo.
 */


process ANNOTATE {
	conda "${projectDir}/envs/prokka.yaml"

	input:
	path clusteredSeqsDB
	val threads
	tuple path(fasta), path(gb)
	val parallel

	output:
	path 'annotated_files/*gff', emit: prokka_gff

	script:
	"""
	#!/bin/bash

	mv $clusteredSeqsDB ./tmp_db.fasta

	seqkit translate tmp_db.fasta > clusteredSeqsDB.faa
	cp clusteredSeqsDB.faa ${params.output}/GENE_DATABASE/clustered_gene_sequences.fasta

	rm tmp_db.fasta

	ls -l *gb | awk 'NR==2{print \$NF}' > first.txt
	name=\$(cat first.txt | awk -F'/' '{print \$NF}')
	species=\$(awk '/ORGANISM/ { print \$2"_"\$3; exit }' "\${name}")

	mkdir -p annotated_files

	#annotate() function will generate gff annotation files for panaroo input.
	annotate() {
	fasta_file=\$1
	species_name=\$2

		name=\$(basename "\${fasta_file%.fasta}")
		prokka --outdir "\${name%.fna}_prokka" --species "\$species_name" --proteins clusteredSeqsDB.faa --rawproduct --cpus $threads "\${fasta_file}"
		mv "\${name%.fna}_prokka"/*gff annotated_files/"\${name}.gff"
	}
	export -f annotate
	find ./ -name "*.fasta" | parallel -j $parallel annotate {} "\$species


	cat .command.out >> ANNOTATE.log
	"""
}
