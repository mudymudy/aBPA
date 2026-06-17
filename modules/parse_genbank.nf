/*
 * PARSE_GENBANK() will perform genbank format quality check.
 */

 
process PARSE_GENBANK {
	conda "${projectDir}/envs/biopython.yaml"
	
	input:
	path gffFiles
	path fastaFiles
		
	output:
	tuple path('cleaned/*fasta'), path('cleaned/*gb'), emit: validFiles

	script:
	"""
	#!/bin/bash

	mkdir -p cleaned

	parseTest.py > parseTest.txt
	grep "is not a valid GenBank file" parseTest.txt | awk '{print \$1}' > blackListed.txt

	if [[ -s blackListed.txt ]] ;then

		while read -r removeMe; do
			rm "\${removeMe%gb}fasta" *"\$removeMe"	
  			echo -e "\${removeMe%.gb} sample has been removed due to format problems."
		done < blackListed.txt
	else
		echo -e "Every file passed the test. Moving on."
		mv *.gb cleaned/
		mv *.fasta cleaned/
	fi
	"""
}
