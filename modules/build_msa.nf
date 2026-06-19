/*
 * BUILD_MSA{} will concatenate individual gene MSA.
 */

process BUILD_MSA {
	conda "${projectDir}/envs/normalization.yaml"


	input:
	path genesAlnSeq, stageAs: 'genes/*'
	path maskedMatrixGenesNoUbiquitous, stageAs: 'filesList_maskedMatrixGenesNoUbiquitous.txt'
	path maskedMatrixGenesOnlyAncient, stageAs: 'filesList_maskedMatrixGenesOnlyAncient.txt'
	path maskedMatrixGenesUbiquitous, stageAs: 'filesList_maskedMatrixGenesUbiquitous.txt'
	path genesAbovePercentSeries, stageAs: 'filesList_genesAbovePercentSeries.txt'
	path sampleNames, stageAs: 'sampleNames.txt'
	val parallel

	output:
	path 'genesAbovePercentMSA.fasta', emit: genesAbovePercentMSA
	path 'maskedMatrixGenesNoUbiquitousMSA.fasta', emit: maskedMatrixGenesNoUbiquitousMSA
	path 'maskedMatrixGenesOnlyAncientMSA.fasta', emit: maskedMatrixGenesOnlyAncientMSA
	path 'maskedMatrixGenesUbiquitousMSA.fasta', emit: maskedMatrixGenesUbiquitousMSA
	
	script:
	"""
	#!/bin/bash

	mkdir -p ${params.output}/GENE_MSA
	mkdir -p ${params.output}/TREE

	touch maskedMatrixGenesNoUbiquitous.fasta
	touch maskedMatrixGenesOnlyAncient.fasta
	touch maskedMatrixGenesUbiquitous.fasta
	touch genesAbovePercentSeries.fasta

	sed -i -e 's/~/_/g' $genesAbovePercentSeries
	sed -i -e 's/~/_/g' $maskedMatrixGenesNoUbiquitous
	sed -i -e 's/~/_/g' $maskedMatrixGenesOnlyAncient
	sed -i -e 's/~/_/g' $maskedMatrixGenesUbiquitous

    #quality check step to remove sites that contain only non-nucleotide characters: n OR N or -.
    scanning_and_trimming() {
    fasta_file=\$1

    awk -v fasta_filename="\$fasta_file" '
        /^>/ {
            headers[++count] = \$0
            next
        }
        {
            seqs[count] = \$0
            if (count == 1) #store the alignment length
                aln_len = length(\$0)
        }

        END {
            #for each column, check if ALL sites have n OR N OR -
            for (col = 1; col <= aln_len; col++) {
                bad = 1
                for (i = 1; i <= count; i++) {
                    c = substr(seqs[i], col, 1)
                    if (c != "n" && c != "N" && c != "-") {
                        bad = 0
                        break
                    }
                }
                keep[col] = !bad   #1 = survives, 0 = trim
            }

            #count trimmed columns
            dropped = 0
            for (col = 1; col <= aln_len; col++)
                if (!keep[col]) dropped++

            if (dropped > 0)
                print dropped " columns were trimmed from" fasta_filename " alignment as they had only non nucleotide characters" >  "/dev/stderr"
            else
                print "No artifact columns containing only non nucleotide characters were found in" fasta_filename " after scanning" >  "/dev/stderr"

            #print each sequence using only kept columns
            for (i = 1; i <= count; i++) {
                print headers[i]
            out     = ""
                for (col = 1; col <= aln_len; col++)
                    if (keep[col])
                        out = out substr(seqs[i], col, 1)
                print out
            }
    }' "\$fasta_file" > "\${fasta_file}_tmp" && mv "\${fasta_file}_tmp" "\$fasta_file"
    }
    export -f scanning_and_trimming
    find genes/ -name "*fasta" | parallel -j $parallel scanning_and_trimming

	build_msa() {
	inputfile=\$1
	filename=\$(basename "\${inputfile}")
	filename="\${filename#filesList_}"
	filename="\${filename%.txt}"

		# Build a list of existing files
		files=()
		while read -r gene; do
	    	file="genes/\${gene}.fasta"
    		if [[ -f "\$file" ]]; then
	       		files+=("\$file")
    		else
	       		echo "Warning: File \$file not found, skipping..."
    		fi
		done < "\$inputfile"
	
		#append genes array to MSA into array
		files=("\$filename".fasta "\${files[@]}")
	
		#paste everything inside array
		paste "\${files[@]}" | tr -d '\t' > TMP_"\${filename}"
		mv TMP_"\${filename}" "\${filename}".fasta
	
		#Clean headers (if needed)
		awk -F'>' '/^>/ {print ">" \$2} !/^>/' "\${filename}".fasta > TMPg_"\${filename}"
		mv TMPg_"\${filename}" "\${filename}".fasta
	}
	export -f build_msa
	find ./ -name "filesList_*" | parallel -j $parallel build_msa

	#rename for now
	mv maskedMatrixGenesNoUbiquitous.fasta maskedMatrixGenesNoUbiquitousMSA.fasta
	mv maskedMatrixGenesOnlyAncient.fasta maskedMatrixGenesOnlyAncientMSA.fasta
	mv maskedMatrixGenesUbiquitous.fasta maskedMatrixGenesUbiquitousMSA.fasta
	mv genesAbovePercentSeries.fasta genesAbovePercentMSA.fasta

	cp -r ./genes/*fasta ${params.output}/GENE_MSA
	cp genesAbovePercentMSA.fasta ${params.output}/TREE/threshold.fasta
	cp maskedMatrixGenesUbiquitousMSA.fasta ${params.output}/TREE/core.fasta
	cp maskedMatrixGenesOnlyAncientMSA.fasta ${params.output}/TREE/ancient.fasta
	cp maskedMatrixGenesNoUbiquitousMSA.fasta ${params.output}/TREE/accessory.fasta
	"""
}
