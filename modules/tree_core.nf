

process TREE_CORE {
    conda "${projectDir}/envs/iqtree.yaml"
        
    input:
    path ubiquitousMSA, stageAs: 'maskedMatrixGenesUbiquitousMSA.fasta'
    val threads

    output:
    path 'core_genome.contree', optional: true
    path 'core_genome.iqtree', optional: true
    path 'core_genome.log', optional: true
	path 'core_genome.treefile', optional: true

    script:
    """
    #!/bin/bash

	mkdir -p ${params.output}/TREE/


	awk '/^>/ {
    	header = \$0
    	getline seq
    	if (seq ~ /^[nN-]+\$/) {
		    print header
    	}
	}' maskedMatrixGenesUbiquitousMSA.fasta  > to_remove

    if [[ -s to_remove ]]; then
        awk 'NR==FNR {remove[\$0]; next}
            /^>/ {keep = !(\$0 in remove)}
            keep' to_remove maskedMatrixGenesUbiquitousMSA.fasta > core_genome.fasta
    else
        mv maskedMatrixGenesUbiquitousMSA.fasta core_genome.fasta
    fi

    iqtree -s core_genome.fasta --prefix core_genome -T $threads -B 1000 -m MFP

	cp core_genome* ${params.output}/TREE/
    """
}
