process TREE_ACCESSORY {
    conda "${projectDir}/envs/iqtree.yaml"
        
    input:
    path noUbiquitousMSA, stageAs: 'maskedMatrixGenesNoUbiquitousMSA.fasta'
    val threads

    output:
    path 'maskedMatrixGenesNoUbiquitousMSA.contree', emit: maskedMatrixGenesNoUbiquitousMSAContree , optional: true
    path 'maskedMatrixGenesNoUbiquitousMSA.iqtree', emit: maskedMatrixGenesNoUbiquitousMSAIqtree , optional: true
    path 'maskedMatrixGenesNoUbiquitousMSA.log', emit: maskedMatrixGenesNoUbiquitousMSALog , optional: true
    path 'maskedMatrixGenesNoUbiquitousMSA.treefile', emit: maskedMatrixGenesNoUbiquitousMSATreefile , optional: true

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
	}' maskedMatrixGenesNoUbiquitousMSA.fasta  > to_remove

    if [[ -s to_remove ]]; then
        awk 'NR==FNR {remove[\$0]; next}
            /^>/ {keep = !(\$0 in remove)}
            keep' to_remove maskedMatrixGenesNoUbiquitousMSA.fasta > accessory_genome.fasta
    else
        mv maskedMatrixGenesNoUbiquitousMSA.fasta accessory_genome.fasta
    fi


    iqtree -s accessory_genome.fasta --prefix accessory_genome -T $threads -B 1000 -m MFP

	cp accessory_genome* ${params.output}/TREE/
    """
}
