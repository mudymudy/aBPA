
process TREE_ANCIENT {
    conda "${projectDir}/envs/iqtree.yaml"
       
    input:
    path ancientMSA, stageAs: 'maskedMatrixGenesOnlyAncientMSA.fasta'
    val threads

    output:
    path 'maskedMatrixGenesOnlyAncientMSA.contree', emit: maskedMatrixGenesOnlyAncientMSAConqtree , optional: true
	path 'maskedMatrixGenesOnlyAncientMSA.iqtree', emit: maskedMatrixGenesOnlyAncientMSAIqtree , optional: true
	path 'maskedMatrixGenesOnlyAncientMSA.log', emit: maskedMatrixGenesOnlyAncientMSALog , optional: true
	path 'maskedMatrixGenesOnlyAncientMSA.treefile', emit: maskedMatrixGenesOnlyAncientMSATreefile , optional: true

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
	}' maskedMatrixGenesOnlyAncientMSA.fasta  > to_remove

    if [[ -s to_remove ]]; then
        awk 'NR==FNR {remove[\$0]; next}
        /^>/ {keep = !(\$0 in remove)}
        keep' to_remove maskedMatrixGenesOnlyAncientMSA.fasta > ancient_set.fasta
    else
        mv maskedMatrixGenesOnlyAncientMSA.fasta ancient_set.fasta
    fi

    iqtree -s ancient_set.fasta --prefix ancient_set -T $threads -B 1000 -m MFP

	cp ancient_set* ${params.output}/TREE/
    """
}
