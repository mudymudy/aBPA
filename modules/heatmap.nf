/*
 * HEATMAP{} will generate an index file with genes presence/absence list after cutoffs, a set of lists with genes based on conditions and plot heatmaps.
 */

process HEATMAP {
        conda "${projectDir}/envs/heatmap.yaml"
        label 'demand_4'

        input:
        path fCSV, stageAs: 'fCSV/*'
        path INDEX, stageAs: 'INDEX/*'
        path matrix, stageAs: 'matrix/*'
        path names, stageAs: 'names/*'
        val threshold_value
        val samples_heatmap
        val plot_heatmap

        output:
        path 'final_matrix.tab', emit: finalMatrix
        path 'presenceAbsence*.png', emit: presenceAbsence, optional: true
        path 'maskedMatrixGenesOnlyAncient.txt', emit: maskedMatrixGenesOnlyAncient
        path 'maskedMatrixGenesUbiquitous.txt', emit: maskedMatrixGenesUbiquitous
        path 'maskedMatrixGenesNoUbiquitous.txt', emit: maskedMatrixGenesNoUbiquitous
        path 'sampleOrdernoUbiquitous*.txt', emit: sampleOrdernoUbiquitous, optional: true
        path 'sampleOrderonlyAncient*.txt', emit: sampleOrderonlyAncient, optional: true
        path 'genesAbovePercentSeries.txt', emit: genesAbovePercentSeries
        path 'blackListedQualityChecked.txt', emit: blackListed
        path '*_presence_absence_genes.index' , emit: genesIndex

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/MATRIX
        mkdir -p ${params.output}/PLOTS


        formatting_fcsv() {
        i=\$1
                name=\$(basename "\${i%_index.tmp_final.csv}")
                sed  -e 's/,/\t/g' "\$i" | awk 'NR>1{print \$0}' > "\${name}"_present_genes.txt
        }
        export -f formatting_fcsv
        find ./fCSV/ -name "*_final.csv" | parallel -j $task.cpus formatting_fcsv


        add_missing_genes() {
        i=\$1
                name=\$(basename "\${i%_present_genes.txt}")
                #Create the presence_absence_genes.index file
                echo "\${name}" > "\${name}"_presence_absence_genes.index

                #Process the INDEX file
                while read -r gene; do
                        toprint=\$(echo "\$gene 0")
                        if grep -wq "\$gene" "\$i"; then
                                grep -w "\$gene" "\$i" >> "\${name}"_presence_absence_genes.index
                        else
                                echo "\$toprint" >> "\${name}"_presence_absence_genes.index
                        fi
                done < INDEX/INDEX

                #Extract the last column
                awk '{print \$NF}' "\${name}"_presence_absence_genes.index > "\${name}"_last_column.txt
        }
        export -f add_missing_genes
        find ./ -name "*present_genes.txt" | parallel -j $task.cpus add_missing_genes


        # I need to add a quality control step right here. User samples can be false positives sometimes or just super low quality and have 0 genes after filtering
    # Then, black list unwanted samples and exclude them from the final_matrix.tab document.

        qual_control() {
        checkSample=\$1

        sampleName=\$(basename "\${checkSample%_last_column.txt}")
                counts=\$(grep -c "0" "\$checkSample")
        totalLines=\$(wc -l "\$checkSample" | awk '{print \$1 - 1}')
        proportion=\$(awk -v absence="\$counts" -v record="\$totalLines" 'BEGIN { print (absence / record ) }')

                if  (( \$(awk -v p="\$proportion" 'BEGIN { print (p > 0.95) }' ) )); then
                echo "\$sampleName" >> blackListedQualityChecked.txt
            fi
                }
        export -f qual_control
        find ./ -name "*_last_column.txt" | parallel -j $task.cpus qual_control

        # Do this only if blacklisted
        if [[ -s blackListedQualityChecked.txt ]]; then
                while read -r removeMe; do
                        mv "\${removeMe}_last_column.txt" "\${removeMe}LowQualitySample"
                        grep -v "\${removeMe}" names/sample_names > names/sample_names.tmp
                        mv names/sample_names.tmp names/sample_names
                done < blackListedQualityChecked.txt
        else
                touch blackListedQualityChecked.txt
        fi

        #making an empty matrix to avoid sample duplication
        awk 'BEGIN {FS=OFS="\t"} {print \$1}' matrix/matrix.tab > empty_matrix.tab
        paste empty_matrix.tab *_last_column.txt > final_matrix.tab
        tr '\n' ' ' < names/sample_names > names_heatmap

        heatmap.py final_matrix.tab names_heatmap $threshold_value $samples_heatmap $plot_heatmap

        cp final_matrix.tab ${params.output}/MATRIX/gene_presence_absence_matrix.tab

        if [[ $plot_heatmap == 1 ]]; then
                cp *png ${params.output}/PLOTS
        fi
        """
}
