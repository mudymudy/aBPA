/*
 * ALIGNMENT_SUMMARY{} will produce simple statistics of raw coverage.
 */


process ALIGNMENT_SUMMARY {
        conda "${projectDir}/envs/alignment.yaml"
        label 'demand_3'

        input:
        path configFile
        path bamfiles
        val parallel
        val extend

        output:
		path 'postPangenomeAlignment*bam' , emit: postAlignmentFiles
        path 'completenessSummary.tab', emit: completenessSummary
        path '*_rawCoverage.txt' , emit: rawCoverage


        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/STATS

        awk '{print \$3}' $configFile | uniq > groups.txt

        while read -r groupID; do
                groupName=\$(echo "\$groupID")
                awk -v group_id="\$groupID" '\$3 == group_id { print \$1 }' $configFile > "\$groupName"ID
        done < groups.txt

        merging() {
        file=\$1
                IDs=\$(basename "\${file%ID}")

                bamFiles=()

                while read -r sampleName; do
                        bamFiles+=("\${sampleName%.fastq*}_aligned.bam")
                        echo "Adding BAM file: \${sampleName%.fastq*}_aligned.bam"
                done < "\$file"

                if [ \${#bamFiles[@]} -eq 1 ]; then
                        cp "\${bamFiles[0]}" postPangenomeAlignment_"\${IDs}".bam
                elif [ \${#bamFiles[@]} -gt 1 ]; then
                        samtools merge postPangenomeAlignment_mergedGroup"\${IDs}".bam "\${bamFiles[@]}"
                fi
        }
        export -f merging
        find ./ -name "*ID" | parallel -j $parallel merging



        #Getting some stats AND Re group the data so genotyping is simpler
        rawStats() {
        bam_file=\$1

                samplename=\$(basename "\${bam_file%.bam}")
                samtools index "\$bam_file"
                samtools depth -a "\$bam_file" > "\${samplename}_PreRawCoverage.txt"

                                #Remove the extended positions from raw coverage
                                awk -v ext=$extend '
                                {
                                        last[\$1] = \$2   # remember last position seen for each gene
                                        data[NR] = \$0   # store all lines
                                        id[NR] = \$1
                                        pos[NR] = \$2
                                }
                                END {
                                        for (g in last) {
                                        start[g] = ext + 1
                                        stop[g]  = last[g] - ext
                                        }

                                        for (i = 1; i <= NR; i++) {
                                        g = id[i]
                                        p = pos[i]
                                        if (p >= start[g] && p <= stop[g])
                                                print data[i]
                                        }
                                }'  "\${samplename}_PreRawCoverage.txt" > "\${samplename}_rawCoverage.txt"

                samtools coverage "\$bam_file" | awk -v samplename="\$samplename" 'NR>1 {print samplename, \$1, \$6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t \$'\t' >> "\${samplename}"_completenessSummary.tab
                                mv "\$bam_file" ./"\${samplename}_TMP.bam"
                                mv "\$bam_file".bai ./"\${samplename}_TMP.bam.bai"
                                picard AddOrReplaceReadGroups I="\${samplename}_TMP.bam" O="\${samplename}.bam" RGLB="\${samplename}" RGSM="\${samplename}" RGPU=Illumina RGPL=ILLUMINA RGID="\${samplename}" RGDS="\${samplename}"
                                samtools index "\${samplename}.bam"
        }
        export -f rawStats
        find ./ -name "postPangenomeAlignment*bam" | parallel -j $parallel rawStats

        rm *TMP.bam*

        cat *_completenessSummary.tab > completenessSummary.tab

        cat .command.out >> alignmentSummary.log
        """
}
