process ALIGNMENT {
        conda "${projectDir}/envs/alignment.yaml"

        input:
        path data
        path panRef
        path configFile
        tuple val(threadsGlobal), val(missingProb), val(seedAlignment), val(gapFraction),val(minReadLength),val(maxReadLength),val(parallel), val(quality)
        val rescale

        output:
        path '*_aligned.bam', emit: postAlignedBams
        path '*_final.fastq', emit: postAlignedReads, optional: true
        path 'extended_pangenome_reference.fasta.*' , emit: pan_index
        path 'results*', emit: mapdamage_results, optional: true

        script:
        """
        #!/bin/bash

        mkdir -p ${params.output}/ALIGNMENT
        mkdir -p ${params.output}/MAPDAMAGE

        bwa index $panRef

        #align() will align the data and perform several post-alignment computations
        align() {
        sample=\$1

                name=\$(basename "\${sample%.fastq*}")
                filename=\$(basename "\${sample}")

                # if --test
                if [ -f "./test.fastq.gz" ]; then
                        softClip=5
                else
                       softClip=\$(awk -v file_name="\$filename" '\$1 == file_name { print \$2 }' $configFile)
                fi

                status=\$(awk -v f="\$filename" '\$1==f {print \$4}' $configFile)
                if [[ "\$status" == "A" ]]; then
                echo "Using bwa aln for \$filename" >> LOGFILE

                            # Making read groups
                            rg_id="\${name}"  # sample name as id
                            rg_sm="\${name}" # sample name again
                            rg_pl="illumina"        # I dont think this is very important for this pipeline so its going to be just illumina because why not
                            rg_lb="lib1"            # group id
                            rg_pu="unit1"           # not sure what Ill put here

                            bwa aln -l $seedAlignment -n $missingProb -o $gapFraction -t $threadsGlobal $panRef "\$sample" > "\${name}.sai"
                            bwa samse -r "@RG\\tID:\$rg_id\\tSM:\$rg_sm\\tPL:\$rg_pl\\tLB:\$rg_lb\\tPU:\$rg_pu" \
                            $panRef "\${name}.sai" "\$sample" > "\${name}.sam"

                            samtools view -bS "\${name}.sam" > "\${name}.bam"
                            samtools quickcheck "\${name}.bam"
                            samtools sort -o "\${name}_sorted.bam" -O bam -@ $threadsGlobal "\${name}.bam"
                            samtools index "\${name}_sorted.bam"
                            rm "\${name}.bam"
                            samtools view -b -@ $threadsGlobal -F 4 "\${name}_sorted.bam" > "\${name}_sorted_mappedreads.bam"
                            samtools index "\${name}_sorted_mappedreads.bam"

                        if [[ $rescale -eq 1 ]]; then
                            mapDamage --rescale --merge-reference-sequences -i "\${name}_sorted_mappedreads.bam" -r $panRef
                            mv results_"\${name}_sorted_mappedreads"/*.rescaled.bam "\${name}_softclipped.bam"
                            picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=TRUE I="\${name}_softclipped.bam" O="\${name}_deduped.bam" M="\${name}_deduped.stats"
                            samtools view -q $quality -o "\${name}_qc.bam" "\${name}_deduped.bam"
                            samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name}_lg.bam" "\${name}_qc.bam"
                            samtools sort -o "\${name}_aligned.bam" -O bam -@ $threadsGlobal "\${name}_lg.bam"
                            samtools coverage "\${name}_aligned.bam" > "\${name}_genomicsMetrics.txt"
                            samtools fastq -@ $threadsGlobal "\${name}_aligned.bam" > "\${name}_final.fastq"
                        else
                            bam trimBam "\${name}_sorted_mappedreads.bam" "\${name}_softclipped.bam" -L "\$softClip" -R "\$softClip" --clip
                            mapDamage --merge-reference-sequences -i "\${name}_sorted_mappedreads.bam" -r $panRef
                            picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=TRUE I="\${name}_softclipped.bam" O="\${name}_deduped.bam" M="\${name}_deduped.stats"
                            samtools view -q $quality -o "\${name}_qc.bam" "\${name}_deduped.bam"
                            samtools view -e 'length(seq)>$minReadLength && length(seq)<$maxReadLength' -O BAM -o "\${name}_lg.bam" "\${name}_qc.bam"
                            samtools sort -o "\${name}_aligned.bam" -O bam -@ $threadsGlobal "\${name}_lg.bam"
                            samtools coverage "\${name}_aligned.bam" > "\${name}_genomicsMetrics.txt"
                            samtools fastq -@ $threadsGlobal "\${name}_aligned.bam" > "\${name}_final.fastq"
                        fi

                    else
                        echo "Using bwa mem for \$filename" >> LOGFILE
                        bwa mem -B 1 -E 1 $panRef "\$sample" -t $threadsGlobal > "\$name".sam
                        samtools view -bS "\$name".sam > "\$name".bam
                        samtools quickcheck "\$name".bam
                        samtools sort -o "\$name"Sorted.bam -O bam -@ $threadsGlobal "\$name".bam
                        picard MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 REMOVE_DUPLICATES=TRUE I="\$name"Sorted.bam O="\${name}_deduped.bam" M="\${name}_deduped.stats"
                        samtools index "\${name}_deduped.bam"
                        samtools view -b -@ 10 -F 4 "\${name}_deduped.bam" > "\$name"SortedMappedreads.bam
                        samtools index "\$name"SortedMappedreads.bam
                        samtools sort -o "\$name"_aligned.bam -O bam -@ $threadsGlobal "\$name"SortedMappedreads.bam
                    fi
        }
        export -f align


        #if --test
        if [ -f "./test.fastq.gz" ]; then
            find ./ -name "test.fastq.gz" | parallel -j $parallel align
        else
            find $data/* -name "*.fastq*" | parallel -j $parallel align
        fi


        shopt -s nullglob
        #remove files that can be optional
        sam_files=( *.sam )
        (( \${#sam_files[@]} )) && rm "\${sam_files[@]}"

        sai_files=( *.sai )
        (( \${#sai_files[@]} )) && rm "\${sai_files[@]}"

        lg_files=( *_lg.bam )
        (( \${#lg_files[@]} )) && rm "\${lg_files[@]}"

        qc_files=( *_qc.bam )
        (( \${#qc_files[@]} )) && rm "\${qc_files[@]}"

        ancient_fastq_files=( *_final.fastq )
        mapdamage_results=( results_* )

        cp *_aligned.bam ${params.output}/ALIGNMENT

        if (( \${#ancient_fastq_files[@]} )); then
                cp *_final.fastq ${params.output}/ALIGNMENT
        fi

        if (( \${#mapdamage_results[@]} )); then
                cp -r results_* ${params.output}/MAPDAMAGE
        fi

        shopt -u nullglob

        cat .command.out >> alignment.log
        """
}
