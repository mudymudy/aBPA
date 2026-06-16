/*
 * FILTER_GENE_ALIGNMENTS{} will complete each gene MSA and parse them to find artifacts.
 */


process FILTER_GENE_ALIGNMENTS {
        conda "${projectDir}/envs/filterGeneAlignments.yaml"

        input:
        path genesAln, stageAs: 'panaroo_genes/*'
        path extractedSequencesFasta, stageAs: 'user_genes/*'
        path fFiles, stageAs: 'modern_data/*'
        val genomes
        path outgroupSeq, stageAs: 'outgroup'
        path blackListed, stageAs: 'blackListed.txt'
        val parallel
        path final_list_genes
        path genes_2_mask
        path config
        path global_coverage, stageAs: 'synthetic/*'

        output:
        path 'sorted/*.fasta', emit: genesAlnSeq
        path 'sampleNames.txt', emit: sampleNames

        script:
        """
        #!/bin/bash

        shopt -s nullglob  #Prevent literal interpretation of wildcards if there are no matchings
        #sanity check before starting the process

        if ! find panaroo_genes/ -name "*.aln.fas" | grep -q .; then
                echo -e "No files with .aln.fas extension were found in panaroo_genes/. Check previous process. Stopping the pipeline."
                exit 1
        else
                echo -e "Fasta files with expected extension .aln.fas were found. Proceeding with the process."
        fi

        if [[ -e outgroup ]]; then
                sed -i -e 's/~/_/g' outgroup
        else
                echo -e "outgroup file was not found. Stopping the pipeline."
                exit 1
        fi

        if ! find modern_data/ -name "*fasta" | grep -q .; then
                echo -e "No modern sequences with .fasta extension were found in modern_data/. Check previous process. Stopping the pipeline."
                exit 1
        else
                echo -e "Fasta files with expected extension .fasta were found. Proceeding with the process."
        fi

        #Remove duplicated genes from the dataset
        mkdir ./redundant_genes

        #valid genes for matching
        grep -F -f $final_list_genes /dev/null >/dev/null 2>&1

        #move all files NOT in the list
        for file in panaroo_genes/*; do
            base=\$(basename "\$file")
            gene=\${base%%.*}

            # If gene not in list, move it
        if ! grep -Fxq "\$gene" $final_list_genes; then
                        mv "\$file" ./redundant_genes/
        fi
        done

        renaming() {
        file=\$1
        new_name=\$(basename "\${file}")
        new_name="\${new_name#extractedSequences}"
        new_name="\${new_name#mergedGroup}"

                mv "\$file" user_genes/"\$new_name"
        }
        export -f renaming
        find ./user_genes -name "*fasta" | parallel -j $parallel renaming


        #if user input data contains modern
        awk '\$4 == "M" { print \$3 }' $config > modern_group_as_input.txt

        input_as_modern_activator=0
        if [[ -s modern_group_as_input.txt ]]; then
                mkdir -p input_modern
                input_as_modern_activator=1

                while read -r sample;do
			ls -l ./user_genes/"\${sample}.fasta"

                        mv ./user_genes/"\${sample}.fasta" input_modern/
                done < modern_group_as_input.txt
        fi


        mkdir -p blacklisted
        echo -e "Checking if there are low quality samples"

        if [[ -s blackListed.txt ]]; then

                while read -r removeMe; do

			echo "Checking file existence just before mv:"
			ls -l ./input_modern/"\${removeMe}.fasta"

                		status=\$(awk -v sample_to_remove="\$removeMe" '\$3 == sample_to_remove { print \$4 }' config_test.tab)
                		if [[ "\$status" == "A" ]]; then
				
					echo "Sample to move: \$removeMe"
                        		mv ./user_genes/*"\${removeMe}.fasta" ./blacklisted/
	                        	echo "\${removeMe} has been removed from analysis due to low quality."

         	               elif [[ "\$status" == "M" ]]; then
				
					echo "Sample to move: \$removeMe"
					ls -l ./input_modern/"\${removeMe}.fasta"
	                        	mv ./input_modern/"\${removeMe}.fasta" ./blacklisted/
	                        	echo "\${removeMe} has been removed from analysis due to low quality.
                	        fi
                done < blackListed.txt
        else
                echo -e "Every sample passed quality checks."
        fi

        synthetics_activator=0
        #Check if there are synthetic reads
        if [[ -f ./synthetic/panchronos_synthetic_reads_global_statistics.tab ]]; then
                synthetics_activator=1
                export synthetics_activator
                mkdir -p ./synthetic_reads/
                awk 'NR>1 {print \$1}' ./synthetic/panchronos_synthetic_reads_global_statistics.tab > synthetic_samples.txt
                while read -r sample;do
                        if [[ -f "user_genes/\$sample".fasta ]]; then
                                mv user_genes/"\$sample".fasta ./synthetic_reads/
                        else
                                echo -e "File \$sample was not found in user_genes"
                        fi
                done < synthetic_samples.txt
        fi

        #mask genes that should not be present
        mkdir -p genes_to_mask
        mv *presence_absence_genes.index genes_to_mask

        genes_2_mask() {
        file=\$1

                name=\$(basename "\${file%_presence_absence_genes.index}")
                name="\${name#mergedGroup}"

                awk 'NR>1 && \$2 == 0 {print \$1}' "\$file" > "\$name"_tmp && mv "\$name"_tmp ./genes_to_mask/"\${name}_presence_absence_genes.index"

        }
        export -f genes_2_mask
        find ./genes_to_mask -name "*_presence_absence_genes.index" | parallel -j $parallel genes_2_mask

        mask_genes() {
        user_gene_seqs=\$1
        synth=\$2

        name=\$(basename "\${user_gene_seqs}")
        name="\${name%.fasta}"
        cleaner_name="\${name#extractedSequences}"


                if [[ "\$synth" -eq 1 ]]; then
                        mask_char="-"
                else
                        mask_char="N"
                fi

                awk -v mask="\$mask_char" 'FNR==NR{ #list of genes to mask
                        genes[\$1] = 1
                        next
                }

                #second file
                /^>/  {


                    print # print header
                    gene = substr(\$0,2) #strip > from header for regex matching
                    getline seq #store sequence for header

                    if (gene in genes) {
                        # mask sequence
                        gsub(/./, mask, seq)
                        print seq
                    } else {
                        # print sequence as-is
                        print seq
                    }

                    next
                }' ./genes_to_mask/"\${cleaner_name}_presence_absence_genes.index" "\$user_gene_seqs" > "\${name}"_tmp_masked && mv "\${name}"_tmp_masked "\${user_gene_seqs}"
        }
        export -f mask_genes

        if [[ "\$synthetics_activator" -eq 1 && "\$input_as_modern_activator" -eq 1 ]]; then
                find ./input_modern/ -name "*fasta" | parallel -j $parallel mask_genes {} 1
                find ./synthetic_reads/ -name "*fasta" | parallel -j $parallel mask_genes {} 1
                find ./user_genes/ -name "*fasta" | parallel -j $parallel mask_genes {} 0
        elif [[ "\$synthetics_activator" -eq 1 && "\$input_as_modern_activator" -eq 0 ]]; then
                find ./synthetic_reads/ -name "*fasta" | parallel -j $parallel mask_genes {} 1
                find ./user_genes/ -name "*fasta" | parallel -j $parallel mask_genes {} 0
        elif [[ "\$synthetics_activator" -eq 0 && "\$input_as_modern_activator" -eq 1 ]]; then
                find ./user_genes/ -name "*fasta" | parallel -j $parallel mask_genes {} 0
                find ./input_modern/ -name "*fasta" | parallel -j $parallel mask_genes {} 1
        else
                find ./user_genes/ -name "*fasta" | parallel -j $parallel mask_genes {} 0
        fi

        # modern_samples_list() will create a text file with modern genomes names
        modern_samples_list() {
        fasta=\$1
                name=\$(basename "\${fasta%.fasta}")
                echo "\${name}" > "\${name}"_modern_TMP
        }
        export -f modern_samples_list
        find modern_data/ -name "*.fasta" | parallel -j $parallel modern_samples_list

        cat *_modern_TMP >> modernSampleNames.txt
        rm *_modern_TMP


        # Adding the outgroup to this as it is modern too
        echo outgroup >> modernSampleNames.txt
        echo -e "Standardise fasta suffixes"

        panaroo_fasta_suffix() {
        fasta_file=\$1
        filename=\$(basename "\${fasta_file%.aln.fas}")

                mv "\${fasta_file}" panaroo_genes/"\${filename}.fasta"
        }
        export -f panaroo_fasta_suffix
        find panaroo_genes/ -name "*.aln.fas" | parallel -j $parallel panaroo_fasta_suffix

        echo -e "Done"

        echo -e "Fixing FASTA headers and sequences formatting with seqtk in existing gene alignments"

        mkdir -p panaroo_parsed
        parsing_panaroo_msa() {
        fasta_file=\$1
                name=\$(basename "\${fasta_file%.fasta}" | sed -e 's/~/_/g') # Also replace  ~ characters with _, with double % to remove two dots
                seqtk seq "\${fasta_file}" | awk '/^>/ {sub(/;.*/, "", \$0)} {print}' > panaroo_parsed/"\${name}_parsing_panaroo.fasta"
        }
        export -f parsing_panaroo_msa
        find panaroo_genes/ -name "*.fasta" | parallel -j $parallel parsing_panaroo_msa

        echo -e "Done"

        #index_and_formatting() send user samplenames to userSampleNames.txt and replaces ~ from gene names with _
        index_and_formatting() {
        sample=\$1
        synth=\$2

                sampleName=\$(basename "\${sample%.fasta}")
                if [[ "\$synth" -eq 0 ]]; then
                        echo "\${sampleName}" >> "\${sampleName}"_userSampleNames_tmp.txt
                fi

                sed -i -e 's/~/_/g' "\${sample}"
        }
        export -f index_and_formatting


        if [[ "\$synthetics_activator" -eq 1 && "\$input_as_modern_activator" -eq 1 ]]; then
                find synthetic_reads/ -name "*.fasta" | parallel -j $parallel index_and_formatting {} 1
                find input_modern/ -name "*fasta" | parallel -j $parallel index_and_formatting {} 1
                find user_genes/ -name "*.fasta" | parallel -j $parallel index_and_formatting {} 0
        elif [[ "\$synthetics_activator" -eq 1 && "\$input_as_modern_activator" -eq 0 ]]; then
                find synthetic_reads/ -name "*fasta" | parallel -j $parallel index_and_formatting {} 1
                find user_genes/ -name "*fasta" | parallel -j $parallel index_and_formatting {} 0
        elif [[ "\$synthetics_activator" -eq 0 && "\$input_as_modern_activator" -eq 1 ]]; then
                find user_genes/ -name "*fasta" | parallel -j $parallel index_and_formatting {} 0
                find input_modern/ -name "*fasta" | parallel -j $parallel index_and_formatting {} 1
        else
                find user_genes/ -name "*.fasta" | parallel -j $parallel index_and_formatting {} 0
        fi

        #concatenating userSampleNames if synth eq 0
        if compgen -G "*_userSampleNames_tmp.txt" > /dev/null; then
                cat *_userSampleNames_tmp.txt >> userSampleNames.txt
                rm *_userSampleNames_tmp.txt
        fi

		#before making the files with filenames I need to exclude blacklisted samples
        grep -vxFf blackListed.txt userSampleNames.txt | sed '/^\$/d' > tmp1 && mv tmp1 userSampleNames.txt      
        grep -vxFf blackListed.txt modern_group_as_input.txt | sed '/^\$/d' > tmp2 && mv tmp2 modern_group_as_input.txt      

	
        #make a file with every sample combined
        cat modernSampleNames.txt userSampleNames.txt > sampleNames.txt

        if [[ -s modern_group_as_input.txt ]]; then
                cat modern_group_as_input.txt >> sampleNames.txt
                cat modernSampleNames.txt modern_group_as_input.txt > panaroo_modern_and_input_modern.txt
        fi

        #including synthetic samples
        if [[ -s synthetic_samples.txt ]]; then
                cat userSampleNames.txt synthetic_samples.txt > sample_names_withSynthetics.txt

                if [[ -s modern_group_as_input.txt ]]; then
                        cat modern_group_as_input.txt >> sample_names_withSynthetics.txt
                fi

                echo "outgroup" >> sample_names_withSynthetics.txt
        fi

        echo -e "Adding user sample genes sequences to each particular gene MSA and replace gene name with sample name"

        #add_user_sample_sequences() adds user sample gene sequences into each panaroo msa and replaces gene name with user sample name
        add_user_sample_sequences() {
        fasta_file=\$1
        modern_switch=\$2


                if [[ -e "\$fasta_file" ]]; then
                        name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")

                        while read -r sampleName; do
                                grep -w -A 1 "\$name" "user_genes/\${sampleName}.fasta" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${fasta_file}"
                        done < userSampleNames.txt

                        if [[ \$modern_switch -eq 1 ]]; then

                                while read -r sampleName; do
                                        grep -w -A 1 "\$name" "input_modern/\${sampleName}.fasta" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${fasta_file}"
                                done < modern_group_as_input.txt

                        fi
                else
                        echo -e "There are no files with .fasta extension in panaroo_parsed/ folder. Stopping the process.\n"
                        exit 1
                fi
        }
        export -f add_user_sample_sequences
        if [[ "\$input_as_modern_activator" -eq 1 ]]; then
                find panaroo_parsed/ -name "*.fasta" | parallel -j $parallel add_user_sample_sequences {} 1
        else
                find panaroo_parsed/ -name "*.fasta" | parallel -j $parallel add_user_sample_sequences {} 0
        fi
        echo -e "Done"


        #add_synthetic_sample_sequences() adds synthetic samples gene sequences into each panaroo msa and replaces gene name with user sample name
        add_synthetic_sample_sequences() {
        fasta_file=\$1

                if [[ -e "\$fasta_file" ]]; then
                        name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")

                        while read -r sampleName; do
                                grep -w -A 1 "\$name" "synthetic_reads/\${sampleName}.fasta" | awk -v newHeader="\$sampleName" '/^>/ {sub(/^>.*/, ">" newHeader, \$0)} {print}' >> "\${fasta_file}"
                        done < synthetic_samples.txt
                else
                        echo -e "There are no files with .fasta extension in panaroo_parsed/ folder. Stopping the process.\n"
                        exit 1
                fi
        }
        export -f add_synthetic_sample_sequences
        if [[ "\$synthetics_activator" -eq 1 ]]; then
                find panaroo_parsed/ -name "*.fasta" | parallel -j $parallel add_synthetic_sample_sequences
        fi
        echo -e "Done"




        echo -e "Adding outgroup gene sequences into each gene MSA file"
        #similarly as add_user_sample_sequences(), add_outgroup_sequences() will add outgroup aligned gene sequences into each gene panaroo msa file
        add_outgroup_sequences() {
        fasta_file=\$1

                name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
                grep -w -A 1 "\$name" outgroup | awk -v outgroup="outgroup" '/^>/ {sub(/^>.*/, ">" outgroup, \$0)} {print}' >> "\${fasta_file}"
        }
        export -f add_outgroup_sequences
        find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel add_outgroup_sequences

        echo -e "Done"


        echo -e "Padding incomplete sequences to the right"

        panaroo_parsed_padding() {
        MSA=\$1
                name=\$(basename "\${MSA}")
                seqLength=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$MSA")

                awk -v seq_length="\$seqLength" '{
                        if (\$0 ~ /^>/) {
                                print
                        } else {
                        while ( length(\$0) < seq_length) {
                                \$0 = \$0 "n"
                        }
                        print
                        }
                }' "\${MSA}" > tmp_"\${name}" && mv tmp_"\${name}" "\${MSA}"
        }
        export -f panaroo_parsed_padding
        find panaroo_parsed/ -name "*parsing_panaroo.fasta" | parallel -j $parallel panaroo_parsed_padding

        echo -e "Done"


        # here, if synthetic reads were used, then it's time to remove the panaroo input sequences from each gene MSA before finding artifacts.
        # samples in sample_names_withSynthetics.txt needs to be kept and the rest dumped for file in panaroo_parsed/
        removing_input_panaroo() {
        fasta_file=\$1

                touch "\$fasta_file"_synth_tmp
                while read -r sample_to_keep;do
                        grep -A 1 -w "\$sample_to_keep" "\$fasta_file" >> "\$fasta_file"_synth_tmp
                done < sample_names_withSynthetics.txt

                mv "\$fasta_file"_synth_tmp "\$fasta_file"
        }
        export -f removing_input_panaroo
        if [[ "\$synthetics_activator" -eq 1 ]]; then
                find ./panaroo_parsed -name "*fasta" | parallel -j  $parallel removing_input_panaroo
        fi


        echo -e "Checking if there are Panaroo headers artifacts"

        check_panaroo_headers_artifacts() {
        file_msa=\$1

                while read -r sampleName; do
                        header_name=">\$sampleName"
                        matched_header=\$(grep "\$sampleName" "\$file_msa")

                        if [[ -n "\$matched_header" ]]; then

                                while IFS= read -r matched; do

                                        if [[ "\$matched" != "\$header_name" ]]; then
                                                sed -i -e "s/\${matched}/\${header_name}/g" "\$file_msa"
                                                echo -e "\$matched name was found in \$file_msa but it should have been \$header_name instead. Fixed"
                                        fi
                                done <<< "\$matched_header"
                        fi

                done < modernSampleNames.txt
        }
        export -f check_panaroo_headers_artifacts
        if [[ "\$synthetics_activator" -eq 0 ]]; then
                find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel check_panaroo_headers_artifacts
        fi
        echo -e "Done"



        modern_append_gaps() {
        fasta_file=\$1
        index_file=\$2

                seq_length=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$fasta_file")

                while read -r strain; do
                        if ! grep -wq "\$strain" "\$fasta_file"; then
                                echo ">\$strain" >> "\$fasta_file"
                                gaps_length=\$(printf '%*s' "\$((seq_length))" | tr ' ' '-')
                                echo "\$gaps_length" >> "\$fasta_file"
                        fi
                done < "\${index_file}"
        }
        export -f modern_append_gaps
        if [[ "\$synthetics_activator" -eq 0 ]]; then
                if [[ "\$input_as_modern_activator" -eq 1 ]]; then
                        find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel modern_append_gaps {} modernSampleNames.txt
                else
                        find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel modern_append_gaps {} panaroo_modern_and_input_modern.txt
                fi
        fi


        #synth_append_gaps() will check if there are modern strains from synthetic reads missing in panaroo_parsed alignments. If yes, append the sample and padding with -
        echo -e "Padding missing samples. Gaps for modern strains and n for ancient samples"

        synth_append_gaps() {
        fasta_file=\$1

                seq_length=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$fasta_file")

                while read -r strain; do
                        if ! grep -wq "\$strain" "\$fasta_file"; then
                                echo ">\$strain" >> "\$fasta_file"
                                gaps_length=\$(printf '%*s' "\$((seq_length))" | tr ' ' '-')
                                echo "\$gaps_length" >> "\$fasta_file"
                        fi
                done < synthetic_samples.txt
        }
        export -f synth_append_gaps
        if [[ "\$synthetics_activator" -eq 1 ]]; then
                find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel synth_append_gaps
        fi


        #ancient_append_missingness() will check if there are ancient strains missing in panaroo_parsed alignments. If yes, append the sample and fill it with n as we don't know the ancient status for that gene.
        ancient_append_missingness() {
        fasta_file=\$1

                seq_length=\$(awk '!/^>/ { if (length > max) max = length } END { print max }' "\$fasta_file")

                while read -r strain; do
                        if ! grep -wq "\$strain" "\$fasta_file"; then
                                echo ">\$strain" >> "\$fasta_file"
                                fakeSeq=\$(printf '%*s' "\$((seq_length))" | tr ' ' 'n')
                                echo "\$fakeSeq" >> "\$fasta_file"
                        fi
                done < userSampleNames.txt
        }

        export -f ancient_append_missingness
        find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel ancient_append_missingness

        echo -e "Done"


        # At this point there will be many msa that are finished: number of headers == number of samples AND headers names == all samples names. I need to find them based on this logic and send them
        # to the final folder. There will be some msa with header duplication problems and they need to be send to special_cases.

        mkdir -p sanitised_msa

        sanitised_msa() {
        msa_file=\$1
        synth=\$2

                name=\$(basename "\${msa_file%_parsing_panaroo.fasta}.fasta")
                sample_count=\$(grep -c '^>' "\$msa_file")

                if [[ "\$synth" -eq 0 ]]; then
                        input_file=sampleNames.txt
                else
                        input_file=sample_names_withSynthetics.txt
                fi

                total_sample_count=\$(wc -l < "\$input_file")

                if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
                        echo "Samples from index and MSA do not match"
                        echo "Samples in msa file"
                        cat \$msa_file
                        echo "Samples in index file"
                        cat \$input_file
                        return
                fi



                missing_samples=0
                while read -r strain; do
                        if ! grep -wq "\$strain" "\$msa_file"; then
                                missing_samples=1
                                break
                        fi
                done < "\$input_file"

                if [[ "\$missing_samples" -eq 0 ]]; then
                        mv "\$msa_file" "sanitised_msa/\${name}"
                fi

        }
        export -f sanitised_msa
        find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel sanitised_msa {} "\$synthetics_activator"



        # Dealing with fragmented Panaroo gene alignments multi-entries
        mkdir -p special_cases

        echo -e "Checking if there are leftovers after sanitising"

        unsanitised=\$(ls panaroo_parsed/*fasta | wc -l )
        if (( unsanitised != 0 )); then
                echo -e "There are \$unsanitised samples that have problems in panaroo_parsed/ directory. Inspecting them. ."
        else
                echo -e "It seems every gene msa was succesfully cleaned. Moving on."
        fi

        fix_duplicated_entries() {
        fasta_file=\$1

                if [[ -e "\$fasta_file" ]] ; then
                        name=\$(basename "\${fasta_file%_parsing_panaroo.fasta}")
                        # Identifying repeated headers and save them to .dup
                        awk '/^>/ {count[\$0]++} END {for (header in count) if (count[header] > 1) print substr(header, 2)}' "\$fasta_file" > special_cases/"\${name}.dup"

                        # Extracting sequences for each repeated entry into temporary files
                        while read -r entry; do
                                awk -v sampleName="\$entry" '
                                \$0 ~ "^>" sampleName {print_header=1; next}
                                /^>/ {print_header=0}
                                print_header {print}' "\$fasta_file" > special_cases/"\${name}_duplicated_\${entry}"
                       done < special_cases/"\${name}.dup"
               # Remove repeated entries and their sequences from the original FASTA file
                       awk -v dpdFile=special_cases/"\${name}.dup" '
                       BEGIN {
                       while (getline < dpdFile) {
                               repeated[\$0] = 1
                               }
                       }
                       /^>/ { header = substr(\$0, 2)
                       if (repeated[header]) {
                               skip = 1
                       } else {
                               skip = 0}
                       }
                       !skip' "\$fasta_file" > special_cases/"\${name}.fasta"

                       for indexSeqs in special_cases/"\${name}_duplicated_"*; do
                               geneName=\$(basename "\${indexSeqs%_duplicated_*}")
                               sampleName=\$(basename "\${indexSeqs##*_duplicated_}")

                              sequence=\$(awk '
                                {
                                    seq_length = length(\$0)
                                    line[NR] = \$0
                                    num_lines = NR
                                }

                                END {
                                    final = ""
                                    for (pos = 1; pos <= seq_length; pos++) {
                                        merged_char = "-"
                                        for (l = 1; l <= num_lines; l++) {
                                            c = substr(line[l], pos, 1)
                                            if (c ~ /[aAtTgGcC]/) {
                                                merged_char = c
                                                break
                                            }
                                        }
                                        final = final merged_char
                                    }

                                    # Remove dashes
                                    gsub(/-/, "", final)
                                    print final
                                }' "\$indexSeqs")



                       # Finally add the selected sequence back to the cleaned original FASTA file with the header as well
                               echo ">\${sampleName}" >> special_cases/"\${name}.fasta"
                               echo "\$sequence" >> special_cases/"\${name}.fasta"
                       done

               else
                       echo -e "No files found with fasta extension in special_cases folder"
                # Cleaning temporary files
               fi
       }

        export -f fix_duplicated_entries

        if (( unsanitised != 0 )); then
                find panaroo_parsed/ -name "*_parsing_panaroo.fasta" | parallel -j $parallel fix_duplicated_entries
        fi

        #finally checking that every seq has the same length

        checking_alignment_lengths() {
        msa_file=\$1
        name=\$(basename "\${msa_file}")

               seq_length=\$(awk 'NR%2 == 0 && length > max { max = length } END { print max }' "\$msa_file")

               awk -v numCols="\$seq_length" '{
                       if (\$0 ~ /^>/) {
                               print
                       } else {
                       while ( length(\$0) < numCols) {
                               \$0 = \$0 "-"
                       }
                       print
                       }
               }' "\$msa_file" > special_cases/tmp_"\${name}" && mv special_cases/tmp_"\${name}" "\${msa_file}"
       }

        #export -f checking_alignment_lengths
        #find special_cases/ -name "*.fasta" | parallel -j $parallel checking_alignment_lengths
        #echo -e "Done"
        #NOTE: We may not need this function anymore as we are re-aligning everything with mafft.


        #now check if these cleaned fasta passed the sanitised_msa() test.
        sanitised_msa_special_cases() {
        msa_file=\$1
        synth=\$2

                name=\$(basename "\${msa_file}")
                sample_count=\$(grep -c '^>' "\$msa_file")
                if [[ "\$synth" -eq 0 ]]; then
                        input_file=sampleNames.txt
                else
                        input_file=sample_names_withSynthetics.txt
                fi

                total_sample_count=\$(wc -l < "\$input_file")


                if [[ "\$sample_count" -ne "\$total_sample_count" ]]; then  #if they are equal then keep going. if not skip current file.
                        return
                fi

                missing_samples=0
                while read -r strain; do
                        if ! grep -wq "\$strain" "\$msa_file"; then
                                missing_samples=1
                                break
                        fi
                done < "\$input_file"

                if [[ "\$missing_samples" -eq 0 ]]; then
                        mv "\$msa_file" "sanitised_msa/\${name}"
                fi

        }
        export -f sanitised_msa_special_cases
        find special_cases/ -name "*.fasta" | parallel -j $parallel sanitised_msa_special_cases {} "\$synthetics_activator"

        # final steps: sort everything.
        # Make INDEX file so we can sort every file in the same order before building MSA based on first file.
        firstFile=\$(ls -1 sanitised_msa/ | awk 'NR==1 {print \$0}') && awk '/^>/ {print \$0}' sanitised_msa/"\$firstFile" > sorting_index

        mkdir -p sorted

        echo -e "Sorting MSAs"

        sortingMSA() {
        fasta=\$1

                name=\$(basename "\${fasta%.fasta}")

                while read -r header; do
                        awk -v headerName="\$header" '
                                \$0 == headerName {
                                print \$0     # Print the matched header
                                getline      # Get the sequence line after the header and store it in line
                                print \$0     # Print the sequene
                                }
                        ' "\$fasta" >> sorted/"\${name}.fasta"
                done < sorting_index
        }

        export -f sortingMSA
        find sanitised_msa/ -name "*.fasta" | parallel -j $parallel sortingMSA

        echo -e "Done"

        cat .command.out >> filterGeneAlignments.log
        """
}
