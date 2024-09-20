#!/bin/bash

#Values
data=""
output=""
lcompleteness=50
coverage=0.5
threads=10
tax_id=""
genomes=100
clustering=0.95
core=0.01
clean="strict"
config=""


# Function to print help
helpf() {
echo -e "\n\e[1;31mSYNOPSIS\e[0m"
        echo -e "\n\e[1;33mUSAGE\e[0m \n\n$ bash normalization.sh -d <sample> -o <OUTPUT PATH> [-t <INT>]"
        echo -e "\n\e[1;33mOPTIONS\e[0m"
        echo -e "\nMandatory"
        echo -e "  -d, --data <PATH>			Set data file PATH"
        echo -e "  -o, --output <PATH>			Set output directory PATH"
	echo -e "  -n, --taxid <INT>			Set taxonomical ID value <INT>"
	echo -e "  -f, --config <PATH>                  Set config file PATH"

        echo -e "\nOptional"
        echo -e "  -t, --threads <INT>			Set number of threads (default: 10)"
	echo -e "  -b, --completeness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
	echo -e "  -C, --coverage <INT/FLOAT>		Set mean depth of coverage threshold (default: 0.5)"
	echo -e "  -g, --genomes <INT> 			Set number of genomes to download (default: 100)"
	echo -e "  -c, --clustering <INT/FLOAT>		Set clustering threshold (default 0.95)"
	echo -e "  -p, --core-threshold <FLOAT>		Set core genome threshold (default: 0.01)"	
 	echo -e "  -m, --clean-mode <STRING>		Set pangenome mode (default: strict)"
        echo -e "  -h, --help				Print this help message and exit."
	
echo -e "\n\e[1;31mDESCRIPTION\e[0m"
	echo -e "\n\e[1;33m-d, --data <PATH>\e[0m:\nPlease specify the full PATH of your data. Example: /home/user/mydata/data"
 	echo -e "\n\e[1;33m-o, --output <DIR>\e[0m:\nPlease specify the full PATH of your output folder. You need to make the folder first before running the program."
	exit 1
}


while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -d|--data) data="$2"; shift ;;
        -t|--threads) threads="$2"; shift ;;
        -o|--output) output="$2"; shift ;;
	-b|--lcompleteness) lcompleteness="$2"; shift ;;
 	-C|--coverage) coverage="$2"; shift ;;
   	-n|--taxid) tax_id="$2"; shift ;;
    	-g|--genomes) genomes="$2"; shift ;;
	-c|--clustering) clustering="$2"; shift ;;
 	-p|--core-threshold) core="$2"; shift ;;
  	-m|--clean-mode) clean="$2"; shift ;;
        -f|--config) config="$2"; shift ;;
        -h|--help) helpf ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done



echo -e "Setting up directory structure\n"

mkdir -p "$output"/NCBI/FASTA
mkdir -p "$output"/NCBI/GFF
mkdir -p "$output"/CLUSTERING
mkdir -p "$output"/PROKKA
mkdir -p "$output"/PROKKA/GFF
mkdir -p "$output"/PANGENOME
mkdir -p "$output"/ALIGNMENTS
mkdir -p "$output"/NORMALIZATION
mkdir -p "$output"/MATRIX
mkdir -p "$output"/PLOTS
mkdir -p "$output"/HETEROPLASMY
mkdir -p "$output"/HETEROPLASMY/intermediate_files
mkdir -p "$output"/HETEROPLASMY/distributions

echo -e "Done\n"

#PLEASE NOTE THAT EVERY ENVIRONMENT FILE IS IN A FOLDER CALLED envs/
#THIS STEP USES ENVIRONMENT CALLED entrez.yaml

echo -e "Downloading GenBank and FASTA files based on taxonomic ID\n"
counter=0
#get GenBank and FASTA files based on taxonomic ID and stop when reaching genomes variable
esearch -db assembly -query "txid${tax_id}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank | while read -r url; do
  if [ "$counter" -ge "$genomes" ]; then
    break
  fi
  
  fname=$(basename "$url")
  
  wget -P "$output/NCBI/GFF" "$url/${fname}_genomic.gbff.gz"
  wget -P "$output/NCBI/FASTA" "$url/${fname}_genomic.fna.gz"
  counter=$((counter + 1))
  echo -e "Strains downloaded: $counter of $genomes\n\n"
done

echo -e "Done\n"


#THIS STEP USES ENVIRONMENT CALLED biopython.yaml
echo -e "Parsing and building FASTA database"

gzip -d "$output"/NCBI/GFF/*
gzip -d "$output"/NCBI/FASTA/*

#Rename sequences because people want better names

for fasta in "$output"/NCBI/FASTA/*; do

        nameFASTA=$(basename "${fasta}")
        whole=$(head -n 1 "$fasta" | awk -F',' '{print $1}' | sed -e 's/chromosome//g' | sed -e 's/://g' | awk '{$1=""; sub(/^ /, "") ;print $0}' | sed -e 's/ /_/g')
        mv "$output"/NCBI/FASTA/"$nameFASTA" "$output"/NCBI/FASTA/"${whole}.fna"

done


for gff in "$output"/NCBI/GFF/*; do

        nameGFF=$(basename "${gff}")
        wholegff=$(awk 'NR==2 {$1=""; sub(/^ /, "");  print $0}' "$gff" | awk -F',' '{print $1}' | sed -e 's/chromosome//g' -e 's/://g' -e 's/ /_/g' -e 's/_$//' -e 's/_.$//g')
        mv "$output"/NCBI/GFF/"$nameGFF" "$output"/NCBI/GFF/"${wholegff}.gbff"
done

python parsing_and_contatenating.py "$output"/NCBI/GFF

mv clustered_sequences.fasta "$output"/CLUSTERING/

echo -e "Done\n"

#THIS STEP USES ENVIRONMENT CALLED cdhit.yaml
echo -e "Clustering gene sequences"

cd-hit-est -i "$output"/CLUSTERING/clustered_sequences.fasta -o "$output"/CLUSTERING/clustered_non_redundant_genes.fasta -c "$clustering" -n 10 -M 0 -T "$threads"

echo -e "Done\n"

#THIS STEP USES ENVIRONMENT CALLED prokka.yaml

echo -e "Annotating FASTA sequences\n"

ls -l "$output"/NCBI/GFF/*gbff | awk 'NR==1{print $NF}' > "$output"/NCBI/GFF/first
name=$(cat "$output"/NCBI/GFF/first | awk -F'/' '{print $NF}')
echo -e "$name"
species=$(head -n 20 "$output"/NCBI/GFF/"$name" | grep "ORGANISM" | awk '{print $2, $3}' | sed -e 's/ /_/g')
echo -e "$species"
for i in "$output"/NCBI/FASTA/*; do
        name=$(basename "$i")
        prokka --outdir "$output"/PROKKA/"${name%.fna}" --addgenes  --addmrna --species "$species" --proteins "$output"/CLUSTERING/clustered_non_redundant_genes.fasta --force --cpus "$threads" "$i"
done


for sample in "$output"/PROKKA/*; do
        name=$(basename "$sample")
        mv "$sample"/*gff "$output"/PROKKA/GFF/"${name}.gff"
done



echo -e "Done\n"

#THIS STEP USES ENVIRONMENT CALLED panaroo.yaml

echo -e "Building pangenome from annotated files\n"

panaroo -i "$output"/PROKKA/GFF/*.gff -o "$output"/PANGENOME/ --clean-mode "$clean" -a core --core_threshold "$core" -t "$threads"

echo -e "Done\n"

#THIS STEP USES ENVIRONMENT CALLED alignment.yaml

echo -e "Formatting pangenome sequence"

seqtk seq "$output"/PANGENOME/pan_genome_reference.fa > "$output"/PANGENOME/pan_genome_reference.fasta

echo -e "Done\n"

echo -e "Aligning samples against pangenome reference"

for sample in "$data"/*; do

        echo -e "\nGenerating the reference genome index files . . ."
        bwa index "$output"/PANGENOME/pan_genome_reference.fasta
        name=$(basename "$sample")
        softclipping=$(grep "$name" "$config" | awk '{print $2}')
        echo -e "You have selected ${softclipping} as the number of bases to be trimmed for ${name}"

	echo -e "\nRunning alignment against reference . . ."
	bwa aln -l 16500 -n 0.01 -o 2 -t "$threads" "$output"/PANGENOME/pan_genome_reference.fasta "$sample" > "$output"/ALIGNMENTS/"${name%.fastq*}.sai"

	echo -e "\nConverting SAI to SAM . . ."
	bwa samse "$output"/PANGENOME/pan_genome_reference.fasta "$output"/ALIGNMENTS/"${name%.fastq*}.sai" "$sample" > "$output"/ALIGNMENTS/"${name%.fastq*}.sam"

	echo -e "\nConverting SAM to BAM and sorting . . ."
	samtools view -bS "$output"/ALIGNMENTS/"${name%.fastq*}.sam" > "$output"/ALIGNMENTS/"${name%.fastq*}.bam"

	echo -e "\nChecking sanity of BAM file . . ."
	samtools quickcheck "$output"/ALIGNMENTS/"${name%.fastq*}.bam"

	echo -e "\nSorting BAM . . ."
	samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}.bam"

	echo -e "\nGenerating BAM index . . ."
	samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam"

	echo -e "\nGetting only mapped reads . . . "
	samtools view -b -@ 10 -F 4 "$output"/ALIGNMENTS/"${name%.fastq*}_sorted.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"
	samtools index "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam"

	echo -e "\nRemoving 5 bases at each end of every reads . . ."
	~/miniforge3/envs/alignment/bin/bam trimBam "$output"/ALIGNMENTS/"${name%.fastq*}_sorted_mappedreads.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam" -L "$softclipping" -R "$softclipping" --clip

	echo -e "\nGetting only reads with 25 mapping quality or more . . ."
	samtools view -q 25 -o "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_softclipped.bam"

	echo -e "\nRemoving reads smaller than 34bp . . ."
	samtools view -e 'length(seq)>34' -O BAM -o "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam" "$output"/ALIGNMENTS/"${name%.fastq*}_qc.bam"

	echo -e "\nSorting the output . . ."
	samtools sort -o "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" -O bam -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_lg.bam"

	echo -e "\nComputing basic statistics . . ."
	samtools coverage "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name}"_genomicsMetrics.txt

	echo -e "\nConverting to FASTQ . . ."
	samtools fastq -@ "$threads" "$output"/ALIGNMENTS/"${name%.fastq*}_DMC_P.bam" > "$output"/ALIGNMENTS/"${name%.fastq*}_final.fastq"

	rm "$output"/ALIGNMENTS/*_softclipped.bam "$output"/ALIGNMENTS/*_qc.bam "$output"/ALIGNMENTS/*_lg.bam "$output"/ALIGNMENTS/*sai "$output"/ALIGNMENTS/*sam "$output"/ALIGNMENTS/"${name%.fastq*}.bam"
	
done

echo -e "Done\n"



#THIS PROCESS IS STILL IN DEVELOPTMENT, PLEASE DO NOT INCLUDE

: ' THIS NEEDS FIXING
for sample in "$data"/*; do
	samtools mpileup
	awk '{print $2,$4,$5,$6,$7}' "$output"/HETEROPLASMY/"$file" > "$output"/HETEROPLASMY/"${file}_filtered_pileup"
	awk '$3 ~ /[acgtACTG]/' "$output"/HETEROPLASMY/"${file}_filtered_pileup" > "$output"/HETEROPLASMY/"${file}_het1"
	sed 's/\([atcg]\)/\U\1/g' "$output"/HETEROPLASMY/"${file}_het1" > "$output"/HETEROPLASMY/"${file}_het2"
	sed 's/,/./g' "$output"/HETEROPLASMY/"${file}_het2" > "$output"/HETEROPLASMY/"${file}_het3"
	sed 's/\^F//g' "$output"/HETEROPLASMY/"${file}_het3" > "$output"/HETEROPLASMY/"${file}_het4"
	sed 's/\$//g' "$output"/HETEROPLASMY/"${file}_het4" > "$output"/HETEROPLASMY/"${file}_het5"
	awk '$3 !~ /-[0-9]+/ && $3 !~ /\+[0-9]+/' "$output"/HETEROPLASMY/"${file}_het5" > "$output"/HETEROPLASMY/"${file}_het6"
	#TODO WHAT THE HECK ARE THESE LINES
	awk '$3 ~ /-[0-9]+/ || $3 ~ /\+[0-9]+/' "$output"/HETEROPLASMY/"${file}_het5" > "$output"/HETEROPLASMY/"${file}_weird_lines"
	awk '$2 != 1' "$output"/HETEROPLASMY/"${file}_het6" > "$output"/HETEROPLASMY/"${file}_het7"
	awk '{ dot_count = gsub(/\./, "", $3); $5 = dot_count - $2; print }' "$output"/HETEROPLASMY/"${file}_het7" > "$output"/HETEROPLASMY/"${file}_het8"
	awk '{print $5}' "${file}_het8" | paste "$output"/HETEROPLASMY/"${file}_het7" - > "$output"/HETEROPLASMY/"${file}_het9"
	awk '{ $5 = ($5 < 0) ? -$5 : $5; print }' "$output"/HETEROPLASMY/"${file}_het9" > "$output"/HETEROPLASMY/"${file}_het10"
	awk '$5 != $2' "$output"/HETEROPLASMY/"${file}_het10" > "$output"/HETEROPLASMY/"${file}_het11"
	awk '{new_column = ($5 / $2) * 100; print $0, new_column}' "$output"/HETEROPLASMY/"${file}_het11" > "$output"/HETEROPLASMY/heterozygosis.txt
	awk '$6<=10' "$output"/HETEROPLASMY/heterozygosis.txt | wc -l > "$output"/HETEROPLASMY/less_than_or_equal_to_10%
	awk '$6>=90' "$output"/HETEROPLASMY/heterozygosis.txt | wc -l > "$output"/HETEROPLASMY/bigger_than_or_equal_to_90%
#Removing sites where there is only 1 read from the main file (OPTIONAL)
#awk 2 != 1 filtered_pileup > mpile_filtered.txt
	mv "$output"/HETEROPLASMY/*_than_* "$output"/HETEROPLASMY/distributions/
	mv "$output"/HETEROPLASMY/*_het* "$output"/HETEROPLASMY/intermediate_files/
done

'







#THIS STEP USES ENVIRONMENT CALLED normalization.yaml

echo -e "Extracting raw coverage per gene\n"


for i in "$output"/ALIGNMENTS/*_DMC_P.bam; do

	samplename=$(basename "${i%_DMC_P.bam}")
	samtools index "$i"
	samtools depth -a "$i" > "$output"/NORMALIZATION/"${samplename}_rawCoverage.txt"
	samtools idxstats "$i" | awk '{sum += $2} END {print sum}' > "$output"/NORMALIZATION/"${samplename}_refLength.txt"
	samtools coverage "$i" | awk -v samplename="$samplename" 'NR>1 {print samplename, $1, $6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t $'\t' >> "$output"/NORMALIZATION/completenessSummary.tab

done

echo -e "Done\n"


echo -e "Normalizing coverage per gene\n"

echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > "$output"/NORMALIZATION/geneNormalizedSummary.txt
echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > "$output"/NORMALIZATION/globalMeanCoverage.txt

for i in "$output"/NORMALIZATION/*_rawCoverage.txt; do
    name=$(basename "${i%_rawCoverage.txt}")
    echo -e "Sample being processed: $name\n"

    # Compute global mean coverage
    globalMean=$(awk -v name="$name" '{sum += $3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check " name}' "$i")
    finalCount=$(awk -v name="$name" '{fcount++} END {print fcount}' "$i")
    refCount=$(cat "$output"/NORMALIZATION/"${name}_refLength.txt")
    echo -e "$name\t$finalCount\t$refCount\t$globalMean" >> "$output"/NORMALIZATION/globalMeanCoverage.txt

    # Normalize coverage per gene
    awk -v globalMean="$globalMean" -v name="$name" -v sampleCoverage="$finalCount" -v refCount="$refCount" '
        {
            if ($2 > geneLength[$1]) {
                geneLength[$1] = $2
            }
            sumgene[$1] += $3
            countgene[$1]++
        }
        END {
            for (gene in geneLength) {
                geneMean = sumgene[gene] / countgene[gene]
                normalizedGeneScaled = (geneMean / globalMean) * geneLength[gene]
		normalizedGeneSimple = (geneMean / globalMean)
		normalizedGenomeSimple = (geneMean / globalMean) * (geneLength[gene] / sampleCoverage)
		normalizedGenomeScaled = (geneMean / globalMean) * (geneLength[gene] / refCount)
                print name"\t"gene"\t"normalizedGeneSimple"\t"normalizedGeneScaled"\t"normalizedGenomeSimple"\t"normalizedGenomeScaled
            }
        }
    ' "$i" >> "$output"/NORMALIZATION/geneNormalizedSummary.txt

done


echo -e "Done\n"

echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > "$output"/NORMALIZATION/geneNormalizedUpdated.tab

sed -i -e 's/~/_/g' "$output"/NORMALIZATION/geneNormalizedSummary.txt

awk 'NR>1{print $1"XYZ"$2, $3, $4, $5, $6}' "$output"/NORMALIZATION/geneNormalizedSummary.txt > "$output"/NORMALIZATION/TMP1

awk '{print $1"XYZ"$2, $3}' "$output"/NORMALIZATION/completenessSummary.tab > "$output"/NORMALIZATION/TMP2

while read -r ID completeness;do

	if grep -wq "${ID}" "$output"/NORMALIZATION/TMP1; then
		oldLine=$(grep -w "${ID}" "$output"/NORMALIZATION/TMP1)
		specificCompleteness=$(grep -w "${ID}" "$output"/NORMALIZATION/TMP2 | awk '{print $NF}')
		echo -e "${oldLine}\t${specificCompleteness}" >> "$output"/NORMALIZATION/geneNormalizedUpdated.tab
	fi

done < "$output"/NORMALIZATION/TMP2

sed -i -e 's/XYZ/\t/g' "$output"/NORMALIZATION/geneNormalizedUpdated.tab
sed -i -e 's/ /\t/g' "$output"/NORMALIZATION/geneNormalizedUpdated.tab


rm "$output"/NORMALIZATION/TMP1 "$output"/NORMALIZATION/TMP2 "$output"/NORMALIZATION/geneNormalizedSummary.txt "$output"/NORMALIZATION/completenessSummary.tab
mv "$output"/NORMALIZATION/geneNormalizedUpdated.tab "$output"/NORMALIZATION/geneNormalizedSummary.tab

python plot_cvg_vs_completeness.py "$output"/NORMALIZATION/geneNormalizedSummary.tab "$lcompleteness" "$coverage"


mv plotCoverage_vs_Completeness.png "$output"/PLOTS/plotCoverage_vs_Completeness.png



awk 'NR==1{print $0}' "$output"/PANGENOME/gene_presence_absence.Rtab > "$output"/MATRIX/matrix.tab
awk 'NR>1 {print $0}' "$output"/PANGENOME/gene_presence_absence.Rtab | sort -k 1 -t $'\t' >> "$output"/MATRIX/matrix.tab
awk 'NR>1 {print $1}' "$output"/PANGENOME/gene_presence_absence.Rtab | sort -k 1 -t $'\t' > "$output"/MATRIX/INDEX

awk 'NR>1 {print $1}' "$output"/NORMALIZATION/globalMeanCoverage.txt > "$output"/MATRIX/sample_names

while read -r name; do

	echo -e "Gene\tnormalizedCoverage\tcompleteness" > "$output"/MATRIX/"${name}"_index.tmp
	grep -w "$name" "$output"/NORMALIZATION/geneNormalizedSummary.tab | awk '{print $2, $3, $NF}' >> "$output"/MATRIX/"${name}"_index.tmp

done < "$output"/MATRIX/sample_names


for i in "$output"/MATRIX/*_index.tmp; do

	sed -i -e 's/ /\t/g' "$i"
	python lambda.py "$i"

done


mv *_final.csv "$output"/MATRIX/

for i in "$output"/MATRIX/*_final.csv; do

	name=$(basename "$i")
	sed  -e 's/,/\t/g' "$i" | awk 'NR>1{print $0}' > "$output"/MATRIX/"${name%_index.tmp_final.csv}"_INDEX.Z

done


for i in "$output"/MATRIX/*_INDEX.Z; do
    name=$(basename "$i")
    # Create the FINAL_INDEX file
    echo "${name%_INDEX.Z}" > "$output"/MATRIX/"${name}"_FINAL_INDEX
    
    # Process the INDEX file
    while read -r gene; do
        toprint=$(echo "$gene 0")
        if grep -wq "$gene" "$i"; then
            grep -w "$gene" "$i" >> "$output"/MATRIX/"${name}"_FINAL_INDEX
        else
            echo "$toprint" >> "$output"/MATRIX/"${name}"_FINAL_INDEX
        fi
    done < "$output"/MATRIX/INDEX
    
    # Extract the last column
    awk '{print $NF}' "$output"/MATRIX/"${name}"_FINAL_INDEX > "$output"/MATRIX/"${name}"_FINALCOLUMN

done


paste "$output"/MATRIX/matrix.tab "$output"/MATRIX/*_FINALCOLUMN > "$output"/MATRIX/final_matrix.tab


rm "$output"/MATRIX/*_INDEX.Z "$output"/MATRIX/*_final.csv "$output"/MATRIX/*_FINAL_INDEX "$output"/MATRIX/*INDEX.Z_FINALCOLUMN "$output"/MATRIX/*_index.tmp "$output"/MATRIX/INDEX "$output"/MATRIX/matrix.tab

mv "$output"/MATRIX/final_matrix.tab "$output"/MATRIX/matrix.tab

tr '\n' ' ' < "$output"/MATRIX/sample_names > "$output"/MATRIX/names_heatmap

python heatmap.py "$output"/MATRIX/matrix.tab "$output"/MATRIX/names_heatmap



mv *png "$output"/PLOTS
