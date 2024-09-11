#!/bin/bash

#Values
data=""
output=""
completeness="50"
coverage="0.5"
softclipping="5"
threads="10"
tax_id=""
genomes="100"

# Function to print help
helpf() {
echo -e "\n\e[1;31mSYNOPSIS\e[0m"
        echo -e "\n\e[1;33mUSAGE\e[0m \n\n$ bash normalization.sh -d <sample> -o <OUTPUT PATH> [-t <INT>]"
        echo -e "\n\e[1;33mOPTIONS\e[0m"
        echo -e "\nMandatory"
        echo -e "  -d, --data <PATH>			Set data file PATH"
        echo -e "  -o, --output <PATH>			Set output directory PATH"
	echo -e "  -n, --taxid <INT>			Set taxonomical ID value <INT>"
 
        echo -e "\nOptional"
        echo -e "  -t, --threads <INT>			Set number of threads (default: 10)"
	echo -e "  -b, --completeness <INT/FLOAT>	Set gene completeness/breadth of coverage threshold (default: 50)"
	echo -e "  -C, --coverage <INT/FLOAT>		Set mean depth of coverage threshold (default: 0.5)"
	echo -e "  -s, --soft-clipping <INT>		Set soft-clipping value (default: 5)"
	echo -e "  -g, --genomes <INT> 			Set number of genomes to download (default: 100)"

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
	-b|--completeness) completeness="$2"; shift ;;
 	-C|--coverage) coverage="$2"; shift ;;
  	-s|--soft-clipping) softclipping="$2"; shift ;;
   	-n|--taxid) taxid="$2"; shift ;;
    	-g|--genomes) genomes="$2"; shift ;;
        -h|--help) helpf ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done




echo -e "Extracting raw coverage per gene"

mkdir -p ./NORMALIZATION

for i in *bam; do

	samtools index "$i"
	samtools depth -a "$i" > ./NORMALIZATION/"${i%.bam}_rawCoverage.txt"
	samtools idxstats "$i" | awk '{sum += $2} END {print sum}' > ./NORMALIZATION/"${i%.bam}_refLength.txt"
	samplename=$(basename "${i%.bam}")
	samtools coverage "$i" | awk -v samplename="$samplename" 'NR>1 {print samplename, $1, $6}' | sed -e 's/~/_/g' | sed -e 's/ /\t/g' | sort -k 1 -t $'\t' >> ./NORMALIZATION/completenessSummary.tab

done


echo -e "\nDone"

echo -e "\nNormalizing coverage per gene"

echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled" > ./NORMALIZATION/geneNormalizedSummary.txt
echo -e "sampleID \t sampleCoverage \t refCount \t globalMean"  > ./NORMALIZATION/globalMeanCoverage.txt

for i in ./NORMALIZATION/*_rawCoverage.txt; do
    name=$(basename "${i%_rawCoverage.txt}")
    echo -e "Sample being processed: $name\n"

    # Compute global mean coverage
    globalMean=$(awk -v name="$name" '{sum += $3; count++} END {if (count > 0) print sum / count; else print "Something went wrong, check " name}' "$i")
    finalCount=$(awk -v name="$name" '{fcount++} END {print fcount}' "$i")
    refCount=$(cat ./NORMALIZATION/"${name}_refLength.txt")
    echo -e "$name\t$finalCount\t$refCount\t$globalMean" >> ./NORMALIZATION/globalMeanCoverage.txt

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
    ' "$i" >> ./NORMALIZATION/geneNormalizedSummary.txt

done


echo -e "\nDone"

echo -e "sampleID\tgene\tnormalizedGeneSimple\tnormalizedGeneScaled\tnormalizedGenomeSimple\tnormalizedGenomeScaled\tgeneCompleteness" > ./NORMALIZATION/geneNormalizedUpdated.tab

sed -i -e 's/~/_/g' ./NORMALIZATION/geneNormalizedSummary.txt

awk 'NR>1{print $1"XYZ"$2, $3, $4, $5, $6}' ./NORMALIZATION/geneNormalizedSummary.txt > ./NORMALIZATION/TMP1

awk '{print $1"XYZ"$2, $3}' ./NORMALIZATION/completenessSummary.tab > ./NORMALIZATION/TMP2

while read -r ID completeness;do

	if grep -wq "${ID}" ./NORMALIZATION/TMP1; then
		oldLine=$(grep -w "${ID}" ./NORMALIZATION/TMP1)
		specificCompleteness=$(grep -w "${ID}" ./NORMALIZATION/TMP2 | awk '{print $NF}')
		echo -e "${oldLine}\t${specificCompleteness}" >> ./NORMALIZATION/geneNormalizedUpdated.tab
	fi

done < ./NORMALIZATION/TMP2

sed -i -e 's/XYZ/\t/g' ./NORMALIZATION/geneNormalizedUpdated.tab
sed -i -e 's/ /\t/g' ./NORMALIZATION/geneNormalizedUpdated.tab


rm ./NORMALIZATION/TMP1 ./NORMALIZATION/TMP2 ./NORMALIZATION/geneNormalizedSummary.txt ./NORMALIZATION/completenessSummary.tab
mv ./NORMALIZATION/geneNormalizedUpdated.tab ./NORMALIZATION/geneNormalizedSummary.tab

python plot_cvg_vs_completeness.py ./NORMALIZATION/geneNormalizedSummary.tab

mkdir -p plots

mv plotCoverage_vs_Completeness.png ./plots/plotCoverage_vs_Completeness.png

#Managing the matrix

mkdir -p matrix

awk 'NR==1{print $0}' gene_presence_absence.Rtab > ./matrix/matrix.tab
awk 'NR>1 {print $0}' gene_presence_absence.Rtab | sort -k 1 -t $'\t' >> ./matrix/matrix.tab
awk 'NR>1 {print $1}' gene_presence_absence.Rtab | sort -k 1 -t $'\t' > ./matrix/INDEX

awk 'NR>1 {print $1}' ./NORMALIZATION/globalMeanCoverage.txt > ./matrix/sample_names

while read -r name;do

	echo -e "Gene\tnormalizedCoverage\tcompleteness" > ./matrix/"${name}"_index.tmp
	grep -w "$name" ./NORMALIZATION/geneNormalizedSummary.tab | awk '{print $2, $3, $NF}' >> ./matrix/"${name}"_index.tmp

done < ./matrix/sample_names


for i in ./matrix/*_index.tmp; do

	sed -i -e 's/ /\t/g' "$i"
	python lambda.py "$i"

done


mv *_final.csv ./matrix/

for i in ./matrix/*_final.csv; do

	name=$(basename "$i")
	sed  -e 's/,/\t/g' "$i" | awk 'NR>1{print $0}' > ./matrix/"${name%_index.tmp_final.csv}"_INDEX.Z

done


for i in ./matrix/*_INDEX.Z; do
    name=$(basename "$i")
    # Create the FINAL_INDEX file
    echo "${name%_INDEX.Z}" > ./matrix/"${name}"_FINAL_INDEX
    
    # Process the INDEX file
    while read -r gene; do
        toprint=$(echo "$gene 0")
        if grep -wq "$gene" "$i"; then
            grep -w "$gene" "$i" >> ./matrix/"${name}"_FINAL_INDEX
        else
            echo "$toprint" >> ./matrix/"${name}"_FINAL_INDEX
        fi
    done < ./matrix/INDEX
    
    # Extract the last column
    awk '{print $NF}' ./matrix/"${name}"_FINAL_INDEX > ./matrix/"${name}"_FINALCOLUMN

done


paste ./matrix/matrix.tab ./matrix/*_FINALCOLUMN > ./matrix/final_matrix.tab


rm ./matrix/*_INDEX.Z ./matrix/*_final.csv ./matrix/*_FINAL_INDEX ./matrix/*INDEX.Z_FINALCOLUMN ./matrix/*_index.tmp ./matrix/INDEX ./matrix/matrix.tab

mv ./matrix/final_matrix.tab ./matrix/matrix.tab

tr '\n' ' ' < ./matrix/sample_names > ./matrix/names_heatmap

python heatmap.py ./matrix/matrix.tab ./matrix/names_heatmap



mv *png ./plots
