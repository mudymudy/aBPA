# aBPA: ancient Bacterial Pangenome Analysis

ancient Bacterial Pangenome Analysis is a `nextflow` package.

# 1/ Installation


First you need to have `nextflow` and `conda/mamba` installed. Please visit: https://www.nextflow.io/docs/latest/install.html and follow the instructions. In the case of conda, if you don't have any conda version installed I would recomend miniforge (https://github.com/conda-forge/miniforge).


After installing `nextflow`, if you wish to be able to call the executable without explicitly define the full PATH while not having to use sudo:


### Open your bashrc file:


>`$ nano .bashrc`


### Then add this line at the end: 


>`$ export PATH=$PATH:/home/users/myuser/Softwares/`

### Finally refresh the environment to activate the changes


>`$ source ~/.bashrc`


### Note that the PATH will be unique in your system. In my case I installed nextflow in the folder /Softwares/


### Then git clone this repository:



>`$ git clone https://github.com/mudymudy/aBPA/`


After downloading the repository you should see the folders `bin/` `config/` `envs/` and the files `aBPA.nf` and `nextflow.config`.
If you can't be bothered to type nextflow run aBPA.nf every time, you can do (assuming you already exported nextflow $PATH into your .bashrc file):

>`$ echo "alias aBPA='nextflow run aBPA.nf'" >> ~/.bashrc`


>`$ source ~/.bashrc`



Then you can just type aBPA in the same directory where aBPA.nf is located.


# 2/ First steps


First things first. aBPA assumes you already know what bacteria is in your data (if you have metagenomic data). If you don't know what is in your data you need to do metagenomic profiling first. There a a number of tools that can perform this such as eager, aMeta or mapache.

Once you have a bacteria in mind, then you need to know the taxonomic ID and the taxonomic ID of another bacteria that you want to use as outgroup for phylogenetic reconstruction. 


### Making the `config.tab` file

The way the pipeline reads your data is through the config.tab file, which has to be located in the `config/` folder. The structure is as follows:




# 3/ Documentation



globalMeanCoverage.txt contains:

sampleID

sampleCoverage = Total number of bases of the reference genome genes where there is at least 1 read covering.

refCount = Total number of bases of the reference genome, including genes that are not being covered by reads.

globalMean = Mean depth of coverage of every gene where there is at least 1 read covering.

geneNormalizedSummary.txt contains:

sampleID
gene = Gene name

normalizedGeneSimple = (geneMeanDepth / globalMean)
-	geneMeanDepth -> Mean Depth of coverage of a particular gene.

normalizedGeneScaled = (geneMeanDepth / globalMean) * geneLength[gene]
-	geneLength[gene] -> Length of a particular gene.

normalizedGenomeSimple = (geneMeanDepth / globalMean) * (geneLength[gene] / sampleCoverage)

normalizedGenomeScaled = (geneMeanDepth / globalMean) * (geneLength[gene] / refCount)




# 4/ TODO



Add paired end reads options.



Add UnifiedGenotyper option for --genotype and classify the current as --genotype bcftools. Improve the way bcftools is behaving. DONE, BUT STILL NEED TO MAKE BCFTOOLS A BIT BETTER



Add the option to run parsnp instead of pMauve with --aligner



Make alignment parameters with bwa aln variables that can be tuned. DONE



Add a antimicrobial resistance pipe and link it with genes.



Add pre-procesing step such as alignment against human reference , adapter removal and deduplication.



Add heteroplasmy process before getting genotypes and update heterozygosis values into the updatedNormalized table. WORKING ON IT



Add an option to filter reads based on length. DONE




Maybe add MapDamage to get more metrics after alignment.




Make mapping quality threshold a variable in alignment process. DONE




Make documentation and improve --help message.




Add a small test run with tiny dataset.




Fix python heatmap labels related to number of samples. Make labels a bit bigger.




Improve the diagram by making letters bigger.




Evaluate the option to work with pangenome graph as well.






