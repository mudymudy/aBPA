# aBPA: ancient Bacterial Pangenome Analysis

ancient Bacterial Pangenome Analysis is a `nextflow` package.

# 1/ Installation


First you need to have `nextflow` and `conda/mamba` installed. Please visit: https://www.nextflow.io/docs/latest/install.html and follow the instructions. In the case of conda, if you don't have any conda version installed I would recomend miniforge (https://github.com/conda-forge/miniforge).


After installing `nextflow`, if you wish to be able to call the executable without explicitly define the full PATH while not having to use sudo:


### Open your bashrc file:


>`nano .bashrc`


### Then add this line at the end: 


>`export PATH=$PATH:/home/users/myuser/Softwares/`

### Finally refresh the environment to activate the changes


>`source ~/.bashrc`


### Note that the PATH will be unique in your system. In my case I installed nextflow in the folder /Softwares/


### Then git clone this repository:



>`git clone https://github.com/mudymudy/aBPA/`


After downloading the repository you should see the folders `bin/` `config/` `envs/` and the files `aBPA.nf` and `nextflow.config`.



# 2/ First steps


First things first. aBPA assumes you already know what bacteria is in your data (if you have metagenomic data). If you don't know what is in your data you need to do metagenomic profiling first. There a a number of tools that can perform this such as eager, aMeta or mapache.

Once you have a bacteria in mind, then you need to know the taxonomic ID and the taxonomic ID of another bacteria that you want to use as outgroup for phylogenetic reconstruction. 


### config.tab file

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



