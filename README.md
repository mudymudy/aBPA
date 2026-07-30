[![DOI](https://zenodo.org/badge/855715084.svg)](https://doi.org/10.5281/zenodo.19684212)

### Panchronos workflow
`panchronos` is a pipeline designed for microbial phylogenetic reconstruction based on pangenome alignment. It provides an end-to-end workflow, with particular emphasis on handling low-quality and fragmented data, as commonly encountered in ancient DNA studies.


## Table of Contents

- [Installation](#1-installation)
- [Perform test run](#2-perform-test-run)
- [Input](#3-input)
- [Configuration](#4-configuration)
- [Running the pipeline](#5-running-the-pipeline)
- [Output](#6-output)
- [Tools installed by panchronos](#7-tools-installed-by-panchronos)
- [Documentation](#8-documentation)
- [Troubleshooting](#9-troubleshooting)




# 1/ Installation
First you need to have `Nextflow` and `Conda/Mamba` installed. To install `Nextflow`, please follow the instructions: https://www.nextflow.io/docs/latest/install.html. For `Conda`, if you don't already have any version installed, we recommend using Miniforge: https://github.com/conda-forge/miniforge

NOTE: `panchronos` has been tested with nextflow versions >= 24 and < 26. If you have nextflow version >=26, you may need to execute `export NXF_SYNTAX_PARSER=v1` before running the pipeline.

To verify that Nextflow is properly installed and available in your environment, please run: 
```
nextflow info
```

### Then git clone this repository:

```
git clone https://github.com/brgonzlez/panchronos/
```

After cloning the repository, you should see the directories bin/ config/ envs/, along with the files `main.nf` and `nextflow.config`.

# 2/ Perform test run

To perform a test run, please go to the main directory where `panchronos` was cloned (the folder containing the `main.nf` file) and execute the following command:

**IMPORTANT: When `Nextflow` is installing a `Conda` environment, do not interrupt or cancel the process, as it can result in a corrupted environment that must be deleted and rebuilt.**

```
dir=$(pwd)
mkdir -p "$dir"/test/results_test/
nextflow run main.nf --test --config "$dir"/test/config_test.tab --output "$dir"/test/results_test --tax_id 10255 --outgroup_tax_id 28871 --genomes 5
```

This test run will download five modern genomes of variola virus (`--genomes 5`, `--tax_id 10255`) and one genome of taterapox virus (`--outgroup_tax_id 28871`) as an outgroup. It will then use one of the variola virus genomes as a template to generate short-read data with damage patterns by `gargammel`. The `Conda` environments created during this test run will be reused for future executions. Running the test is optional—these environments will also be created automatically the first time you run the pipeline with real data.

# 3/ Input

The `panchronos` workflow assumes that you already know which bacteria or virus is present in your metagenomic dataset. If you are unsure, you should perform metagenomic screening beforehand to identify candidate organisms.

To run `panchronos`, the workflow requires the following four components:
(1) `.fastq` files
(2) `config.tab` file
(3) Taxonomic IDs
(4) Paths to your data and output directories

# (1) FASTQ files
- The workflow accepts compressed (`.gz`) or uncompressed `.fastq` or `.fq` files.
- All input FASTQ files should be placed in the same directory if you are analysing multiple datasets.
- `panchronos` does not support paired-end data. If you have paired-end reads, you must collapse/merge them beforehand.
- Multiple single-end/collapsed libraries from the same individual can be included. Please see `config.tab` below for grouping instructions.

# (2) `config.tab` file
`config.tab` is a tab-separated text file with four fields:
| Sample filename | Trimming value | Sample ID | Status |
|----------|----------|----------|----------|
| Sample_A.fastq | 2 | Individual_1  | A |
| Sample_B.fastq | 2 | Individual_1  | A |
| Sample_C.fastq | 5 | Individual_1  | A |
| Sample_D.fastq | 0 | Individual_2  | M |

The workflow uses this file to:
- Apply sample-specific trimming*
- Merge aligned data for samples belonging to the same individual (i.e., those sharing the same Sample ID) after alignment.
- Apply different data processing strategies for ancient (`A`) and modern (`M`) datasets

**Note:**
If the `--rescale` option is set to 0, panchronos trims _n_ bases from both ends of each read after alignment, based on the specified trimming value, instead of performing soft-clipping.
The config.tab file must include all four fields described above, even when the `--rescale` parameter is set to 1.


You can see an example file in the `config/` directory of the repository.

# (3) Required taxonomic IDs
You need to provide two NCBI taxonomic IDs:


- Target species taxonomic ID — used for pangenome construction. The workflow will automatically download the necessary genomic data based on this ID. You can control how many genomes are downloaded for pangenome construction using the `--genomes` parameter.

- Outgroup taxonomic ID — used only for rooting the phylogenetic tree.


If you already have your own curated dataset (FASTA + GenBank files), you can provide its path using `--trusted_data`. If this flag is declared, `panchronos` will not download genomic data and will ignore the value for `--tax_id`. However, you still need to provide the outgroup taxonomic ID.


****How to obtain a NCBI taxonomic ID****


Go to the NCBI Taxonomy database:
https://www.ncbi.nlm.nih.gov/taxonomy
Enter the name of your organism (e.g., _Escherichia coli_) in the search bar.
Click on the relevant entry in the results list.
The taxonomic ID (TaxID) is displayed on the organism’s page (e.g., _Escherichia coli_ has Taxonomic ID 562).


Important requirements:
- FASTA and GenBank filenames must match exactly (aside from extensions).
- Required extensions: `.fasta` and `.gb`

# (4) Data and output directories
Before running the workflow, prepare two directories:
- A folder containing only your `FASTQ/FASTQ.gz` files. **Only store the data that will be included in the analysis.**
- A separate folder where pipeline outputs will be written.

Provide their paths with the parameters:
- `--data` (input FASTQ directory)
- `--output` (output directory)

See the Configuration section below.

# 4/ Configuration
The workflow's configuration is controlled through the `nextflow.config` file. This file contains all available parameters, allowing you to fine-tune the pipeline according to your analysis needs. You can adjust the default values directly to match your setup and computational resources.

Many parameters are associated with specific workflow modules, and CPU/thread usage is particularly important to configure. In addition to threads used by individual tools, the pipeline may spawn parallel jobs, multiplying the total number of threads in use.
For example:
- `--alignment_threads 10`
- `--alignment_parallel 5`
- 5 samples in the `config.tab`

This combination means the alignment process will use 10 × 5 = 50 threads simultaneously.

**Important NOTE:**
**Do not assign more than 3 threads to `--get_data_parallel`.**

NCBI may reject download requests if the number of parallel queries exceeds 3.

# 5/ Running the pipeline
Once you have configured your settings in the `nextflow.config` file (including all required paths), you can execute the workflow with: 

>`nextflow run main.nf -resume` 

It is recommended to include the `-resume` flag for every run, as this enables Nextflow’s caching system and prevents unnecessary recomputation.

**Overriding parameters from the command line**

You can override any parameter directly from the terminal by prefixing it with `--`.
Values provided this way take precedence over those defined in `nextflow.config`:

>`nextflow run main.nf --data /new/path/ --panaroo_alignment_type core -resume`

In this example, the values `/new/path/` and `core` will be used instead of the corresponding entries in `nextflow.config`.

When the workflow starts, `panchronos` prints in the terminal a summary of all parameters used for the current run.
This allows you to quickly verify that your settings are being applied correctly.

# 6/ Output
As the workflow runs, the `--output` directory will begin to populate with results. Files are written as soon as they are generated and automatically organized into structured subdirectories, making it easy to track each stage of the analysis.

| Folder | Description | Format |
| :--- | :--- | :---: | 
|**ALIGNMENT**| Aligned reads against the pangenome reference sequences | `.bam`, `.fastq` |
|**DOWNLOADED**| Input genomes used for pangenome construction. If `--trusted_data` is provided, this folder stores the user-supplied curated datasets | `.fasta`, `.gb` |
|**GENE_DATABASE**| Guide file for `Prokka`, containing clustered genes extracted from `GenBank` files prior to annotation | `.fasta` |
|**GENE_MSA**| Multiple sequence alignments for each individual gene | `.fasta` |
|**GENOTYPING**| Variant-calling results, including consensus sequences per sample per gene, and associated VCF files | `.fasta`, `.vcf*` |
|**MAPDAMAGE**| Output from DNA damage pattern assessment using `mapDamage` | Multiple files |
|**MATRIX**| Final gene presence/absence matrix | `.tab` |
|**PANGENOME**| Outputs from pangenome construction and extension. Includes the original Panaroo reference genome, extended/unextended versions, and a BLAST database | `.fasta`, BLAST database |
|**PLOTS**| Visualizations including presence/absence heatmaps and coverage vs. completeness plots | `.png` |
|**STATS**| Summary statistics, including gene-level coverage and completeness normalization per sample | `.tab`, `.txt` |
|**TREE**| Results from IQ-TREE, including the concatenated MSAs used for reconstruction | `.fasta`, `.treefile`, other outputs from IQ-TREE |

# 7/ Tools installed by `panchronos`
`panchronos` has only been tested with the specific tool versions listed below, all of which are installed automatically when the workflow is executed. Compatibility with newer or different releases is not guaranteed.

- Entrez Direct v24.0
- Biopython v1.85
- fastANI v1.34
- CD-HIT v4.8.1
- Prokka v1.14.6
- Panaroo v1.5.2
- bedtools v2.31.1
- seqtk v1.5
- BWA v0.7.18
- bamUtil v1.0.15
- SAMtools v1.19
- mapDamage v2.2.3
- art_illumina v2.5.8
- bcftools v1.21
- MAFFT v7.525
- IQ-TREE v3.0.1
- picard v3.4.0
- BLAST v2.16.0+
- gargammel v1.1.4


# 8/ Documentation
If you are unsure about the available workflow settings or parameters, you can display all options directly from the terminal using:

```
nextflow run main.nf --help
```

For a more comprehensive documentation, please click on this [link](https://shigekinakagomelab.com/s/panchronos_v10-beta15_Documentation.pdf)




# 9/ Troubleshooting

# (1) Errors related to NCBI servers

Sometimes you may encounter errors such as:


>`HTTP/1.1 502 Bad Gateway`


This is not related to a malfunction of `panchronos`. It indicates a temporary issue with the NCBI servers (e.g., proxy or upstream service problems). In such cases, there is usually nothing to fix on your side. The issue is transient and typically resolves on its own. If the error persists consistently, it may be worth checking your query or network connection, but most occurrences are short-lived server-side issues.


If you get an error that is not listed in this section, please raise a new issue desccribing it.
