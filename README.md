# Predicting Pathogenic Gene Splicing Mutations from DNA Sequencing Data

This repository contains a pipeline for identifying pathogenic gene splicing mutations from DNA sequencing data. 

## Overview
Rare genetic diseases are caused by mutations in a single or a small number of genes and are commonly referred to as Mendelian diseases. Recent advances in genome sequencing have increased the discovery rate of genetic variants that cause Mendelian diseases, yet 50–75% of patients remain undiagnosed for many rare disorders. While recent studies have suggested that disruption of RNA splicing can cause rare diseases, few pipelines have been designed to broadly identify pathogenic variants that impact RNA splicing from clinical sequencing data. Our research aims to leverage modern variant calling pipelines with new splice site prediction models for this task. Specifically, we integrated Pangolin, a deep learning model that accurately identifies loss of function splicing mutations, with the germline variant discovery pipeline from the Genome Analysis Toolkit. Pipeline performance was evaluated using Illumina’s “hap.py” framework to benchmark variant calls against gold standard datasets, enabling optimization of the parameters throughout the pipeline. Our tool holds promise to aid in diagnosing patients with rare genetic disorders and provide insights into previously characterized human DNA sequencing data. Identified genetic variants can serve as targets for patient-specific gene therapies and enable the treatment of otherwise untreatable disorders. In the future, this pipeline can be enhanced to predict more complex pathogenic structural variants and could be tested with new state-of-the-art variant calling pipelines such as DeepVariant.

## Installation and Setup
To install the pipeline, clone the repository and create a conda environment with the necessary dependencies. You can use the provided `create_env.sh` file to create the environment given that you have [`mamba`](https://github.com/mamba-org/mamba) installed. Note that this environment is sufficient for running all steps of the pipeline except Pangolin. For that, you will need to create a separate environment by following the steps [here](https://github.com/tkzeng/Pangolin?tab=readme-ov-file#installation).

To run the pipeline, you first need to create a configuration file called `config.yaml` in the `config` directory. This file should contain paths to the following:
- `ref`: Path to the reference genome fasta file.
- `known_sites`: Path to the known sites vcf file, for example from dbSNP.
- `tmp`: Path to a temporary directory.
- `interval_list`: Path to the interval list file which simply lists chromosomes of interest when variant calling, each on a new line.
- `pangolin_db`: Path to the database used when calling [Pangolin](https://github.com/tkzeng/Pangolin?tab=readme-ov-file#usage-command-line).
- `1000G`, `omni`, `hapmap`: Paths to the corresponding files used in the variant quality score recalibration step. These files can be downloaded from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

## Usage
As the pipeline is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/), you can simply run Snakemake while requesting a given file you would like outputted, and the pipeline will automatically run all necessary steps to generate that file. For example, to generate the a variant call for a single sample from dataset `data1`, population `pop1`, and run `run1`, you can run the following command:

```bash
snakemake --cores 1 results/data1/pop1/variants/run1.g.vcf.gz
```

Snakemake will create a DAG representing the pipeline and run all necessary steps to generate the requested file. Note that you can specify a number of parameters to Snakemake. I recommend creating a Snakemake profile to manage these parameters. You can find more information on how to do this [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

## In-Depth Description of Steps
The pipeline is implemented in Snakemake and consists of the following steps:

### 1. Remove Adapters from Paired End Reads
First, any adapters in the input sequences will be trimmed using [`fastp`](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters).

**Rule:** `RemoveAdapters_paired`  

**Input:** Two paired end fastq files compressed with gzip. These must be stored in `resources/{dataset}/{population}/fastq-dl/{runID}_1.fastq.gz` and `resources/{dataset}/{population}/fastq-dl/{runID}_2.fastq.gz` respectively. 

**Output:** Two paired end fastq files compressed with gzip. These will be stored in `resources/{dataset}/{population}/trimmed-fastq/{runID}_1_trimmed.fastq.gz` and `resources/{dataset}/{population}/trimmed-fastq/{runID}_2_trimmed.fastq.gz` respectively. This rule also outputs a file containing information about the samples, stored in `resources/{dataset}/{population}/trimmed-fastq/{runID}.fastp_report.html`.

### 2. Align Reads to Reference Genome
Next, the trimmed reads will be aligned to the reference genome using [`BWA mem`](https://bio-bwa.sourceforge.net/index.shtml).

**Rule:** `BWA_paired`

**Input:** Two paired end fastq files compressed with gzip. These must be stored in `resources/{dataset}/{population}/trimmed-fastq/{runID}_1_trimmed.fastq.gz` and `resources/{dataset}/{population}/trimmed-fastq/{runID}_2_trimmed.fastq.gz` respectively. 

**Output:** A single bam file containing the aligned reads. This will be stored in `results/{dataset}/{population}/mapped/{runID}.bam`.

### 3. Filter Unmapped Reads
Unmapped reads will be filtered out from the bam file using [`samtools view`](http://www.htslib.org/doc/samtools-view.html).

**Rule:** `FilterAlignment`

**Input:** A single bam file containing the aligned reads. This must be stored in `results/{dataset}/{population}/mapped/{runID}.bam`.

**Output:** A single bam file containing only the mapped reads. This will be stored in `results/{dataset}/{population}/filteredAlignment/{runID}_filtered.bam`.

### 4. Sort Alignments
The filtered bam file will be sorted using [`samtools sort`](https://www.htslib.org/doc/samtools-sort.html).

**Rule:** `SortAlignments`

**Input:** A single bam file containing the aligned reads. This must be stored in `results/{dataset}/{population}/filteredAlignment/{runID}_filtered.bam`.

**Output:** A single sorted bam file. This will be stored in `results/{dataset}/{population}/mapped/{runID}_sorted.bam`.

### 5. Mark Duplicates
Duplicate reads will be marked using [`gatk MarkDuplicates`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard).

**Rule:** `MarkDups`

**Input:** A single sorted bam file. This must be stored in `results/{dataset}/{population}/mapped/{runID}_sorted.bam`.

**Output:** A single bam file with duplicates marked. This will be stored in `results/{dataset}/{population}/MarkDups/{runID}.bam`. This rule also outputs a file containing information about the duplicates, stored in `results/{dataset}/{population}/MarkDups/{runID}_metrics.txt`.

### 6. Base Quality Score Recalibration
Next, the base quality scores will be recalibrated using a two step process. First, a recalibration table is created using [`gatk BaseRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator).

**Rule:** `BaseRecalibrator`

**Input:** A single bam file with duplicates marked. This must be stored in `results/{dataset}/{population}/MarkDups/{runID}.bam`.

**Output:** A single recalibration table. This will be stored in `results/{dataset}/{population}/BQSR/{runID}_covariates.tab`.

### 7. Apply Base Quality Score Recalibration
The recalibration table is then used to recalibrate the base quality scores of the reads using [`gatk ApplyBQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR).

**Rule:** `ApplyBQSR`

**Input:** A single bam file with duplicates marked. This must be stored in `results/{dataset}/{population}/MarkDups/{runID}.bam`. Also, a recalibration table must be present in `results/{dataset}/{population}/BQSR/{runID}_covariates.tab`.

**Output:** A single bam file with recalibrated base quality scores. This will be stored in `results/{dataset}/{population}/BQSR/{runID}_recal.bam`.

### 8. Call Variants
Our data has now been fully processed, and we can call variants per sample using [`gatk HaplotypeCaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) in GVCF mode. Note that an alternate rule, `call_per_sample_vcf`, is also available for calling variants in VCF mode.

**Rule:** `call_per_sample`

**Input:** A single bam file with recalibrated base quality scores. This must be stored in `results/{dataset}/{population}/BQSR/{runID}_recal.bam`.

**Output:** A single g.vcf file containing the variants called for the sample. This will be stored in `results/{dataset}/{population}/variants/{runID}.g.vcf.gz`.

### 9. Making a Sample Map
Next, we will create a sample map that contains the paths to the g.vcf files for each sample. This will be used in the next step to combine the g.vcf files into a single file.

**Rule:** `make_sample_map`

**Input:** There is no explicit input for this step. The corresponding directory of g.vcf files will be scanned to create this map.

**Output:** A single file containing the paths to the g.vcf files for each sample. This will be stored in `resources/{dataset}/{population}/sample_map.txt`.

### 10. Combine GVCFs
The g.vcf files for each sample are combined into a single file using [`gatk GenomicsDBImport`](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport).

**Rule:** `consolidate`

**Input:** A sample map containing the paths to the g.vcf files for each sample. This must be stored in `resources/{dataset}/{population}/sample_map.txt`.

**Output:** A single directory containing the information for the newly created database. This will be stored in `results/{dataset}/{population}/consolidated`.

### 11. Joint Genotyping
The variants are jointly genotyped across all samples using [`gatk GenotypeGVCFs`](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs).

**Rule:** `joint_call`

**Input:** A directory containing the information for the newly created database. This must be stored in `results/{dataset}/{population}/consolidated`.

**Output:** A single vcf file containing the joint genotyped variants. This will be stored in `results/{dataset}/{population}/joint_call/{population}.vcf.gz`.

### 12. Filter Variants
The joint genotyped variant quality scores are recalibrated in a two step process similar to the base quality score recalibration. First, a recalibration table is created using [`gatk VariantRecalibrator`](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator).

**Rule:** `filter_variants`

**Input:** A single vcf file containing the joint genotyped variants. This must be stored in `results/{dataset}/{population}/joint_call/{population}.vcf.gz`.

**Output:** A single recalibration table and a tranches file that shows various metrics of the recalibration callset for slices of the data. These are stored in `results/{dataset}/{population}/VSQR/{population}.recal` and `results/{dataset}/{population}/VSQR/{population}.tranches` respectively.

### 13. Apply Variant Quality Score Recalibration
The recalibration table is then used to recalibrate the variant quality scores of the joint genotyped variants using [`gatk ApplyVQSR`](https://gatk.broadinstitute.org/hc/en-us/articles/360036485072-ApplyVQSR).

**Rule:** `ApplyVQSR`

**Input:** A single vcf file containing the joint genotyped variants. This must be stored in `results/{dataset}/{population}/joint_call/{population}.vcf.gz`. Also, a recalibration table and a tranches file must be present in `results/{dataset}/{population}/VSQR/{population}.recal` and `results/{dataset}/{population}/VSQR/{population}.tranches` respectively.

**Output:** A single vcf file with recalibrated variant quality scores. This will be stored in `results/{dataset}/{population}/filteredVariants/{population}.filtered.vcf.gz`.

### 14. Annotate Variants
The recalibrated vcf file will be annotated using [`SnpEff`](https://pcingola.github.io/SnpEff/#snpeff).

**Rule:** `annotate`

**Input:** A single vcf file with recalibrated variant quality scores. This must be stored in `results/{dataset}/{population}/filteredVariants/{population}.filtered.vcf.gz`.

**Output:** A single vcf file with annotated variants. This will be stored in `results/{dataset}/{population}/annotated/{population}.ann.vcf.gz`. Also, a series of statistics will be outputted to `results/{dataset}/{population}/annotated/{population}.ann.html`.

### 15. Filter Annotated Variants
The annotated vcf file will be filtered to only include variants that are predicted to have a high impact using [`SnpSift](https://pcingola.github.io/SnpEff/#snpsift).

**Rule:** `sift`

**Input:** A single vcf file with annotated variants. This must be stored in `results/{dataset}/{population}/annotated/{population}.ann.vcf.gz`.

**Output:** A single vcf file with filtered annotated variants. This will be stored in `results/{dataset}/{population}/annotated/{population}.filt.vcf`.

### 16. Predict Splicing Mutations
The filtered annotated vcf file will be used to predict pathogenic gene splicing mutations using [`Pangolin`](https://github.com/tkzeng/Pangolin?tab=readme-ov-file). Note: due to conflicts in software versioning, it may be necessary to create a separate conda environment for running Pangolin.

**Rule:** `pangolin`

**Input:** A single vcf file with filtered annotated variants. This must be stored in `results/{dataset}/{population}/annotated/{population}.filt.vcf`.

**Output:** A single file containing the predicted pathogenic gene splicing mutations. This will be stored in `results/{dataset}/{population}/pangolin/{population}.pangolin.vcf.gz`.