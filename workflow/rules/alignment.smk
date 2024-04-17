'''
Go from FASTQ to alignment (BAM)

AUTHOR: Henry Bloom

Li Lab, University of Chicago
This code significantly builds off of Chao Dai's A2I project for RNAseq data
'''

rule makefai:
    message: '''### Indexing reference genome –– get the .fai file ###'''
    input:
        config['ref']
    output:
        config['ref'] + ".fai"
    resources:
        mem_mb=36000, cpu=6, partition="broadwl" # this took 65 minutes & a little under 5GB of memory (~1 SU)
    shell:
        '''
        samtools faidx {input}
        '''

rule index:
    message: '''### Indexing reference genome ###'''
    input:
        config['ref']
    output:
        config['ref'] + ".amb",
        config['ref'] + ".ann",
        config['ref'] + ".bwt",
        config['ref'] + ".pac",
        config['ref'] + ".sa" # indexed using bwa version 0.7.17-r1188
    resources:
        mem_mb=36000, cpu=6, partition="broadwl" # this took 65 minutes & a little under 5GB of memory (~1 SU)
    shell:
        '''
        bwa index {input} -p {input} 
        '''
        # the path to the output is the same as the input (the genome.fa file)


# Remove adapters for paired end reads
rule RemoveAdapters_paired:
    message: '''### Trim adapters ###'''
    input:
        f1 = "resources/{dataset}/{population}/fastq-dl/{runID}_1.fastq.gz",
        f2 = "resources/{dataset}/{population}/fastq-dl/{runID}_2.fastq.gz",
    output:
        f1 = "resources/{dataset}/{population}/trimmed-fastq/{runID}_1_trimmed.fastq.gz",
        f2 = "resources/{dataset}/{population}/trimmed-fastq/{runID}_2_trimmed.fastq.gz",
        report_html = "resources/{dataset}/{population}/trimmed-fastq/{runID}.fastp_report.html",
        report_json = "resources/{dataset}/{population}/trimmed-fastq/{runID}.fastp_report.json"
    threads: 1
    resources: time=600, mem_mb=20000, cpu=1
    shell:
        '''
        fastp --in1 {input.f1} --in2 {input.f2} --out1 {output.f1} --out2 {output.f2}\
            --thread {threads} \
            --dont_overwrite \
            --overrepresentation_analysis \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --disable_quality_filtering \
            --length_required 15 \
            --json {output.report_json} \
            --html {output.report_html}
        '''
    

# Remove adapters for single end reads
# rule RemoveAdapters_single:
#     message: '''### Trim adapters ###'''
#     input:
#         f1 = "resources/{dataset}/{population}/fastq-dl/{runID}.fastq.gz",
#     output:
#         f1 = "resources/{dataset}/{population}/trimmed-fastq/{runID}_trimmed.fastq.gz",
#         report_html = "resources/{dataset}/{population}/trimmed-fastq/{runID}.fastp_report.html",
#         report_json = "resources/{dataset}/{population}/trimmed-fastq/{runID}.fastp_report.json"
#     threads: 1
#     resources: time=200, mem_mb=20000, cpu=1
#     shell:
#         '''
#         fastp --in1 {input.f1} --out1 {output.f1} \
#             --thread {threads} \
#             --dont_overwrite \
#             --overrepresentation_analysis \
#             --detect_adapter_for_pe \
#             --trim_poly_g \
#             --disable_quality_filtering \
#             --length_required 15 \
#             --json {output.report_json} \
#             --html {output.report_html}
#         '''

# BWA alignment for Paired-End Reads
rule BWA_paired:
    message: '''### Alignment using BWA on paired-end reads ###'''
    input:
        ref = config['ref'],
        f1 = "resources/{dataset}/{population}/trimmed-fastq/{runID}_1_trimmed.fastq.gz",
        f2 = "resources/{dataset}/{population}/trimmed-fastq/{runID}_2_trimmed.fastq.gz"
    output: 
        bam = "results/{dataset}/{population}/mapped/{runID}.bam"
    params:
        readgroup = "\'@RG\\tID:{runID}\\tSM:{runID}\\tPL:ILLUMINA\\tLB:{runID}\\tPU:{runID}\'"
    threads: 
        6
    resources: 
        time=2100, mem_mb=40000, cpu=10, partition="broadwl"
    shell: 
        '''
        bwa mem -R {params.readgroup} {input.ref} {input.f1} {input.f2} | samtools view -Sb - > {output}
        '''

# BWA alignment for Single-End Reads
rule BWA_single:
    message: '''### Alignment using BWA on single-end reads ###'''
    input:
        ref = config['ref'],
        fa = "resources/{dataset}/{population}/trimmed-fastq/{runID}_trimmed.fastq.gz",
    output: 
        bam = "results/{dataset}/{population}/mapped/{runID}.bam"
    params:
        readgroup = "\'@RG\\tID:{runID}\\tSM:{runID}\\tPL:ILLUMINA\\tLB:{runID}\\tPU:{runID}\'"
    threads: 
        6
    resources: 
        time=2100, mem_mb=36000, cpu=6, partition="broadwl"
    shell: 
        '''
        bwa mem -R {params.readgroup} {input.ref} {input.fa} | samtools view -Sb - > {output}
        '''

# mem is recommended for a certain length of reads (longer) –– for older (shorter) reads, use aln ––> read carefully about the differences
    # if < ~75, use aln
# less the fastq file and take a count to check if it's the right length
# make sure to understand the output format
# using -R to add read group information

# filtering out unmapped reads and enforcing MAPQ >= 20
rule FilterAlignment:
    message: 
        '''### Filter BAM on flag, mapq, chrom ### '''
    input: 
        "results/{dataset}/{population}/mapped/{runID}.bam"
    output: 
        "results/{dataset}/{population}/filteredAlignment/{runID}_filtered.bam"
    threads: 
        2
    resources: 
        time=500, mem_mb=30000, cpu=4
    shell:
        '''
        samtools view -h -F 4 -q 20 -b {input} > {output}
        '''
# "-F 4" excludes unmapped reads (flag 4)
# "-q 20" will enforce MAPQ >= 20 (can change number value)
# https://davetang.org/wiki/tiki-index.php?page=SAMTools#Filtering_out_unmapped_reads_in_BAM_files

# Sort BAM
rule SortAlignments:
    message: ''' ### Sort BAM using samtools ###'''
    input:
        "results/{dataset}/{population}/filteredAlignment/{runID}_filtered.bam"
    output:
        "results/{dataset}/{population}/mapped/{runID}_sorted.bam"
    threads: 
        2
    resources: 
        time=500, mem_mb=30000, cpu=4
    shell:
        """
        samtools sort -o {output} {input}
        samtools index -@ {threads} {output}   
        """

# Mark duplicates
rule MarkDups:
    message: '''### Mark Duplicates ###'''
    input:
        bam = "results/{dataset}/{population}/mapped/{runID}_sorted.bam"
    output:
        bam = "results/{dataset}/{population}/MarkDups/{runID}.bam",
        metrics = "results/{dataset}/{population}/MarkDups/{runID}_metrics.txt"
    params:
        tmp = "/scratch/midway2/hmbloom/TMP" # added because gatk's tmp directory error
    threads: 
        2
    resources: 
        time=2000, mem_mb=30000, cpu=4
    shell:
        '''
        gatk MarkDuplicates  \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --TMP_DIR {params.tmp} \
            --REMOVE_SEQUENCING_DUPLICATES
        '''

# BASE QUALITY SCORE RECALIBRATION (two steps):
# step 1: compute scores
rule BaseRecalibrator:
    message: '''### Compute covariate matrix of Base Recalibration ###'''
    input:
        bam = "results/{dataset}/{population}/MarkDups/{runID}.bam",
        ref = config['ref'],
        known_sites = config['known_sites']
    output:
        recal_file = "results/{dataset}/{population}/BQSR/{runID}_covariates.tab"
    params:
        tmp = "/scratch/midway2/hmbloom/TMP" # added because gatk's tmp directory errorthreads: 1
    resources: time=2000, mem_mb=30000, cpu=4
    shell:
        '''
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.known_sites} \
            --tmp-dir {params.tmp} \
            -O {output.recal_file}
        '''

"""
-R, --reference              reference file
-I, --input                  input BAM/SAM/CRAM file(s)
-known-sites                 known variants file(s)
-O, --output                 output recalibration table file
--tmp-dir                    tmp directory to use
"""

# step 2: apply to the old BAM to make a new BAM
rule ApplyBQSR:
    message: '''### Apply Base Recalibration, output calibrated BAM ###'''
    input:
        bam = "results/{dataset}/{population}/MarkDups/{runID}.bam",
        ref = config['ref'],
        recal_file = rules.BaseRecalibrator.output.recal_file
    output:
        bam = "results/{dataset}/{population}/BQSR/{runID}_recal.bam"
    params:
        tmp = "/scratch/midway2/hmbloom/TMP" # added because gatk's tmp directory errorthreads: 1
    resources: time=2000, mem_mb=30000, cpu=4
    shell:
        '''
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.recal_file} \
            --tmp-dir {params.tmp} \
            -O {output.bam}
        '''

# NOTE: WHEN DONE WITH THIS STUFF, FIND 3 SAMPLES THAT ARE NOT RELATED AND RUN THEM THROUGH THE PIPELINE 
# -- (restrict to chromosome 22 to make faster)
# use   samtools view    and specify a certain region (chromosome) -- select region
# -- select all alignments from a region and output to a separate file --> working with just chr22 speeds up the process
# BWA should take a while 