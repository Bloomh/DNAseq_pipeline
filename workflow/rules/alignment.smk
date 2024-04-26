'''
Go from FASTQ to alignment (BAM)

AUTHOR: Henry Bloom

Li Lab, University of Chicago
This code significantly builds off of Chao Dai's A2I project for RNAseq data
'''

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
        time=2100, mem_mb=40000, cpu=10
    shell: 
        '''
        bwa mem -R {params.readgroup} {input.ref} {input.f1} {input.f2} | samtools view -Sb - > {output}
        '''

# filtering out unmapped reads
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
        tmp = config['tmp']
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
        tmp = config['tmp']
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