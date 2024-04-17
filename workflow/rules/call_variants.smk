'''
Go from BAM to calling variants (VCF)

AUTHOR: Henry Bloom

Li Lab, University of Chicago
This code significantly builds off of Chao Dai's A2I project for RNAseq data
'''

import os

# Call variants per sample
rule call_per_sample:
    message: "### Calling variants per sample ###"
    input: 
        bam = "results/{dataset}/{population}/BQSR/{runID}_recal.bam"
    output: 
        gvcf = "results/{dataset}/{population}/variants/{runID}.g.vcf.gz",
    params:
        ref = config["ref"]
    threads: 6
    resources: time=2100, mem_mb=30000, cpu=10
    shell:
        '''
        gatk HaplotypeCaller \
            -R {params.ref} \
            -I {input.bam} \
            -ERC GVCF \
            -O {output.gvcf}
        '''
        
rule call_per_sample_vcf:
    message: "### Calling variants per sample directly into vcf format ###"
    input: 
        bam = "results/{dataset}/{population}/BQSR/{runID}_recal.bam"
    output: 
        vcf = "results/{dataset}/{population}/vcfs/{runID}.vcf.gz",
    params:
        ref = config["ref"]
    threads: 6
    resources: time=2100, mem_mb=30000, cpu=10
    shell:
        '''
        gatk HaplotypeCaller \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf}
        '''

# HaplotypeCaller outputs an intermediate file, not a VCF directly (change from gvcf or whatever file it gets outputted as)
# - "Note that running HaplotypeCaller in GVCF mode generates intermediate GVCF files, which need to be combined using the GenomicsDBImport tool or CombineGVCFs tool before further analysis."

# Get the VCFs for a given individual
def getRunIDs(dataset, population):
    '''
    Get runIDs for a given individual
    '''
    variants_dir = os.path.join("results", dataset, population, "variants")
    runIDs = set(path[:path.find('.')] for path in os.listdir(variants_dir))
    return runIDs    

# make the sample map for GenomicsDBImport
rule make_sample_map:
    message: "### Making sample map for GenomicsDBImport ###"
    output:
        sample_map = "resources/{dataset}/{population}/sample_map.txt"
    resources: time=1000, mem_mb=1000, cpu=1
    run:
        dataset = ''.join({wildcards.dataset})
        population = ''.join({wildcards.population})
        runIDs = getRunIDs(dataset, population)
        with open(output.sample_map, "w") as f:
            for runID in runIDs:
                f.write(f"{runID}\tresults/{dataset}/{population}/variants/{runID}.g.vcf.gz\n")

### NEED TO FIGURE OUT A WAY TO CLEAR THIS DIRECTORY EACH TIME
# Consolidate GVCFs
rule consolidate:
    message: "### Consolidating GVCFs ###"
    input: 
        sample_map = "resources/{dataset}/{population}/sample_map.txt"
    output:
        datastore = directory("results/{dataset}/{population}/consolidated")
    params:
        intervals = config["interval_list"]
    threads: 2
    resources: time=1000, mem_mb=30000, cpu=4
    shell:
        '''
        gatk GenomicsDBImport \
            --intervals {params.intervals} \
            --sample-name-map {input.sample_map} \
            --genomicsdb-workspace-path {output.datastore}
        '''

# Joint-Call Cohort
rule joint_call:
    message: "### Joint calling cohort ###"
    input: 
        gvcfs = "results/{dataset}/{population}/consolidated"
    output:
        vcf = "results/{dataset}/{population}/joint_call/{population}.vcf.gz"
    params:
        ref = config["ref"]
    threads: 2
    resources: time=1000, mem_mb=30000, cpu=4
    shell:
        '''
        gatk GenotypeGVCFs \
            -R {params.ref} \
            -V gendb://{input.gvcfs} \
            -O {output.vcf}
        '''

# Snpeff

# Filter Variants by Variant (Quality Score) Recalibration
rule filter_variants:
    message: "### Filtering variants by Variant Recalibration ###"
    input: 
        vcf = "results/{dataset}/{population}/joint_call/{population}.vcf.gz"
    output:
        recal = "results/{dataset}/{population}/VSQR/{population}.recal",
        tranches = "results/{dataset}/{population}/VSQR/{population}.tranches"#,
        # plots = "results/{dataset}/{population}/VSQR/{population}.plots.R"
    params:
        ref = config["ref"],
        known_sites = config["known_sites"],
        hapmap = config["hapmap"],
        omni = config["omni"],
        _1000G = config["_1000G"]
    threads: 2
    resources: time=2100, mem_mb=30000, cpu=8
    shell:
        '''
        gatk VariantRecalibrator \
            -R {params.ref} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 {params.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {params._1000G} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.known_sites} \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
            -mode SNP \
            -O {output.recal} \
            --tranches-file {output.tranches} 
        '''
# Not sure what some of the bottom things are but basically it needs to train on tons of known data to give more nuanced filtering of variants

# Whole genomes and exomes take slightly different parameters, so make sure you adapt your commands accordingly!
# LOOK AT https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator

rule ApplyVQSR:
    message: "### Applying Variant Recalibration ###"
    input: 
        vcf = "results/{dataset}/{population}/joint_call/{population}.vcf.gz",
        recal = "results/{dataset}/{population}/VSQR/{population}.recal",
        tranches = "results/{dataset}/{population}/VSQR/{population}.tranches"
    output:
        filtered_vcf = "results/{dataset}/{population}/filteredVariants/{population}.filtered.vcf.gz"
    params:
        ref = config["ref"]
    threads: 1
    resources: time=2100, mem_mb=30000, cpu=8
    shell:
        '''
        gatk ApplyVQSR \
            -R {params.ref} \
            -V {input.vcf} \
            --recal-file {input.recal} \
            --create-output-variant-index true \
            --tranches-file {input.tranches} \
            -mode SNP \
            -O {output.filtered_vcf}
        '''
# could also add       --truth-sensitivity-filter-level NUMBER \ to the above command