'''
Annotation and downstream analysis of variants

AUTHOR: Henry Bloom
'''

rule annotate:
    message: '''### Annotating variants with SnpEff ###'''
    input: 
        vcf = "results/{dataset}/{population}/joint_call/{population}.vcf.gz"
    output:
        vcf = "results/{dataset}/{population}/annotated/{population}.ann.vcf.gz",
        stats = "results/{dataset}/{population}/annotated/{population}.ann.html"
    params:
        snpeff = "/home/hmbloom/snpEff/snpEff.jar",
        config = "/home/hmbloom/snpEff/snpEff.config"
    shell: 
        '''
        java -Xmx8g -jar {params.snpeff} -c {params.config} -v -stats {output.stats} GRCh38.99 {input.vcf} | bgzip > {output.vcf}
        '''

rule sift:
    message: '''### Filtering annotated variants with SnpSift ###'''
    input: 
        vcf = "results/{dataset}/{population}/annotated/{population}.ann.vcf.gz"
    output:
        vcf = "results/{dataset}/{population}/annotated/{population}.filt.vcf"
    params:
        snpsift = "/home/hmbloom/snpEff/SnpSift.jar",
        config = "/home/hmbloom/snpEff/snpEff.config"
    shell:
        '''
        java -Xmx8g -jar {params.snpsift} filter "ANN[*].IMPACT has 'HIGH'" {input.vcf} > {output.vcf}
        '''

rule pangolin:
    message: '''### Predicting Splice Site Mutations with Pangolin ###'''
    input: 
        vcf = "results/{dataset}/{population}/annotated/{population}.filt.vcf"
    output:
        vcf = "results/{dataset}/{population}/pangolin/{population}.pangolin.vcf.gz"
    params:
        db = "../../Tools/Pangolin/gencode.v38.annotation.db",
        ref = config["ref"]
    shell:
        '''
        pangolin {input.vcf} {params.ref} {params.db} {output.vcf}
        '''
