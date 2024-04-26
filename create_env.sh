mamba create --name snakemake snakemake samtools bwa fastp gatk4 bcftools tabix 
mamba activate snakemake
pip install snakemake-executor-plugin-cluster-generic