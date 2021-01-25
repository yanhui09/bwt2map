#--------------
configfile: "config.yaml"
#wildcard_constraints:
#        barcode="[0-9][0-9]"
#--------------
INPUT_DIR = config["fq_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
REF_FA = config["ref_fa"].rstrip("/")

BATCH, LIB_ID = glob_wildcards(INPUT_DIR + "/{batch, NXT\d+}_{lib_id, XIC\d+}")
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input: 
        OUTPUT_DIR + "/bwt2_index/ref.fa",
        expand(OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth", batch = BATCH, lib_id = LIB_ID)

rule bwt2_index:
    input:
        REF_FA
    output:
        OUTPUT_DIR + "/bwt2_index/ref.fa"
    conda: 
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/bwt2_index.log"
    shell:
        "bowtie2-build {input} {output}"

rule cleanup_fq_names:
    input: 
        INPUT_DIR + "/{batch}_{lib_id}"
    output:
        r1 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R2.fastq.gz"
    log:
        OUTPUT_DIR + "/logs/cleanup_fq_names/{batch}/{lib_id}.log"
    shell:
        """
        mv {input}/*R1*.fastq.gz {output.r1}
        mv {input}/*R2*.fastq.gz {output.r2}
        """

rule bwt2_map:
    input:  
        r1 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R2.fastq.gz"
    output: 
        OUTPUT_DIR + "/bwt2_map/{batch}/{lib_id}.bam"
    params:
        bwt2_index = OUTPUT_DIR + "/bwt2_index/ref.fa"
    conda:
        "envs/bwt2sam.yaml"
    log: 
        OUTPUT_DIR + "/logs/bwt2_map/{batch}/{lib_id}.log"
    threads: 8
    shell:
        "bowtie2 -p {threads} -x {params.bwt2_index} -1 {input.r1} -2 {input.r2} | samtools view -bS - > {output}"

rule depth_per_base:
    input:
        OUTPUT_DIR + "/bwt2_map/{batch}/{lib_id}.bam"
    output:
        OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth"
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/depth/{batch}/{lib_id}.log"
    shell:
        "samtools depth {input} > {output}"
