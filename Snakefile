#--------------
configfile: "config.yaml"
#--------------
INPUT_DIR = config["fq_dir"].rstrip("/")
OUTPUT_DIR = config["results_dir"].rstrip("/")
REF_FA = config["ref_fa"].rstrip("/")

BATCH, LIB_ID = glob_wildcards(INPUT_DIR + "/{batch, .*\d+}_{lib_id, .*\d+}")
#---------------

# Allow users to fix the underlying OS via singularity.
#singularity: "docker://continuumio/miniconda3"

rule all:
    input: 
        #expand(OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth", batch = BATCH, lib_id = LIB_ID),
        expand(OUTPUT_DIR + "/depth_analysis/{batch}", batch = BATCH)
    
rule bwt2_build:
    input:
        REF_FA
    output:
        o1 = expand(OUTPUT_DIR + "/bwt2_build/ref" + '.{i}.bt2', i=[1,2,3,4]),
        o2 = expand(OUTPUT_DIR + "/bwt2_build/ref" + '.rev.{i}.bt2', i=[1,2])
    params:
        prefix = OUTPUT_DIR + "/bwt2_build/ref"
    conda: 
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/bwt2_build.log"
    shell:
        "bowtie2-build {input} {params.prefix} >> {log} 2>&1"

rule cleanup_fq_names:
    input: 
        INPUT_DIR + "/{batch}_{lib_id}"
    output:
        r1 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R2.fastq.gz"
    shell:
        """
        mv {input}/*R1*.fastq.gz {output.r1}
        mv {input}/*R2*.fastq.gz {output.r2}
        """

rule bwt2_map:
    input:  
        r1 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R1.fastq.gz",
        r2 = INPUT_DIR + "/{batch}_{lib_id}/{batch}_{lib_id}_R2.fastq.gz",
        build_out = rules.bwt2_build.output
    output: 
        OUTPUT_DIR + "/bwt2_map/{batch}/{lib_id}.bam"
    params:
        bwt2_build = OUTPUT_DIR + "/bwt2_build/ref"
    conda:
        "envs/bwt2sam.yaml"
    log: 
        OUTPUT_DIR + "/logs/bwt2_map/{batch}/{lib_id}.log"
    threads: 8
    shell:
        "(bowtie2 -p {threads} -x {params.bwt2_build} -1 {input.r1} -2 {input.r2} | samtools view -bS - > {output}) 2>> {log}"

rule samtools_sort:
    input:
        OUTPUT_DIR + "/bwt2_map/{batch}/{lib_id}.bam"
    output:
        OUTPUT_DIR + "/sorted_bams/{batch}/{lib_id}.bam"
    params:
        temp_prefix = OUTPUT_DIR + "/sorted_bams/{batch}/{lib_id}"
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/sorted_bams/{batch}/{lib_id}.log"
    shell:
        "(samtools sort -T {params.temp_prefix} -O bam {input} > {output}) 2>> {log}"
        
rule depth_per_base:
    input:
        OUTPUT_DIR + "/sorted_bams/{batch}/{lib_id}.bam"
    output:
        OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth"
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/depth/{batch}/{lib_id}.log"
    shell:
        "(samtools depth -a {input} > {output}) 2>> {log}"

rule depth_analysis:
    input:
        expand(OUTPUT_DIR + "/depth/{batch}/{lib_id}.depth", batch = BATCH, lib_id = LIB_ID)
    output:
        directory(OUTPUT_DIR + "/depth_analysis/{batch}")
    params:
        depth_dir = OUTPUT_DIR + "/depth/{batch}",
        ref = OUTPUT_DIR + "/depth/{batch}/XIC16.depth"
    conda:
        "envs/bwt2sam.yaml"
    log:
        OUTPUT_DIR + "/logs/depth_analysis_{batch}.log"
    shell:
        """
        mkdir {output}
        scripts/plot_depth.py -d {params.depth_dir} -r {output} -f {params.ref} -w 500
        """
