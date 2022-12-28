rule HLA_LA:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai",
        index="resources/graphs/PRG_MHC_GRCh38_withIMGT/serializedGRAPH",
    output:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt",
    threads: 16
    log:
        "logs/HLA-LA/{sample}.log",
    params:
        graph=lambda w, input: os.path.basename(os.path.dirname(input.index)),
        graphdir=lambda w, input: os.path.dirname(os.path.dirname(input.index)),
    conda:
        "../envs/hla_la.yaml"
    shell:
        "HLA-LA.pl --bam {input.bam} --sampleID {wildcards.sample} --graph {params.graph} --customGraphDir {params.graphdir} --workingDir results/HLA-LA/output --maxThreads {threads} > {log} 2>&1"


rule parse_HLA_LA:
    input:
        "results/HLA-LA/output/{sample}/hla/R1_bestguess_G.txt",
    output:
        report(
            "results/HLA-LA/hlaI_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
        report(
            "results/HLA-LA/hlaII_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(HLA-LA)",
        ),
    log:
        "logs/parse-HLA-LA/{sample}.log",
    script:
        "../scripts/parse_HLA_types.py"

rule map_hla_reads: #OG
    input:
        # reads=get_map_reads_input,
        reads=get_quant_reads_input, # RNA
        idx=rules.HLA_index.output,
    output:
        "results/HLA_mapped/{sample}.hla.sorted.bam",
    log:
        "logs/bwa_mem_hla/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 16
    wrapper:
        "0.56.0/bio/bwa/mem"

# rule map_hla_reads: #OG
#     input:
#         reads=get_map_reads_input,
#         idx=rules.HLA_index2.output,
#     output:
#         "results/HLA_mapped/{sample}.hla.sorted.bam",
#     log:
#         "logs/bwa_mem_hla/{sample}.log",
#     params:
#         index=lambda w, input: os.path.splitext(input.idx[0])[0],
#         extra=get_read_group,
#         sort="samtools",
#         sort_order="coordinate",
#     threads: 8
#     wrapper:
#         "v1.21.1/bio/bwa-mem2/mem"

# rule filter_hla_reads: #OG
#     input:
#         "results/HLA_mapped/{sample}.hla.sorted.bam",
#     output:
#         bam="results/HLA_mapped/{sample}.hla.filtered.sorted.bam",
#         idx="results/HLA_mapped/{sample}.hla.filtered.sorted.bai",
#     log:
#         "logs/bwa_mem_hla/{sample}.filter.log",
#     params:
#         extra="-F 0x4",  
#     threads: 2
#     wrapper:
#         "v1.21.0/bio/samtools/view"

rule separate_hla_reads: # OG
    input:
        "results/HLA_mapped/{sample}.hla.sorted.bam",
    output:
        temp("results/HLA_mapped/{sample}.R1.fq.gz"),
        temp("results/HLA_mapped/{sample}.R2.fq.gz"),
    params:
        sort = "-m 4G",
        bam2fq = "-n -F 2308" #OG removed
    threads:  # Remember, this is the number of samtools' additional threads
        16     # At least 2 threads have to be requested on cluster sumbission.
              # Thus, this value - 2 will be sent to samtools sort -@ argument.
    wrapper:
        "0.61.0/bio/samtools/bam2fq/separate"


rule razers3:
    input:
        reads="results/merged/DNA/{sample}_{fq}.fastq.gz",
    output:
        bam="results/razers3/bam/{sample}_{fq}.bam",
    threads: 16
    log:
        "logs/razers3/{sample}_{fq}.log",
    params:
        genome=config["HLAtyping"]["optitype_data"],
        extra=config["params"]["razers3"],
    wrapper:
        "0.61.0/bio/razers3"


rule bam2fq:
    input:
        "results/razers3/bam/{sample}_{fq}.bam",
    output:
        "results/razers3/fastq/{sample}_{fq}.fished.fastq",
    params:
        "",
    log:
        "logs/razers3-bam2fq/{sample}-{fq}.log",
    threads: 1
    wrapper:
        "0.61.0/bio/samtools/bam2fq/interleaved"


rule OptiType: #og w rna
    input:
        reads=get_optitype_reads_input,
    output:
        multiext(
            "results/optitype/{sample}/{sample}", "_coverage_plot.pdf", "_result.tsv"
        ),
    log:
        "logs/optitype/{sample}.log",
    params:
        extra=config["params"]["optitype"],
        sequencing_type="dna", #og
    # resources:
        # mem_mb=30000
    conda:
        "../envs/optitype.yaml"
    shell:
        "OptiTypePipeline.py --input {input} --outdir results/optitype/{wildcards.sample} --prefix {wildcards.sample} --dna"

# rule OptiType: #og original
#     input:
#         reads=get_optitype_reads_input,
#     output:
#         multiext(
#             "results/optitype/{sample}/{sample}", "_coverage_plot.pdf", "_result.tsv"
#         ),
#     log:
#         "logs/optitype/{sample}.log",
#     params:
#         extra=config["params"]["optitype"],
#         sequencing_type="dna",
#         sequencing_type="rna", #og
#     wrapper:
#         "0.63.0/bio/optitype"


rule parse_Optitype:
    input:
        "results/optitype/{sample}/{sample}_result.tsv",
    output:
        report(
            "results/optitype/{sample}/hla_alleles_{sample}.tsv",
            caption="../report/HLA_Types.rst",
            category="HLA-Typing(Optitype)",
        ),
    log:
        "logs/parse-optitype/{sample}.log",
    shell:
        "cut {input} -f2-7 | awk 'NR == 1 {{print}} NR>1 {{for (i = 1; i<=6; ++i) sub(/^/, \"&HLA-\", $i); print}}' "
        '| sed -e s/[*,:]//g | sed "s/ /\t/g" > {output} 2> {log}'
