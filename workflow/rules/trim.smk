rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "0.56.0/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{seqtype}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/cutadapt/{sample}-{seqtype}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 0  # this does not need CPU
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed/{sample}/{seqtype}/{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}/{seqtype}/{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}/{seqtype}/{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}-{seqtype}-{unit}.log",
    params:
        others=config["params"]["cutadapt"],
        adapters=lambda w: str(
            units.loc[w.sample].loc[w.seqtype].loc[w.unit, "adapters"]
        ),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed/{sample}/{seqtype}/{unit}.single.fastq.gz",
        qc="results/trimmed/{sample}/{seqtype}/{unit}.single.qc.txt",
    log:
        "logs/cutadapt/{sample}-{seqtype}-{unit}.log",
    params:
        others=config["params"]["cutadapt"],
        adapters_r1=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "0.59.2/bio/cutadapt/se"


rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "results/merged/{seqtype}/{sample}_{read}.fastq.gz",
        # OG Razers3 can't accept fake gzip files. They really should be gzip.
        # "results/merged/{seqtype}/{sample}_{read}.fastq", #OG 
    log:
        "logs/merge-fastqs/{seqtype}_{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
        seqtype="DNA|RNA",
    shell:
        "cat {input} > {output} 2> {log}"
