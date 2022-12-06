import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

##### config file #####

validate(config, schema="../schemas/config.schema.yaml")

##### sample sheets #####

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], dtype=str, sep="\t").set_index(
    ["sample", "sequencing_type", "unit"], drop=False
)
validate(units, schema="../schemas/units.schema.yaml")

contigs = [c for c in range(1, 23)]
contigs.extend(["X", "Y"])

##### constraining wildcards #####


wildcard_constraints:
    pair="|".join(samples[samples.type == "tumor"]["sample"]),
    sample="|".join(samples["sample"]),
    caller="|".join(["freebayes", "delly"]),
    event="somatic|germline|complete",


### Output generation ###


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_final_output():
    if config["epitope_prediction"]["activate"]:
        final_output = expand(
            "results/neoantigens/{mhc}/{S.sample}.{S.sequencing_type}.xlsx",
            S=units.loc[samples[samples.type == "tumor"]["sample"]]
            .drop_duplicates(["sample", "sequencing_type"])
            .itertuples(),
            mhc=list(
                filter(
                    None,
                    [
                        "netMHCpan" if is_activated("affinity/netMHCpan") else None,
                        "netMHCIIpan" if is_activated("affinity/netMHCIIpan") else None,
                    ],
                )
            ),
        )
    else:
        if config["HLAtyping"]["HLA_LA"]["activate"]:
            final_output = expand(
                [
                    "results/optitype/{sample}/hla_alleles_{sample}.tsv",
                    "results/HLA-LA/hlaI_{sample}.tsv",
                    "results/HLA-LA/hlaII_{sample}.tsv",
                ],
                sample=samples["sample"],
            )
        else:
            final_output = expand(
                "results/optitype/{sample}/hla_alleles_{sample}.tsv",
                sample=samples["sample"],
            )
    return final_output


def get_fusion_output():
    if config["fusion"]["arriba"]["activate"]:
        fusion_output = expand(
            "results/fusion/arriba/{sample}.fusions.tsv",
            sample=units[units["sequencing_type"] == "RNA"]["sample"],
        )
    else:
        fusion_output = []
    return fusion_output


def get_tmb_targets():
    if is_activated("tmb"):
        return expand(
            "results/plots/tmb/{group}.{mode}.svg",
            group=samples[(samples.type == "tumor")]["sample"],
            mode=config["tmb"].get("mode", "curve"),
        )
    else:
        return []


caller = list(
    filter(None, ["freebayes" if is_activated("calling/freebayes") else None])
)

### helper functions ###

## alignment ##


def get_cutadapt_input(wildcards):
    # unit = units.loc[wildcards.sample].loc[wildcards.unit].loc[wildcards.seqtype]
    unit = units.loc[wildcards.sample].loc[wildcards.seqtype,wildcards.unit] #modified OG

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{T}/{U}.fq1.fastq{E}".format(
            S=unit.sample, U=unit.unit, T=unit.sequencing_type, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{T}/{U}.{{read}}.fastq{E}".format(
                S=wildcards.sample, U=unit.unit, T=unit.sequencing_type, E=ending
            ),
            # read=["fq1", "fq2"],
            read=["R1", "R2"], #og
        )


def get_cutadapt_pipe_input(wildcards):
    # pattern = (
    #     units.loc[wildcards.sample]
    #     .loc[wildcards.seqtype]
    #     .loc[wildcards.unit, wildcards.fq]
    # )
    # files = list(sorted(glob.glob(pattern)))
    # assert len(files) > 0, "no files found at {}".format(pattern)
    # return files
    return (units.loc[wildcards.sample].loc[wildcards.seqtype].loc[wildcards.unit, ["fq1","fq2"]].to_list())


def is_paired_end(sample, seqtype):
    sample_units = units.loc[sample].loc[seqtype]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {} and sequencing type {}, must be all paired end or all single end".format(
        sample, seqtype
    )
    return all_paired


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}/{seqtype}/{unit}_{read}.fastq.gz",
            # unit=units.loc[wc.seqtype].loc[wc.sample, "unit_name"],
            unit=units.loc[wc.sample].loc[wc.seqtype, "unit"], #modified by OG
            sample=wc.sample,
            read=wc.read,
            seqtype=wc.seqtype,
        )
    unit = units.loc[wc.sample].loc[wc.seqtype]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand(
            "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return units.loc[wc.sample].loc[wc.seqtype, fq].tolist()


def get_map_reads_input(wildcards):
    if is_paired_end(wildcards.sample, "DNA"):
        return [
            "results/merged/DNA/{sample}_R1.fastq.gz",
            "results/merged/DNA/{sample}_R2.fastq.gz",
        ]
    return "results/merged/DNA/{sample}_single.fastq.gz"


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample, platform=samples.loc[wildcards.sample, "platform"]
    )


def get_recalibrate_quality_input(wildcards, bai=False):
    ext = ".bai" if bai else ""
    if is_activated("remove_duplicates"):
        return "results/dedup/{}.sorted.bam{}".format(wildcards.sample, ext)
    else:
        return "results/mapped/{}.sorted.bam{}".format(wildcards.sample, ext)


## HLA Typing ##


# def get_optitype_reads_input(wildcards):
#     if is_activated("HLAtyping/optitype_prefiltering"):
#         if is_paired_end(wildcards.sample, "DNA"):
#             return expand(
#                 "results/razers3/fastq/{sample}_{fq}.fished.fastq",
#                 sample=wildcards.sample,
#                 fq=["R1", "R2"],
#             )
#         return "results/razers3/fastq/{sample}_single.fastq"
#     else:
#         return get_map_reads_input(wildcards)


# def get_oncoprint_batch(wildcards):
#     if wildcards.batch == "all":
#         groups = samples[samples["type"] == "tumor"]["sample"].unique()
#     else:
#         groups = samples.loc[
#             samples[config["oncoprint"]["stratify"]["by-column"]] == wildcards.batch,
#             "group",
#         ].unique()
#     return expand(
#         "results/merged-calls/{group}.{{event}}.fdr-controlled.bcf", group=groups
#     )


# ## variant calls ##


# def get_annotated_bcf(wildcards):
#     selection = ".annotated"
#     return "results/calls/{pair}.{scatteritem}{selection}.bcf".format(
#         pair=wildcards.pair, selection=selection, scatteritem=wildcards.scatteritem
#     )


# def get_scattered_calls(ext=".bcf"):
#     def inner(wildcards):
#         return expand(
#             "results/calls/{{pair}}.{caller}.{{scatteritem}}.sorted{ext}",
#             caller=caller,
#             ext=ext,
#         )

#     return inner


# def get_fdr_control_params(wildcards):
#     query = config["calling"]["fdr-control"]["events"][wildcards.event]
#     threshold = query.get(
#         "threshold", config["calling"]["fdr-control"].get("threshold", 0.05)
#     )
#     events = query["varlociraptor"]
#     return {"threshold": threshold, "events": events}


# def get_pair_variants(wildcards, index):
#     if index:
#         ext = ".csi"
#     else:
#         ext = ""
#     variants = [
#         "results/strelka/somatic/{}/results/variants/somatic.complete.tumor.bcf{}".format(
#             wildcards.sample, ext
#         )
#     ]
#     variants.append(
#         "results/strelka/germline/{}/results/variants/variants.reheader.bcf{}".format(
#             get_normal(wildcards), ext
#         )
#     )
#     return variants


# def get_pair_observations(wildcards):
#     return expand(
#         "results/observations/{pair}/{sample}.{caller}.{scatteritem}.bcf",
#         caller=wildcards.caller,
#         pair=wildcards.pair,
#         scatteritem=wildcards.scatteritem,
#         sample=get_paired_samples(wildcards),
#     )


# def get_merge_input(ext=".bcf"):
#     def inner(wildcards):
#         return expand(
#             "results/calls/{{pair}}.{vartype}.{{event}}.fdr-controlled{ext}",
#             ext=ext,
#             vartype=["SNV", "INS", "DEL", "MNV"],
#             filter=config["calling"]["fdr-control"]["events"][wildcards.event],
#         )

#     return inner


# def get_pair_aliases(wildcards):
#     return [
#         samples.loc[samples.loc[wildcards.pair, "matched_normal"], "type"],
#         samples.loc[wildcards.pair, "type"],
#     ]


def get_tabix_params(wildcards):
    if wildcards.format == "vcf":
        return "-p vcf"
    if wildcards.format == "txt":
        return "-s 1 -b 2 -e 2"
    raise ValueError("Invalid format for tabix: {}".format(wildcards.format))


# ## RNA ##


# def get_quant_reads_input(wildcards):
#     if is_paired_end(wildcards.sample, "RNA"):
#         return [
#             "results/merged/RNA/{sample}_R1.fastq.gz",
#             "results/merged/RNA/{sample}_R2.fastq.gz",
#         ]
#     return "results/merged/RNA/{sample}_single.fastq.gz"


# def kallisto_params(wildcards, input):
#     extra = config["params"]["kallisto"]
#     if len(input.fastq) == 1:
#         extra += " --single"
#         extra += (
#             " --fragment-length {unit.fragment_len_mean} " "--sd {unit.fragment_len_sd}"
#         ).format(unit=units.loc[(wildcards.sample, wildcards.unit)])
#     else:
#         extra += " --fusion"
#     return extra


## helper functions ##


def get_paired_samples(wildcards):
    return [
        samples.loc[(wildcards.pair), "matched_normal"],
        samples.loc[wildcards.pair, "sample"],
    ]


def get_paired_bams(wildcards):
    return expand(
        "results/recal/{sample}.sorted.bam", sample=get_paired_samples(wildcards)
    )


def get_paired_bais(wildcards):
    return expand(
        "results/recal/{sample}.sorted.bam.bai", sample=get_paired_samples(wildcards)
    )


def get_normal(wildcards):
    return samples.loc[(wildcards.sample), "matched_normal"]


def get_reads(wildcards):
    return get_seperate(wildcards.sample, wildcards.group)


def get_seperate(sample, group):
    return units.loc[(sample, "DNA"), "fq{}".format(str(group))]


def get_proteome(wildcards):
    return expand(
        "results/microphaser/fasta/germline/{normal}/{mhc}/reference_proteome.bin",
        normal=get_normal(wildcards),
        mhc=wildcards.mhc,
    )


def get_alleles_MHCI(wildcards):
    if wildcards.group == "wt":
        return "results/optitype/{S}/hla_alleles_{S}.tsv".format(
            S=get_normal(wildcards)
        )
    else:
        return "results/optitype/{S}/hla_alleles_{S}.tsv".format(S=wildcards.sample)


def get_alleles_MHCII(wildcards):
    if wildcards.group == "wt":
        return "results/HLA-LA/hlaI_{S}.tsv".format(S=get_normal(wildcards))
    else:
        return "results/HLA-LA/hlaI_{S}.tsv".format(S=wildcards.sample)


def get_normal_bam(wildcards):
    return expand("results/recal/{normal}.sorted.bam", normal=get_normal(wildcards))


def get_normal_bai(wildcards):
    return expand("results/recal/{normal}.sorted.bam.bai", normal=get_normal(wildcards))
