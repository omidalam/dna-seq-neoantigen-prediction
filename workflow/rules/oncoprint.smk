rule build_oncoprint_table:
    input:
        bcf=get_oncoprint_batch,
    output:
        "plots/oncoprint/{batch}.{event}.tsv",
    log:
        "logs/oncoprint/{batch}.{event}.table.log",
    conda:
        "../envs/oncoprinttable.yaml"
    script:
        "../scripts/build_oncoprint_matrix.py"


rule plot_oncoprint:
    input:
        "plots/oncoprint/{batch}.{event}.tsv",
    output:
        report(
            "plots/oncoprint/{batch}.{event}.pdf",
            category="Oncoprint",
            caption="../report/oncoprint.rst",
        ),
    log:
        "logs/oncoprint/{batch}.{event}.plot.log",
    conda:
        "../envs/oncoprint.yaml"
    script:
        "../scripts/oncoprint.R"
