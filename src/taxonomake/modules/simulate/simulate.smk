import os.path
from taxonomake.modules.common import (
    config_sample_reads1, config_sample_reads2, config_sample_names,
    config_truth, config_genomes_list,
    config_taxonomy, config_readsim_bin
)

SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')
MANIFEST_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(workflow.snakefile))), 'pixi.toml')

THREADS = 8

rule simulate_paired_reads_rename:
    input:
        r1 = ["readsim_" + config["readsim"]["tool"] + f"/{sample}_1.fq.gz" for sample in config_sample_names(config)],
        r2 = ["readsim_" + config["readsim"]["tool"] + f"/{sample}_2.fq.gz" for sample in config_sample_names(config)]
    output:
        r1 = config_sample_reads1(config),
        r2 = config_sample_reads2(config)
    shell:
        f"python3 {SIM_SCRIPTS_DIR}/rename_all.py " \
        "-i {input.r1} {input.r2} " \
        "-o {output.r1} {output.r2}"

rule simulate_art_paired_reads_sample:
    output:
        r1 = "readsim_art/{sample}_1.fq.gz",
        r2 = "readsim_art/{sample}_2.fq.gz"
    input:
        truth=config_truth(config),
        genomes_list=config_genomes_list(config),
        taxonomy=config_taxonomy(config)
    params:
        art_bin=config_readsim_bin(config)
    threads: THREADS
    log: "logs/{sample}.log"
    shell:
        "mkdir -p readsim_art && " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e art " \
        f"python3 {SIM_SCRIPTS_DIR}/simulate_art.py " \
        "--art {params.art_bin} " \
        "--threads {threads} " \
        "--coverage-file {input.truth} " \
        "--genome-list {input.genomes_list} " \
        "--taxonomy {input.taxonomy} " \
        "--sample {wildcards.sample} " \
        "-1 {output.r1} " \
        "-2 {output.r2} " \
        "2> {log}"