import os.path

SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')
MANIFEST_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(workflow.snakefile))), 'pixi.toml')

def make_absolute(*paths):
    return os.path.normpath(os.path.join(config["configfiledir"], *paths))

def make_absolute_if_path(path):
    return make_absolute(path) if os.path.dirname(path) != "" else path

config:
    threads: 1

rule simulate_paired_reads_rename:
    input:
        r1 = ["readsim_" + config["readsim"]["tool"] + f"/{sample}_1.fq.gz" for sample in config["samples"]["names"]],
        r2 = ["readsim_" + config["readsim"]["tool"] + f"/{sample}_2.fq.gz" for sample in config["samples"]["names"]]
    output:
        r1 = [make_absolute(s) for s in config["samples"]["reads1"]],
        r2 = [make_absolute(s) for s in config["samples"]["reads2"]]
    shell:
        f"python3 {SIM_SCRIPTS_DIR}/rename_all.py " \
        "-i {input.r1} {input.r2} " \
        "-o {output.r1} {output.r2}"

rule simulate_art_paired_reads_sample:
    output:
        r1 = "readsim_art/{sample}_1.fq.gz",
        r2 = "readsim_art/{sample}_2.fq.gz"
    input:
        truth=make_absolute(config["truth"]),
        genomes_list=make_absolute(config["genomes_list"]),
        taxonomy=make_absolute(config["taxonomy"])
    params:
        art_bin=make_absolute_if_path(config["readsim"]["bin"])
    threads: config["threads"]
    log: "logs/{sample}.log"
    shell:
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