SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')

rule simulate_art_paired_reads:
    input:
        config["reads1"],
        config["reads2"]

rule simulate_art_paired_reads_rename:
    input:
        r1 = [f"{sample}_1.fq.gz" for sample in config["samples"]],
        r2 = [f"{sample}_2.fq.gz" for sample in config["samples"]]
    output:
        r1 = config["reads1"],
        r2 = config["reads2"]
    shell:
        f"python3 {SIM_SCRIPTS_DIR}/rename_all.py " \
        "-i {input.r1} {input.r2} " \
        "-o {output.r1} {output.r2}"

rule simulate_art_paired_reads_sample:
    output:
        r1 = "{sample}_1.fq.gz",
        r2 = "{sample}_2.fq.gz"
    input:
        coverage_file=config["coverage_file"],
        genomes_list=config["genomes_list"]
    params:
        art_bin=config["art_bin"]
    threads: config["threads"]
    log: "logs/{sample}.log"
    shell:
        "pixi run -e art " \
        f"python3 {SIM_SCRIPTS_DIR}/simulate_art.py " \
        "--art {params.art_bin} " \
        "--threads {threads} " \
        "--coverage-file {input.coverage_file} " \
        "--genome-list {input.genomes_list} " \
        "--sample {wildcards.sample} " \
        "-1 {output.r1} " \
        "-2 {output.r2} " \
        "2> {log}"