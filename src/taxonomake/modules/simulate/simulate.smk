SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')

rule simulate_art_reads:
    output:
        r1=config["reads1"]
        r2=config["reads2"]
    input:
        coverage_file=config["coverage_file"],
        genomes_file=config["genomes_file"]
    params:
        art_bin=config["art_bin"]
    threads: config["threads"]
    shell:
        "mkdir -p 9_related/truths 9_related/local_reads && " \
        "pixi run -e art " \
        "python3 {SIM_SCRIPTS_DIR}/simulate_art.py --art {params.art_bin} " \
        "--threads {threads} " \
        "--coverage-file {input.coverage_file} " \
        "--genome-list {input.genome_list} " \
        "-1 {output.r1} " \
        "-2 {output.r2} " \
        "2> {log}"