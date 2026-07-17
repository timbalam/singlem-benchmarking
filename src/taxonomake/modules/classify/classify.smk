import os.path
from taxonomake.modules.common import (
    config_gtdbtk_dir, config_gtdbtk_data, config_genomes_list
)

SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')
MANIFEST_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(workflow.snakefile))), 'pixi.toml')

# release 207
rule download_gtdbtk_r207_data:
    output:
        config_gtdbtk_dir(config, 'gtdbtk_r207_v2_data.tar.gz')
    log:
        "logs/gtdbtk_r207_data-download.log"
    shell:
        "bash -c "\
        "'cd " + config_gtdbtk_dir(config) + " && "\
        "wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz' &> {log}"

rule extract_gtdbtk_r207_data:
    input:
        config_gtdbtk_dir("gtdbtk_r207_v2_data.tar.gz")
    output:
        directory(config_gtdbtk_data(config))
    log:
        "log/gtdbtk_r207_v2_data-extract.log"
    shell:
        "bash -c " \
        "'cd " + config_gtdbtk_dir(config) + " && " \
        "tar -xzf gtdbtk_r207_v2_data.tar.gz && " +
        "mv release207_v2 " + config_gtdbtk_data(config) + "' &> {log}"

rule gtdbtk_r207_identify:
    input:
        genomes_list=config_genomes_list(config),
        data_path=config_gtdbtk_data(config)
    output:
        output_dir=directory(config_gtdbtk_dir(config, "identify"))
    log:
        "logs/gtdbtk/r207/identify.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e gtdbtk-r207 " \
        "gtdbtk identify --batchfile {input.genomes_list} " \
        "--out_dir {output.output_dir} " \
        "--extension .fasta " \
        "&> {log}"

rule gtdbtk_r207_align:
    input:
        id=config_gtdbtk_dir(config, "identify"),
        data_path=config_gtdbtk_data(config)
    output:
        output_dir=directory(config_gtdbtk_dir(config, "align"))
    log:
        "logs/gtdbtk/r207/align.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e gtdbtk-r207 " \
        "gtdbtk align --identify_dir {input.id} " \
        "--out_dir {output.output_dir} " \
        "&> {log}"

rule gtdbtk_r207_classify:
    input:
        genomes_file = config_genomes_file(config),
        al = config_gtdbtk_dir(config, "align"),
        data_path = config_gtdbtk_data(config)
    output:
        output_dir = directory(config_gtdbtk_dir(config, "classify")),
        done=touch(config_gtdbtk_dir(config, "classify.done"))
    log:
        "logs/gtdbtk/r207/classify.log"
    resources:
        mem_mb=64000
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e gtdbtk-r207 " \
        "gtdbtk classify --batchfile {input.genomes_file} " \
        "--align_dir {input.al} " \
        "--extension .fasta " \
        #"--scratch_data " \
        "--out_dir {output.output_dir} " \
        "&> {log}"
