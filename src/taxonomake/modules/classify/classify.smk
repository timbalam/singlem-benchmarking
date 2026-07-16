import os.path

SIM_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')
MANIFESTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(workflow.snakefile)))), 'manifests')

def make_absolute(*paths):
    return os.path.normpath(os.path.join(config["configfiledir"], *paths))

gtdbtk_dir = make_absolute(config["gtdbtk"]["dir"])

# release 207
rule download_gtdbtk_r207_data:
    output:
        make_absolute(config["gtdbtk"]["dir"], '/gtdbtk_r207_v2_data.tar.gz')
    log:
        "logs/gtdbtk_r207_data-download.log"
    shell:
        "bash -c "\
        "'cd " + config["gtdbtk_dir"] + " && "\
        "wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz' &> {log}"

rule extract_gtdbtk_r207_data:
    input:
        config["gtdbtk_dir"] + "/gtdbtk_r207_v2_data.tar.gz"
    output:
        directory(config[["gtdbtk_data"]])
    log:
        "log/gtdbtk_r207_v2_data-extract.log"
    shell:
        "bash -c " \
        "'cd " + config["gtdbtk_dir"] + " && " \
        "tar -xzf gtdbtk_r207_v2_data.tar.gz && " +
        "mv release207_v2 " + config["gtdbtk_data"] + "' &> {log}"

rule gtdbtk_r207_identify:
    input:
        genomes_list=config[["genomes_list"]],
        data_path=config[["gtdbtk_data"]]
    output:
        output_dir=directory(config[["gtdbtk_dir"]] + "/identify")
    log:
        "logs/gtdbtk_r207_identify.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e gtdbtk-r207 " \
        "gtdbtk identify --batchfile {input.genomes_list} " \
        "--out_dir {output.output_dir} " \
        "--extension .fasta " \
        "&> {log}"

rule gtdbtk_r207_align:
    input:
        id=config[["gtdbtk_dir"]] + "/identify",
        data_path=config[["gtdbtk_data"]]
    output:
        output_dir=directory(config[["gtdbtk_dir"]] + "/align")
    log:
        "logs/gtdbtk_r207_align.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        f"pixi run --manifest-path {MANIFEST_PATH} -e gtdbtk-r207 " \
        "gtdbtk align --identify_dir {input.id} " \
        "--out_dir {output.output_dir} " \
        "&> {log}"

rule gtdbtk_r207_classify:
    input:
        genomes_file = config[["genomes_file"]],
        al = config[["gtdbtk_dir"]] + "/align",
        data_path = config[["gtdbtk_data"]]
    output:
        output_dir = directory(config[["gtdbtk_dir"]] + "/classify"),
        done=touch(config[["gtdbtk_dir'"]] + "/classify.done")
    log:
        "logs/gtdbtk_r207_classify.log"
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
