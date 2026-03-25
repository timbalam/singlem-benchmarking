
singlem_metapackage = "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"

rule all:
    input:
        expand("5_novelty/debug/div{div}top{top}/{sample}_hits.tsv",
               div = [2, 3, 4, 6, 8], top = [2, 5, 10],
               sample = ["marine0"])

rule one_test:
    input:
        "5_novelty/debug/div2top2/marine0_hits.tsv"

rule singlem_dev_renew:
    input:
        archive_otu_table = "5_novelty/output_singlem_dev/singlem_dev/{sample}.sma",
        metapackage = singlem_metapackage
    output:
        archive_otu_table = "5_novelty/debug/div{div}top{top}/{sample}.sma",
        done = touch("5_novelty/debug/div{div}top{top}/{sample}.sma.done")
    shell:
        "mkdir -p 5_novelty/debug/div{wildcards.div}top{wildcards.top} && " \
        "pixi run -e singlem-dev singlem renew " \
        "--input-archive-otu-table {input.archive_otu_table} " \
        "--archive-otu-table {output.archive_otu_table} " \
        "--metapackage {input.metapackage} " \
        "--max-species-divergence {wildcards.div} " \
        "--diamond-top {wildcards.top}"

rule write_best_hits:
    input:
        archive_otu_table = "5_novelty/debug/div{div}top{top}/{sample}.sma",
        metapackage = singlem_metapackage,
        script = "singlem/extras/write_archive_best_hits_table.py"
    output:
        hits_table = "5_novelty/debug/div{div}top{top}/{sample}_hits.tsv",
        done = touch("5_novelty/debug/div{div}top{top}/{sample}.done")
    shell:
        "pixi run -e singlem-dev python3 " \
        "{input.script} " \
        "--input-archive-otu-table {input.archive_otu_table} " \
        "--metapackage {input.metapackage} " \
        "--hits-table {output.hits_table}"