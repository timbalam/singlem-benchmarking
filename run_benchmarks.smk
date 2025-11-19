from os.path import join

tools = ['singlem', 'sylph']
datasets = [f'marine{i}' for i in range(10)]
benchmark_dirs = ['5_novelty']
num_threads = config['benchmarking_threads']

singlem_bin = "singlem"
singlem_metapackage = "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"
sylph_db = "tool_reference_data/gtdb_database.syldb"

#####################################################################

rule all:
    input:
        expand("{bench_dir}/output_{tool}/opal/{sample}.opal_report",
               bench_dir = benchmark_dirs, sample = datasets, tool = tools)

rule generate_communities_5:
    input:
        "5_novelty/generate_communities.done"

rule generate_communities:
    input:
        [join("{bench_dir}", f'truth/marine{i}.finished') for i in range(10)],
        [join("{bench_dir}", f'reads/marine{i}.finished') for i in range(10)],
        [join("{bench_dir}", f'truth/marine{i}.condensed.biobox') for i in range(10)],
    output:
        done=touch("{bench_dir}/generate_communities.done")

rule generate_5_community_and_reads:
    input:
        gtdb_bac_metadata = './bac120_metadata_r207.tsv'
        gtdb_ar_metadata = './ar53_metadata_r207.tsv'
    output:
        r1="5_novelty/reads/{sample}.1.fq.gz",
        r2="5_novelty/reads/{sample}.2.fq.gz",
        condensed = "5_novelty/truth/{sample}.condensed",
        genomewise = "5_novelty/truth/{sample}.genomewise.csv",
        done = touch("5_novelty/truth/{sample}.finished"),
        done2 = touch("5_novelty/reads/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('marine', ''),
    threads: num_threads
    shell:
        "mkdir -p 5_novelty/truth 5_novelty/reads && " \
        "pixi run -e art python3 " \
        "5_novelty/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 5_novelty/coverage_definitions/coverage{params.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} --gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--genome-list 5_novelty/shadow_genome_paths.csv --output-condensed {output.condensed} " \
        "-1 5_novelty/reads/{wildcards.sample}.1.fq.gz " \
        "-2 5_novelty/reads/{wildcards.sample}.2.fq.gz " \
        "--output-genomewise-coverage {output.genomewise}"

rule truth_condensed_to_biobox:
    input:
        condensed = "{bench_dir}/truth/{sample}.condensed",
    output:
        biobox = "{bench_dir}/truth/{sample}.condensed.biobox"
    shell:
        "pixi run -e singlem_deps python3 " \
        "bin/condensed_profile_to_biobox.py --input-condensed-table {input.condensed} " \
        "--output-biobox {output.biobox}"

rule tool_condensed_to_biobox:
    input:
        profile = "{bench_dir}/output_{tool}/{tool}/{sample}.profile",
        truth = "{bench_dir}/truth/{sample}.condensed.biobox",
    output:
        biobox = "{bench_dir}/output_{tool}/biobox/{sample}.biobox"
    shell:
        "pixi run -e singlem_deps python3 " \
        "bin/condensed_profile_to_biobox.py --input-condensed-table {input.profile} " \
        "--output-biobox {output.biobox} --template-biobox {input.truth} "

rule opal:
    input:
        biobox = "{bench_dir}/{tool_output}/biobox/{sample}.biobox")
    params:
        output_dir = "{bench_dir}/{tool_output}",
        output_opal_dir = "{bench_dir}/{tool_output}/opal/{sample}.opal_output_directory",
        truth = join(truth_dir, "{sample}.condensed.biobox"),
    output:
        report="{bench_dir}/{tool_output}/opal/{sample}.opal_report",
        done=touch("{bench_dir}/{tool_output}/opal/{sample}.opal_report.done")
    shell:
        "pixi run -e opal " \
        "opal.py -g {params.truth} -o {params.output_opal_dir} {input.biobox} || echo 'expected opal non-zero exit status'; mv {params.output_opal_dir}/results.tsv {output.report} && rm -rf {params.output_opal_dir}"


###############################################################################################
###############################################################################################
###############################################################################################
#########
######### tool-specific rules - singlem first

rule singlem_run_pipe:
    input:
        r1="{bench_dir}/reads/{sample}_1.fastq.gz",
        r2="{bench_dir}/reads/{sample}_2.fastq.gz",
        db=singlem_metapackage,
    output:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        done=touch("{bench_dir}/output_singlem/singlem/{sample}.sma.done")
    threads:
        num_threads
    log:
        "{bench_dir}/output_singlem/logs/singlem/{sample}.log"
    shell:
        "pixi run -e singlem_deps" \
        "{singlem_bin} pipe --threads {threads} -1 {input.r1} -2 {input.r2} " \
        "--archive-otu-table {output.report} --metapackage {input.db} &> {log}"

rule singlem_run_condense:
    input:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        done="{bench_dir}/output_singlem/singlem/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem/singlem/{sample}.profile",
        done=touch("{bench_dir}/output_singlem/singlem/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem/logs/singlem/{sample}.log"
    shell:
        "pixi run -e singlem_deps " \
        "{singlem_bin} condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --metapackage {input.db} &> {log}"


###############################################################################################
###############################################################################################
###############################################################################################
######### sylph

rule sylph_run:
    input:
        r1 = "{bench_dir}/reads/{sample}_1.fastq.gz",
        r2 = "{bench_dir}/reads/{sample}_2.fastq.gz",
        db = sylph_package
    output:
        report="{bench_dir}/output_sylph/sylph/{sample}.tsv",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.done")
    threads: num_threads
    resources:
        mem_mb=32000
    log:
        "{bench_dir}/output_sylph/logs/sylph/{sample}.log"
    shell:
        "pixi run -e sylph " \
        "sylph profile {input.db} -1 {input.r1} -2 {input.r2} -t {threads} " \
        "> {output.report} 2> {log}"

rule sylph_report_to_condensed:
    input:
        report = "{bench_dir}/output_sylph/sylph/{sample}.tsv",
        gtdb_bac_tax = "./bac120_taxonomy_r207.tsv",
        gtdb_ar_tax = "./ar53_taxonomy_r207.tsv",
    output:
        profile = "{bench_dir}/output_sylph/sylph/{sample}.profile"
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.profile.done")
    shell:
        "pixi shell -e singlem_deps python3 " \
        "bin/sylph_to_condensed.py --sylph-genome {input.report} " \
        "--sample {wildcards.sample} " \
        "--bac-tax {input.gtdb_bac_tax} " \
        "--arc-tax {input.gtdb_ar_tax} > {output.profile}"
