from os.path import join

datasets_bench5 = [f'marine{i}' for i in range(1)]

singlem_metapackage = "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"
sylph_package = "tool_reference_data/gtdb_database.syldb"

datasets_bench6 = ['SRR22388335', 'SRR9650389', 'SRR16352837', 'SRR16352839', 'SRR17498764',
                   'SRR22870123', 'SRR23961386', 'SRR24982124', 'SRR6201989', 'SRR5264410',
                   'SRR6869034', 'SRR5264435']

singlem_r226_metapackage = 'tool_reference_data/S5.4.0.GTDB_r226.metapackage_20250331.smpkg'
sylph_r226_package = 'tool_reference_data/gtdb-r226-c200-dbv1.syldb'
gtdb_r226_bac120_tax = 'tool_reference_data/bac120_taxonomy_r226.tsv'
gtdb_r226_ar53_tax = 'tool_reference_data/ar53_taxonomy_r226.tsv'

#####################################################################

rule all_bench5:
    input:
        expand("5_novelty/output_{tool}/opal/{sample}.opal_report",
               sample = datasets_bench5, tool = ['singlem', 'sylph', 'singlem_dev'])

rule all_bench6:
    input:
        expand("6_host_assocs/output_{tool}/opal/{sample}.opal_report",
               sample = datasets_bench6, tool = ['singlem_r266', 'sylph_r266'])

rule generate_communities_bench5:
    input:
        [f'5_novelty/truths/{sample}.finished' for sample in datasets_bench5],
        [f'5_novelty/local_reads/marine{sample}.finished' for sample in datasets_bench5],
        [f'5_novelty/truths/marine{sample}.condensed.biobox' for sample in datasets_bench5],
    output:
        done=touch("5_novelty/generate_communities.done")

rule download_bench6:
    input:
        [f'6_host_assocs/local_reads/{sample}_1.fastq.gz' for sample in datasets_bench6],
        [f'6_host_assocs/local_reads/{sample}_2.fastq.gz' for sample in  datasets_bench6]

rule generate_community_and_reads_bench5:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="5_novelty/local_reads/{sample}.1.fq.gz",
        r2="5_novelty/local_reads/{sample}.2.fq.gz",
        condensed = "5_novelty/truths/{sample}.condensed",
        #genomewise = "5_novelty/truths/{sample}.genomewise.csv",
        done = touch("5_novelty/truths/{sample}.finished"),
        done2 = touch("5_novelty/local_reads/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('marine', ''),
    log: "5_novelty/local_reads/{sample}.log"
    threads: 8
    shell:
        "mkdir -p 5_novelty/truths 5_novelty/local_reads && " \
        "pixi run -e art " \
        "python3 5_novelty/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 5_novelty/coverage_definitions/coverage{params.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 5_novelty/local_reads/{wildcards.sample}.1.fq.gz " \
        "-2 5_novelty/local_reads/{wildcards.sample}.2.fq.gz " \
        "2> {log}"

rule truth_condensed_to_biobox:
    input:
        condensed = "{bench_dir}/truths/{sample}.condensed",
    output:
        biobox = "{bench_dir}/truths/{sample}.condensed.biobox"
    shell:
        "pixi run -e singlem " \
        "python3 bin/condensed_profile_to_biobox.py --input-condensed-table {input.condensed} " \
        "--output-biobox {output.biobox}"

rule tool_condensed_to_biobox:
    input:
        profile = "{bench_dir}/output_{tool}/{tool}/{sample}.profile",
        truth = "{bench_dir}/truths/{sample}.condensed.biobox",
    output:
        biobox = "{bench_dir}/output_{tool}/biobox/{sample}.biobox"
    shell:
        "pixi run -e singlem " \
        "python3 bin/condensed_profile_to_biobox.py --input-condensed-table {input.profile} " \
        "--output-biobox {output.biobox} --template-biobox {input.truth} "

rule opal:
    input:
        biobox = "{bench_dir}/{tool_output}/biobox/{sample}.biobox"
    params:
        output_dir = "{bench_dir}/{tool_output}",
        output_opal_dir = "{bench_dir}/{tool_output}/opal/{sample}.opal_output_directory",
        truth = "{bench_dir}/truths/{sample}.condensed.biobox",
    output:
        report="{bench_dir}/{tool_output}/opal/{sample}.opal_report",
        done=touch("{bench_dir}/{tool_output}/opal/{sample}.opal_report.done")
    shell:
        "pixi run -e opal " \
        "opal.py -g {params.truth} -o {params.output_opal_dir} {input.biobox} || echo 'expected opal non-zero exit status'; mv {params.output_opal_dir}/results.tsv {output.report} && rm -rf {params.output_opal_dir}"

rule download_fastq:
    output:
        r1="{bench_dir}/local_reads/{sample}_1.fastq.gz",
        r2="{bench_dir}/local_reads/{sample}_2.fastq.gz",
        done=touch("{bench_dir}/local_reads/{sample}.done")
    shell:
        "pixi run -e kingfisher " \
        "kingfisher get -r {wildcards.sample} " \
        "--output_directory {bench_dir}/local_reads " \
        "-m ena-ftp prefetch "


###############################################################################################
###############################################################################################
###############################################################################################
#########
######### tool-specific rules - singlem first

rule singlem_run_pipe:
    input:
        r1="{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2="{bench_dir}/local_reads/{sample}.2.fq.gz",
        db=singlem_metapackage,
    output:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        done=touch("{bench_dir}/output_singlem/singlem/{sample}.sma.done")
    threads:
        8
    log:
        "{bench_dir}/output_singlem/logs/singlem/{sample}.log"
    shell:
        "pixi run -e singlem " \
        "singlem pipe --threads {threads} -1 {input.r1} -2 {input.r2} " \
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
        "pixi run -e singlem " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --metapackage {input.db} &> {log}"

rule singlem_dev_run_condense:
    input:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        done="{bench_dir}/output_singlem/singlem/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_dev/singlem_dev/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_dev/logs/singlem_dev/{sample}.log"
    shell:
        "pixi run -e singlem-dev " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --metapackage {input.db} &> {log}"


###############################################################################################
###############################################################################################
###############################################################################################
######### sylph

rule sylph_run:
    input:
        r1 = "{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2 = "{bench_dir}/local_reads/{sample}.2.fq.gz",
        db = sylph_package
    output:
        report="{bench_dir}/output_sylph/sylph/{sample}.tsv",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.done")
    threads: 8
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
        gtdb_bac_tax = "bac120_taxonomy_r207.tsv",
        gtdb_ar_tax = "ar53_taxonomy_r207.tsv",
    output:
        profile = "{bench_dir}/output_sylph/sylph/{sample}.profile",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.profile.done")
    shell:
        "pixi run -e singlem " \
        "python3 bin/sylph_to_condensed.py --sylph-genome {input.report} " \
        "--sample {wildcards.sample} " \
        "--bac-tax {input.gtdb_bac_tax} " \
        "--arc-tax {input.gtdb_ar_tax} > {output.profile}"
