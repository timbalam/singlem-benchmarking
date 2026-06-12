from os.path import join, dirname

datasets_bench5 = [f'marine{i}' for i in range(1)]

datasets_bench7 = ['SRR8648366', 'SRR29850984']

novelties_bench8 = [0.5, 0.95, 1.25]
datasets_bench8 = [f'marine{i}novelty{nr}' for i, nr in enumerate(novelties_bench8)]

tune_bench8_r = [
    f's{s}g{g}f{f}o{o}c{c}p{p}d{d}r{r}'
    for s in [0.0]
    for g in [0.0]
    for f in [0.0]
    for o in [0.0]
    for c in [0.0]
    for p in [0.0]
    for d in [0.0]
    for r in [0.0, 0.1, 0.5, 1]
]
tune_bench8_d = [
    f's{s}g{g}f{f}o{o}c{c}p{p}d{d}r{r}'
    for s in [0.3]
    for g in [0.3]
    for f in [0.3]
    for o in [0.3]
    for c in [0.3]
    for p in [0.3]
    for d in [0.3, 0.4, 0.5, 0.6]
    for r in [round(3.1 - d - p - c - o - f - g - s, 5)]
]
tune_bench8_g = [
    f's{s}g{g}f{f}o{o}c{c}p{p}d{d}r{r}'
    for s in [0.0]
    for g in [0.0, 0.01, 0.1, 1.0, 10.0]
    for f in [round(5000.0 - g - s, 4)]
    for o in [0.0]
    for c in [0.0]
    for p in [0.0]
    for d in [0.0]
    for r in [0.0]
]
tune_bench8_s = [
    f's{s}g{g}f{f}o{o}c{c}p{p}d{d}r{r}'
    for s in [0.0, 0.01, 0.1, 1.0, 10.0, 100.0, 500.0]
    for g in [round(5000.0 - s, 4)]
    for f in [0.0]
    for o in [0.0]
    for c in [0.0]
    for p in [0.0]
    for d in [0.0]
    for r in [0.0]
]
tune_bench8_a = [
    f's{s}g{g}f{f}o{o}c{c}p{p}d{d}r{r}'
    for s in [500.0, 1000.0, 2000.0, 5000.0, 10000.0]
    for g in [0.0]
    for f in [0.0]
    for o in [0.0]
    for c in [0.0]
    for p in [0.0]
    for d in [0.0]
    for r in [0.0]
]



singlem_metapackage = "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"
sylph_package = "tool_reference_data/gtdb_database.syldb"

#####################################################################

rule bench5:
    input:
        expand("5_novelty/output_{tool}/opal/{sample}.opal_report",
               sample = datasets_bench5, tool = ['singlem', 'sylph', 'singlem_dev']),
        expand("5_novelty/output_{tool}/after_em/{sample}.sma",
               sample = datasets_bench5, tool = ['singlem', 'singlem_dev'])

rule generate_communities_bench5:
    input:
        [f'5_novelty/truths/{sample}.finished' for sample in datasets_bench5],
        [f'5_novelty/local_reads/{sample}.finished' for sample in datasets_bench5],
        [f'5_novelty/truths/{sample}.condensed.biobox' for sample in datasets_bench5],
    output:
        done=touch("5_novelty/generate_communities.done")

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

rule download_shakya:
    input:
        '7_shakya_synthetic/local_reads/SRR606249.1.fq.gz',
        '7_shakya_synthetic/local_reads/SRR606249.2.fq.gz',


rule generate_shakya_truth_condensed_format:
    input:
        sup_xlsx = "7_shakya_synthetic/emi12086-sup-0010-tables1.xlsx",
        gtdb_bac_tax = "bac120_taxonomy_r207.tsv",
        gtdb_ar_tax = "ar53_taxonomy_r207.tsv",
    output:
        condensed = "7_shakya_synthetic/truths/{sample}.condensed",
    log:
        "7_shakya_synthetic/logs/generate_truth_condensed_format-{sample}.log"
    shell:
        "pixi run -e singlem " \
        "python3 bin/shakya_sup_table_to_condensed.py --supplementary-xlsx {input.sup_xlsx} " \
        "--sample {wildcards.sample} " \
        "--bac-tax {input.gtdb_bac_tax} " \
        "--arc-tax {input.gtdb_ar_tax} > {output.profile} 2> {log}"


rule generate_communities_bench8:
    input:
        [f'8_tuning/truths/{sample}.finished' for sample in datasets_bench8],
        [f'8_tuning/local_reads/{sample}.finished' for sample in datasets_bench8],
        [f'8_tuning/truths/{sample}.condensed.biobox' for sample in datasets_bench8],
    output:
        done=touch("8_tuning/generate_communities.done")

rule generate_communities_bench8_test:
    input:
        '8_tuning/truths/marine0novelty0.5.finished',
        '8_tuning/local_reads/marine0novelty0.5.finished',
        '8_tuning/truths/marine0novelty0.5.condensed.biobox'

rule generate_community_and_reads_bench8:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="8_tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.1.fq.gz",
        r2="8_tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.2.fq.gz",
        condensed = "8_tuning/truths/marine{coverage_number}novelty{novelty_ratio}.condensed",
        done = touch("8_tuning/truths/marine{coverage_number}novelty{novelty_ratio}.finished"),
        done2 = touch("8_tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.finished"),
    log: "8_tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.log"
    threads: 8
    shell:
        "mkdir -p 8_tuning/truths 8_tuning/local_reads && " \
        "pixi run -e art " \
        "python3 8_tuning/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 8_tuning/coverage_definitions/coverage{wildcards.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} -2 {output.r2} " \
        "--novelty-ratio {wildcards.novelty_ratio} " \
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
    wildcard_constraints:
        sample="[^/]+"
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
    wildcard_constraints:
        sample="[^/]+"
    shell:
        "pixi run -e opal " \
        "opal.py -g {params.truth} -o {params.output_opal_dir} {input.biobox} || echo 'expected opal non-zero exit status'; mv {params.output_opal_dir}/results.tsv {output.report} && rm -rf {params.output_opal_dir}"

rule download_fastq:
    output:
        r1="{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2="{bench_dir}/local_reads/{sample}.2.fq.gz",
        done=touch("{bench_dir}/local_reads/{sample}.done")
    log:
        "{bench_dir}/local_reads/{sample}.log"
    threads:
        1
    wildcard_constraints:
        bench_dir="^7_.+"
    shell:
        "pixi run -e kingfisher " \
        "kingfisher get -r {wildcards.sample} " \
        "--output_directory {wildcards.bench_dir}/local_reads " \
        "-m ena-ftp prefetch -f fastq.gz " \
	    "&& mv {wildcards.bench_dir}/local_reads/{wildcards.sample}_1.fastq.gz {output.r1} " \
        "&& mv {wildcards.bench_dir}/local_reads/{wildcards.sample}_2.fastq.gz {output.r2} " \
        "&> {log}"



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
    wildcard_constraints:
        tunedir="([^/]+/)?",
        sample="[^/]+"
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
        after_em="{bench_dir}/output_singlem/after_em/{sample}.sma",
        done=touch("{bench_dir}/output_singlem/singlem/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem/logs/singlem/{sample}.log"
    shell:
        "pixi run -e singlem " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "--output-after-em-otu-table {output.after_em} " \
        "-p {output.profile} --metapackage {input.db} &> {log}"

rule singlem_dev_run_renew:
    input:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        db=singlem_metapackage
    output:
        report="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma",
        done=touch("{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma.done")
    threads:
        8
    log:
        "{bench_dir}/output_singlem_dev/logs/singlem_dev/{sample}.log"
    shell:
        "pixi run -e singlem-dev " \
        "singlem renew --threads {threads} --input-archive-otu-table {input.report} " \
        "--archive-otu-table {output.report} --metapackage {input.db} &> {log}"

rule singlem_dev_run_condense:
    input:
        report="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma",
        done="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.profile",
        after_em="{bench_dir}/output_singlem_dev/after_em/{sample}.sma",
        done=touch("{bench_dir}/output_singlem_dev/singlem_dev/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_dev/logs/singlem_dev/{sample}.log"
    wildcard_constraints:
        sample='[^/]+'
    shell:
        "pixi run -e singlem-dev " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --apply-nonneg-matrix-factorisation " \
        "--output-after-em-otu-table {output.after_em} " \
        "--metapackage {input.db} &> {log}"

###############################################################################################
###############################################################################################
###############################################################################################
######### sylph

def sylph_sketch_dir(wildcards, input, output):
    return dirname(output.sp)

rule sylph_sketch:
    input:
        r1 = "{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2 = "{bench_dir}/local_reads/{sample}.2.fq.gz"
    output:
        sp="{bench_dir}/output_sylph/sylph/{sample}.paired.sylsp",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.paired.slysp.done")
    threads: 8
    resources:
        mem_mb=32000
    log:
        "{bench_dir}/output_sylph/logs/sylph/{sample}.sketch.log"
    params:
        dir=sylph_sketch_dir
    shell:
        "pixi run -e sylph " \
        "sylph sketch -1 {input.r1} -2 {input.r2} -t {threads} " \
        "-S {wildcards.sample} -d {params.dir} " \
        "2> {log}"

rule sylph_profile:
    input:
        sp = "{bench_dir}/output_sylph/sylph/{sample}.paired.sylsp",
        db = sylph_package
    output:
        report="{bench_dir}/output_sylph/sylph/{sample}.tsv",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.done")
    threads: 8
    resources:
        mem_mb=32000
    log:
        "{bench_dir}/output_sylph/logs/sylph/{sample}.profile.log"
    shell:
        "pixi run -e sylph " \
        "sylph profile {input.db} {input.sp} -t {threads} " \
        "-u --read-seq-id 99.5 " \
        "-o {output.report} 2> {log}"

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
