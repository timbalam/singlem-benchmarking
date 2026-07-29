from os.path import join, dirname

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



singlem_metapackage = "../tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"
sylph_package = "../tool_reference_data/gtdb_database.syldb"

#####################################################################
      
rule bench7:
    input:
        expand("sra_mostly_novel/output_{tool}/{tool}/{sample}.profile",
               sample = datasets_bench7,
               tool = ['singlem', 'sylph'])

rule download_bench7:
    input:
        [f'sra_mostly_novel/local_reads/{sample}.1.fq.gz' for sample in datasets_bench7],
        [f'sra_mostly_novel/local_reads/{sample}.2.fq.gz' for sample in datasets_bench7]

rule bench8_baseline:
    input:
        expand("tuning/output_singlem/opal/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5)),
        expand("tuning/output_singlem_dev/opal/tune_s0.0g0.0f0.0o0.0c0.0p0.0d0.0r0.0mask{mask}/{sample}.opal_report",
               sample = datasets_bench8,
               mask = range(5)),
        expand("tuning/output_singlem_dev/singlem_dev/tune_s0.0g0.0f0.0o0.0c0.0p0.0d0.0r0.0mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5))
    
rule bench8_r_test:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8[:1],
               mask = range(1), tunestr = tune_bench8_r[:1]),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8[:1],
               mask = range(1), tunestr = tune_bench8_r[:1])

rule bench8_r:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5), tunestr = tune_bench8_r),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5), tunestr = tune_bench8_r)

rule bench8_d:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5), tunestr = tune_bench8_d),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5), tunestr = tune_bench8_d)

rule bench8_g:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5), tunestr = tune_bench8_g),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5), tunestr = tune_bench8_g)

rule bench8_s:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5), tunestr = tune_bench8_s),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5), tunestr = tune_bench8_s)

rule bench8_a:
    input:
        expand("tuning/output_singlem_dev/opal/tune_{tunestr}mask{mask}/{sample}.opal_report",
               sample = datasets_bench8, 
               mask = range(5), tunestr = tune_bench8_a),
        expand("tuning/output_singlem_dev/singlem_dev/tune_{tunestr}mask{mask}/{sample}.loss.tsv",
               sample = datasets_bench8,
               mask = range(5), tunestr = tune_bench8_a)

rule renew_singlem_dev_bench8:
    input:
        [f'tuning/output_singlem_dev/singlem_dev/{sample}.sma' for sample in datasets_bench8]

rule generate_communities_bench8:
    input:
        [f'tuning/truths/{sample}.finished' for sample in datasets_bench8],
        [f'tuning/local_reads/{sample}.finished' for sample in datasets_bench8],
        [f'tuning/truths/{sample}.condensed.biobox' for sample in datasets_bench8],
    output:
        done=touch("tuning/generate_communities.done")

rule generate_communities_bench8_test:
    input:
        'tuning/truths/marine0novelty0.5.finished',
        'tuning/local_reads/marine0novelty0.5.finished',
        'tuning/truths/marine0novelty0.5.condensed.biobox'

rule generate_community_and_reads_bench8:
    input:
        gtdb_bac_metadata = '../bac120_metadata_r207.tsv',
        gtdb_ar_metadata = '../ar53_metadata_r207.tsv',
        known_genome_list = '../1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '../4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '../4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.1.fq.gz",
        r2="tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.2.fq.gz",
        condensed = "tuning/truths/marine{coverage_number}novelty{novelty_ratio}.condensed",
        done = touch("tuning/truths/marine{coverage_number}novelty{novelty_ratio}.finished"),
        done2 = touch("tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.finished"),
    log: "tuning/local_reads/marine{coverage_number}novelty{novelty_ratio}.log"
    threads: 8
    shell:
        "mkdir -p tuning/truths tuning/local_reads && " \
        "pixi run -e art " \
        "python3 tuning/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file tuning/coverage_definitions/coverage{wildcards.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} -2 {output.r2} " \
        "--novelty-ratio {wildcards.novelty_ratio} " \
        "2> {log}"

ranks_percent_id_by_rank = ["genus", "family", "class", "order"]
datasets_percent_id_by_rank = [f'sample{i}' for i in range(2)]
reps_percent_id_by_rank = range(4)

rule renew_singlem_dev_percent_id_by_rank:
    input:
        [f'percent_id_by_rank/output_singlem_dev/singlem_dev/{rank}/{sample}-{rep}.sma'
         for sample in datasets_percent_id_by_rank
         for rank in ranks_percent_id_by_rank
         for rep in reps_percent_id_by_rank]

rule generate_communities_percent_id_by_rank:
    input:
        [f'percent_id_by_rank/truths/{rank}/{sample}-{rep}.finished'
         for sample in datasets_percent_id_by_rank
         for rank in ranks_percent_id_by_rank
         for rep in reps_percent_id_by_rank],
        [f'percent_id_by_rank/local_reads/{rank}/{sample}-{rep}.finished'
         for sample in datasets_percent_id_by_rank
         for rank in ranks_percent_id_by_rank
         for rep in reps_percent_id_by_rank],
        [f'percent_id_by_rank/truths/{rank}/{sample}-{rep}.condensed.biobox'
         for sample in datasets_percent_id_by_rank
         for rank in ranks_percent_id_by_rank
         for rep in reps_percent_id_by_rank],
    output:
        done=touch("percent_id_by_rank/generate_communities.done")

# percent-known-at in rank order sgfocpd 
def percent_known_at_by_rank(wildcards, rank_novelty = 1):
    idx = "gfocpd".index(wildcards.rank[0])
    known_at = ["0"] * 7
    known_at[0] = "1" # species
    known_at[idx+1] = str(rank_novelty)
    return " ".join(known_at)

rule generate_community_and_reads_percent_id_by_rank:
    input:
        gtdb_bac_metadata = '../bac120_metadata_r207.tsv',
        gtdb_ar_metadata = '../ar53_metadata_r207.tsv',
        known_genome_list = '../1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '../4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '../4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="percent_id_by_rank/local_reads/{rank}/sample{coverage_number}-{rep}.1.fq.gz",
        r2="percent_id_by_rank/local_reads/{rank}/sample{coverage_number}-{rep}.2.fq.gz",
        condensed = "percent_id_by_rank/truths/{rank}/sample{coverage_number}-{rep}.condensed",
        done = touch("percent_id_by_rank/truths/{rank}/sample{coverage_number}-{rep}.finished"),
        done2 = touch("percent_id_by_rank/local_reads/{rank}/sample{coverage_number}-{rep}.finished"),
    log: "percent_id_by_rank/local_reads/{rank}/sample{coverage_number}-{rep}.log"
    threads: 8
    params:
        percent_known_at = percent_known_at_by_rank
    shell:
        "mkdir -p percent_id_by_rank/truths/{wildcards.rank} percent_id_by_rank/local_reads/{wildcards.rank} && " \
        "pixi run -e art " \
        "python3 ../11_rankwise_novelty/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file tuning/coverage_definitions/coverage{wildcards.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} -2 {output.r2} " \
        "--percent-known-at {params.percent_known_at} " \
        "2> {log}"

module run_benchmarks:
    snakefile: "../run_benchmarks.smk"
    #prefix: "experiments/"

use rule truth_condensed_to_biobox from run_benchmarks

use rule tool_condensed_to_biobox from run_benchmarks with:
    input:
        profile = "{bench_dir}/output_{tool}/{tool}/{tunedir}{sample}.profile",
        truth = "{bench_dir}/truths/{sample}.condensed.biobox",
    output:
        biobox = "{bench_dir}/output_{tool}/biobox/{tunedir}{sample}.biobox"
    wildcard_constraints:
        tunedir="([^/]+/)?",
        sample="[^/]+"

use rule opal from run_benchmarks with:
    input:
        biobox = "{bench_dir}/{tool_output}/biobox/{tunedir}{sample}.biobox"
    params:
        output_dir = "{bench_dir}/{tool_output}",
        output_opal_dir = "{bench_dir}/{tool_output}/opal/{tunedir}{sample}.opal_output_directory",
        truth = "{bench_dir}/truths/{sample}.condensed.biobox",
    output:
        report="{bench_dir}/{tool_output}/opal/{tunedir}{sample}.opal_report",
        done=touch("{bench_dir}/{tool_output}/opal/{tunedir}{sample}.opal_report.done")
    wildcard_constraints:
        tunedir="([^/]+/)?",
        sample="[^/]+"

use rule download_fastq from run_benchmarks with:
    wildcard_constraints:
        bench_dir="^sra_mostly_novel"



###############################################################################################
###############################################################################################
###############################################################################################
#########
######### tool-specific rules - singlem first

use rule singlem_run_pipe from run_benchmarks

use rule singlem_run_condense from run_benchmarks

use rule singlem_dev_run_renew from run_benchmarks

use rule singlem_dev_run_condense from run_benchmarks

rule singlem_dev_mask_5fold_mask:
    input:
        report="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma",
        done="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma.done"
    output:
        ["{bench_dir}/output_singlem_dev/singlem_dev/{sample}.mask"+str(mask)+".txt" for mask in range(5)],
        done=touch("{bench_dir}/output_singlem_dev/singlem_dev/{sample}.mask.done")
    log:
        "{bench_dir}/output_singlem_dev/logs/singlem_dev/{sample}.mask.log"
    params:
        output_dir=lambda wildcards, input, output: dirname(input.report)
    shell:
        "pixi run -e singlem-dev " \
        "python3 bin/generate_masks.py --input-archive-otu-table {input.report} " \
        "--fold 5 --output-mask-dir {params.output_dir} &> {log}"

rule singlem_dev_run_condense_tune:
    input:
        report="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma",
        mask="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.mask{mask}.txt",
        done="{bench_dir}/output_singlem_dev/singlem_dev/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_dev/singlem_dev/tune_s{ts}g{tg}f{tf}o{to}c{tc}p{tp}d{td}r{tr}mask{mask}/{sample}.profile",
        loss="{bench_dir}/output_singlem_dev/singlem_dev/tune_s{ts}g{tg}f{tf}o{to}c{tc}p{tp}d{td}r{tr}mask{mask}/{sample}.loss.tsv",
        done=touch("{bench_dir}/output_singlem_dev/singlem_dev/tune_s{ts}g{tg}f{tf}o{to}c{tc}p{tp}d{td}r{tr}mask{mask}/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_dev/logs/singlem_dev/tune_s{ts}g{tg}f{tf}o{to}c{tc}p{tp}d{td}r{tr}mask{mask}/{sample}.log"
    params:
        output_dir=lambda wildcards, input, output: dirname(output.profile)
    shell:
        "mkdir -p {params.output_dir} && " \
        "pixi run -e singlem-dev " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --apply-nonneg-matrix-factorisation " \
        "--rank-penalty-steps {wildcards.ts} {wildcards.tg} {wildcards.tf} {wildcards.to} " \
        "{wildcards.tc} {wildcards.tp} {wildcards.td} {wildcards.tr} " \
        "--mask-otus-file {input.mask} " \
        "--max-num-steps 5000 "
        "--output-loss {output.loss} " \
        "--metapackage {input.db} &> {log}"

###############################################################################################
###############################################################################################
###############################################################################################
######### sylph

use rule sylph_sketch from run_benchmarks

use rule sylph_profile from run_benchmarks

use rule sylph_report_to_condensed from run_benchmarks
