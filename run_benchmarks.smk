from os.path import join, dirname, abspath
import polars as pl
import sys
import os
from ruamel.yaml import YAML

datasets_bench5 = [f'marine{i}' for i in range(1)]

singlem_metapackage = join(dirname(workflow.snakefile), "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg")
sylph_package = join(dirname(workflow.snakefile), "tool_reference_data/gtdb_database.syldb")

weebill_binary = abspath(join('weebill', 'target', 'release', 'weebill'))
#weebill_db = '/scratch/microbiome/woodcrob/non_sensitive/weebill_dbs/r207.100.syl2db'
weebill_db = join(dirname(workflow.snakefile), "tool_reference_data/weebill_r207.100.syl2db")
singlem_regime3_git_base_directory = join(dirname(workflow.snakefile), 'singlem_sylph_condense_regime')
regime3_weebill_scripts = join(
    singlem_regime3_git_base_directory, 'extras/singlem_weebill_benchmark/scripts')

singlem_regime3_metapackage_local = singlem_metapackage
regime3_weebill_mem_mb = 16000
singlem_ncbi_metapackage_local = join(dirname(workflow.snakefile), "tool_reference_data", "ncbi_names.GTDB_r232.metapackage_20261106.smpkg.zb.tar.gz")
cami_strain_samples = [f"strain{sample_number}" for sample_number in range(100)]

def parse_clinical_stool_samples(report_file):
    return pl.read_csv(report_file, separator = "\t", columns = "run").to_series(0)

clinical_stool_samples = parse_clinical_stool_samples(workflow.source_path("19_clinical_stool/PRJNA786578_annotate.tsv"))

metax_cami3_package = join(dirname(workflow.snakefile), "tool_reference_data/metax_cami3_db/metax_bavfph/metax_db.json")
metax_cami3_taxdmp = join(dirname(workflow.snakefile), "tool_reference_data/metax_cami3_taxdmp")

metax_gut_samples = [f"sample_{sample_number}" for sample_number in range(10)]

#####################################################################

# data simulation

rule bench5:
    input:
        expand("5_small_1novel/output_{tool}/opal/{sample}.opal_report",
               sample = datasets_bench5, tool = ['singlem', 'sylph', 'singlem_dev']),
        expand("5_small_1novel/output_{tool}/after_em/{sample}.sma",
               sample = datasets_bench5, tool = ['singlem', 'singlem_dev'])

rule generate_communities_bench5:
    input:
        [f'5_small_1novel/truths/{sample}.finished' for sample in datasets_bench5],
        [f'5_small_1novel/local_reads/{sample}.finished' for sample in datasets_bench5],
        [f'5_small_1novel/truths/{sample}.condensed.biobox' for sample in datasets_bench5],
    output:
        done=touch("5_small_1novel/generate_communities.done")

rule generate_community_and_reads_bench5:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="5_small_1novel/local_reads/{sample}.1.fq.gz",
        r2="5_small_1novel/local_reads/{sample}.2.fq.gz",
        condensed = "5_small_1novel/truths/{sample}.condensed",
        #genomewise = "5_small_1novel/truths/{sample}.genomewise.csv",
        done = touch("5_small_1novel/truths/{sample}.finished"),
        done2 = touch("5_small_1novel/local_reads/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('marine', ''),
    log: "5_small_1novel/local_reads/{sample}.log"
    threads: 8
    shell:
        "mkdir -p 5_small_1novel/truths 5_small_1novel/local_reads && " \
        "pixi run -e art " \
        "python3 5_small_1novel/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 5_small_1novel/coverage_definitions/coverage{params.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 5_small_1novel/local_reads/{wildcards.sample}.1.fq.gz " \
        "-2 5_small_1novel/local_reads/{wildcards.sample}.2.fq.gz " \
        "2> {log}"

# ## bench 7
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

# ## 8 cami2 strain madness
rule bench_cami_strain:
    input:
        expand("8_cami2_strain/output_{tool}/opal/{sample}.opal_report",
               sample = cami_strain_samples, tool = ['singlem', 'sylph', 'singlem_joint_a846591', 'singlem-regime3'])

# ## bench 9 related

datasets_bench9 = [f"sample{i}" for i in range(4, 6)]
percent_dominance_bench9 = [10, 50]

rule bench9:
    input:
        expand("9_related/output_{tool}/opal/dominance{percent_dom}/{sample}.opal_report",
               sample = datasets_bench9,
               percent_dom = percent_dominance_bench9,
               tool = ['singlem',
                       'sylph',
                       'singlem_dev',
                       'singlem_joint', 'singlem_inject',
                       'singlem_nnls', 'singlem_truecov'
                       ])

rule generate_communities_bench9:
    input:
        [f'9_related/truths/dominance{percent_dom}/{sample}.finished'
         for sample in datasets_bench9
         for percent_dom in percent_dominance_bench9],
        [f'9_related/local_reads/dominance{percent_dom}/{sample}.finished'
         for sample in datasets_bench9
         for percent_dom in percent_dominance_bench9],
        [f'9_related/truths/dominance{percent_dom}/{sample}.condensed.biobox'
         for sample in datasets_bench9
         for percent_dom in percent_dominance_bench9],
    output:
        done=touch("9_related/generate_communities.done")

rule generate_community_and_reads_bench9:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="9_related/local_reads/dominance{percent_dom}/{sample}.1.fq.gz",
        r2="9_related/local_reads/dominance{percent_dom}/{sample}.2.fq.gz",
        condensed = "9_related/truths/dominance{percent_dom}/{sample}.condensed",
        done = touch("9_related/truths/dominance{percent_dom}/{sample}.finished"),
        done2 = touch("9_related/local_reads/dominance{percent_dom}/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('sample', ''),
    log: "9_related/local_reads/dominance{percent_dom}/{sample}.log"
    threads: 8
    shell:
        "mkdir -p 9_related/truths 9_related/local_reads && " \
        "pixi run -e art " \
        "python3 9_related/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 9_related/coverage_definitions/coverage{params.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--percent-known 100 " \
        "--percent-dominant {wildcards.percent_dom} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} " \
        "-2 {output.r2} " \
        "2> {log}"

# ## bench 10 related and novel

datasets_bench10 = [f"sample{i}" for i in range(4, 6)]
percent_dominance_bench10 = [10, 50]
percent_known_bench10 = [0, 10, 50, 70, 100]
tools_bench10 = ['singlem', 'sylph', 'singlem_dev', 'singlem_joint_ee20bc', 'singlem_inject_ee20bc',
                 'singlem_nnls_ee20bc', 'singlem_truecov_ee20bc']

rule bench10:
    input:
        expand("10_related_and_novel/output_{tool}/opal/known{percent_known}/dominance{percent_dom}/{sample}.opal_report",
               sample = datasets_bench10,
               percent_dom = percent_dominance_bench10,
               percent_known = percent_known_bench10[1:],
               tool = tools_bench10),
        expand("10_related_and_novel/output_{tool}/opal/known{percent_known}/dominance{percent_dom}/{sample}.opal_report",
               sample = datasets_bench10,
               percent_dom = percent_dominance_bench10,
               percent_known = percent_known_bench10[:1],
               tool = tools_bench10[:1]+tools_bench10[2:]), # skip sylph for known0

rule generate_communities_bench10:
    input:
        [f'10_related_and_novel/truths/known{percent_known}/dominance{percent_dom}/{sample}.finished'
         for sample in datasets_bench10
         for percent_dom in percent_dominance_bench10
         for percent_known in percent_known_bench10],
        [f'10_related_and_novel/local_reads/known{percent_known}/dominance{percent_dom}/{sample}.finished'
         for sample in datasets_bench10
         for percent_dom in percent_dominance_bench10
         for percent_known in percent_known_bench10],
        [f'10_related_and_novel/truths/known{percent_known}/dominance{percent_dom}/{sample}.condensed.biobox'
         for sample in datasets_bench10
         for percent_dom in percent_dominance_bench10
         for percent_known in percent_known_bench10],
    output:
        done=touch("10_related_and_novel/generate_communities.done")

rule generate_community_and_reads_bench10:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="10_related_and_novel/local_reads/known{percent_known}/dominance{percent_dom}/{sample}.1.fq.gz",
        r2="10_related_and_novel/local_reads/known{percent_known}/dominance{percent_dom}/{sample}.2.fq.gz",
        condensed = "10_related_and_novel/truths/known{percent_known}/dominance{percent_dom}/{sample}.condensed",
        done = touch("10_related_and_novel/truths/known{percent_known}/dominance{percent_dom}/{sample}.finished"),
        done2 = touch("10_related_and_novel/local_reads/known{percent_known}/dominance{percent_dom}/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('sample', ''),
    log: "10_related_and_novel/local_reads/known{percent_known}/dominance{percent_dom}/{sample}.log"
    threads: 8
    shell:
        "mkdir -p 10_related_and_novel/truths 10_related_and_novel/local_reads && " \
        "pixi run -e art " \
        "python3 9_related/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 10_related_and_novel/coverage_definitions/coverage{params.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--percent-known {wildcards.percent_known} " \
        "--percent-dominant {wildcards.percent_dom} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} " \
        "-2 {output.r2} " \
        "2> {log}"

# ## bench 11 rankwise novelty

datasets_bench11 = [f"sample{i}" for i in [0, 4, 5]]
skew_bench11 = ["low_skew", "high_skew"]

rule bench11:
    input:
        expand("11_rankwise_novelty/output_{tool}/opal/{skew}/{sample}.opal_report",
               sample = datasets_bench11,
               skew = skew_bench11,
               tool = tools_bench10)

rule generate_communities_bench11:
    input:
        [f'11_rankwise_novelty/truths/{skew}/{sample}.finished'
         for skew in skew_bench11
         for sample in datasets_bench11],
        [f'11_rankwise_novelty/local_reads/{skew}/{sample}.finished'
         for sample in datasets_bench10
         for skew in skew_bench11],
        [f'11_rankwise_novelty/truths/{skew}/{sample}.condensed.biobox'
         for sample in datasets_bench10
         for skew in skew_bench11],
    output:
        done=touch("11_rankwise_novelty/generate_communities.done")

# percent-known-at in rank order sgfcopd for low and high skew novelty
skew_known_at = {
    "low_skew": "5 4 3 2 0 0 0",
    "high_skew": "100 10 5 2 1 0 0"
}

rule generate_community_and_reads_bench11:
    input:
        gtdb_bac_metadata = 'bac120_metadata_r207.tsv',
        gtdb_ar_metadata = 'ar53_metadata_r207.tsv',
        known_genome_list = '1_novel_strains/shadow_genome_paths.csv',
        novel_genomes_gtdbtk_output_directory = '4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207',
        novel_genome_list = '4_complex_and_novel/gtdbtk_batchfile.random1000.csv',
    output:
        r1="11_rankwise_novelty/local_reads/{skew}/sample{coverage_number}.1.fq.gz",
        r2="11_rankwise_novelty/local_reads/{skew}/sample{coverage_number}.2.fq.gz",
        condensed = "11_rankwise_novelty/truths/{skew}/sample{coverage_number}.condensed",
        done = touch("11_rankwise_novelty/truths/{skew}/sample{coverage_number}.finished"),
        done2 = touch("11_rankwise_novelty/local_reads/{skew}/sample{coverage_number}.finished"),
    log: "11_rankwise_novelty/local_reads/{skew}/sample{coverage_number}.log"
    threads: 8
    params:
        percent_known_at = lambda wildcards: skew_known_at[wildcards.skew]
    shell:
        "mkdir -p 11_rankwise_novelty/truths/{wildcards.skew} 11_rankwise_novelty/local_reads/{wildcards.skew} && " \
        "pixi run -e art " \
        "python3 11_rankwise_novelty/generate_community.py --art art_illumina --threads {threads} " \
        "--coverage-file 10_related_and_novel/coverage_definitions/coverage{wildcards.coverage_number}.tsv " \
        "--gtdb-bac-metadata {input.gtdb_bac_metadata} " \
        "--gtdb-ar-metadata {input.gtdb_ar_metadata} " \
        "--known-genome-list {input.known_genome_list} " \
        "--novel-genome-gtdbtk-output {input.novel_genomes_gtdbtk_output_directory} " \
        "--novel-genome-list {input.novel_genome_list} " \
        "--output-condensed {output.condensed} " \
        "-1 {output.r1} -2 {output.r2} " \
        "--percent-known-at {params.percent_known_at} " \
        "2> {log}"

# ## bench18 metax gut

# ## bench19 clinical stool
rule bench19:
    input:
        expand("19_clinical_stool/output_{tool}/{tool}/{sample}.profile",
               tool = ["singlem-regime3", "sylph", "metax_cami3"],
               sample = ['SRR17233267'])
               #sample = clinical_stool_samples)


rule bench19_test:
    input:
        expand("19_clinical_stool/output_{tool}/{tool}/{sample}.profile",
               tool = ["singlem-regime3", "sylph", "metax_cami3"][2:],
               sample = ['SRR17298983'])

# ## bench20 metax gut ncbi
rule bench20:
    input:
        expand("20_metax_gut_ncbi/output_{tool}/opal/{sample}.opal_report",
               tool = ["metax-ncbi", "singlem-ncbi"][1:],
               sample = metax_gut_samples[:1])

# utils

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

rule truth_condensed_to_biobox:
    input:
        condensed = "{bench_dir}/truths/{sample}.condensed",
    output:
        biobox = "{bench_dir}/truths/{sample}.condensed.biobox"
    shell:
        "pixi run -e singlem " \
        "python3 " + workflow.source_path("bin/condensed_profile_to_biobox.py") + " " \
        "--input-condensed-table {input.condensed} " \
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
        after_em="{bench_dir}/output_singlem/after_em/{sample}.sma",
        done=touch("{bench_dir}/output_singlem/singlem/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem/logs/singlem/{sample}.log"
    shell:
        "pixi run -e singlem " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "--output-after-em-otu-table {output.after_em} " \
        "-p {output.profile} --metapackage {input.db} &> {log}"

READS_SEP = {
    '18_metax_gut': '_',
    '20_metax_gut_ncbi': '_'
}
def singlem_regime3_input(wildcards):
    reads_path = READS_PATH.get(wildcards.bench_dir, 'local_reads')
    reads_sep = READS_SEP.get(wildcards.bench_dir, '.')
    reads_pref = f"{wildcards.bench_dir}/{reads_path}/{wildcards.sample}{reads_sep}"
    return {
        'r1': f"{reads_pref}1.fq.gz",
        'r2': f"{reads_pref}2.fq.gz",
        'db': singlem_regime3_metapackage_local
    }
rule singlem_regime3_pipe_to_archive:
    input: unpack(singlem_regime3_input)
    output:
        archive="{bench_dir}/output_<tool>/<tool>/{sample}.archive.json",
        done=touch("{bench_dir}/output_<tool>/<tool>/{sample}.archive.done")
    threads: 8
    pathvars:
        tool = 'singlem-regime3'
    resources:
        runtime=180,  # singlem pipe --no-weebill; 3h ceiling
    log:
        "{bench_dir}/output_<tool>/logs/<tool>/{sample}.pipe.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "singlem pipe --threads {threads} -1 {input.r1} -2 {input.r2} "
        "--no-weebill --archive-otu-table {output.archive} --metapackage {input.db} &> {log}"

def singlem_ncbi_input(wildcards):
    inp = singlem_regime3_input(wildcards)
    inp["db"] = singlem_ncbi_metapackage_local
    return inp

use rule singlem_regime3_pipe_to_archive as singlem_ncbi_pipe_to_archive with:
    input: unpack(singlem_ncbi_input)
    pathvars:
        tool = 'singlem_ncbi'

rule singlem_regime3_weebill_profile:
    input:
        r1="{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2="{bench_dir}/local_reads/{sample}.2.fq.gz",
        database=weebill_db
    output:
        profile="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.tsv",
        done=touch("{bench_dir}/output_<tool>/<tool>/{sample}.weebill.done"),
    pathvars:
        tool = "singlem-regime3"
    threads: 8
    resources:
        mem_mb=regime3_weebill_mem_mb,
        runtime=120,  # weebill profile is fast; 2h is ample
    log:
        "{bench_dir}/output_<tool>/logs/<tool>/{sample}.weebill.log",
    shell:
        "{weebill_binary} "
        "profile --two-stage -t {threads} -c 100 "
        "-1 {input.r1} -2 {input.r2} -o {output.profile} {input.database} -u &> {log} && "
        # weebill writes only a header when nothing is detected; ensure >=1 data row
        "awk 'END {{ exit NR < 2 }}' {output.profile}"

use rule singlem_regime3_weebill_profile as singlem_ncbi_weebill_profile with:
    input:
        r1 = "{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2 = "{bench_dir}/local_reads/{sample}.2.fq.gz",
        database = weebill_db
    pathvars:
        tool = 'singlem_ncbi'

rule singlem_regime3_annotate_weebill:
    input:
        profile="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.tsv",
        db=singlem_regime3_metapackage_local,
    output:
        annotated="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.annotated.tsv",
    pathvars:
        tool = "singlem-regime3"
    resources:
        runtime=60,  # trivial annotation step
    log:
        "{bench_dir}/output_<tool>/logs/<tool>/{sample}.annotate.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "PYTHONPATH={singlem_regime3_git_base_directory} python3 {regime3_weebill_scripts}/annotate_weebill.py "
        "--profile {input.profile} --metapackage {input.db} --output {output.annotated} &> {log}"

use rule singlem_regime3_annotate_weebill as singlem_ncbi_annotate_weebill with:
    input:
        profile="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.tsv",
        db = singlem_ncbi_metapackage_local
    pathvars:
        tool = 'singlem_ncbi'

rule singlem_regime3_condense_joint:
    input:
        archive="{bench_dir}/output_<tool>/<tool>/{sample}.archive.json",
        weebill_profile="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.annotated.tsv",
        db=singlem_regime3_metapackage_local,
    output:
        profile="{bench_dir}/output_<tool>/<tool>/{sample}.profile",
        done=touch("{bench_dir}/output_<tool>/<tool>/{sample}.profile.done")
    pathvars:
        tool = "singlem-regime3"
    threads: 1
    resources:
        runtime=60,  # condense is quick
    log:
        "{bench_dir}/output_singlem-regime3/logs/singlem-regime3/{sample}.condense.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "singlem condense --input-archive-otu-table {input.archive} --metapackage {input.db} "
        "--weebill-profile {input.weebill_profile} --joint --taxonomic-profile {output.profile} "
        "--joint-pin-weebill-species --joint-novel-budget --alpha 1 "
        "&> {log}"

use rule singlem_regime3_condense_joint as singlem_ncbi_condense_joint with:
    input:
        archive="{bench_dir}/output_<tool>/<tool>/{sample}.archive.json",
        weebill_profile="{bench_dir}/output_<tool>/<tool>/{sample}.weebill.annotated.tsv",
        db=singlem_ncbi_metapackage_local
    pathvars:
        tool = 'singlem_ncbi'

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
    shell:
        "pixi run -e singlem-dev " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --apply-nonneg-matrix-factorisation " \
        "--output-after-em-otu-table {output.after_em} " \
        "--metapackage {input.db} &> {log}"

rule singlem_joint_run_renew:
    input:
        report="{bench_dir}/output_singlem/singlem/{sample}.sma",
        db=singlem_metapackage
    output:
        report="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma",
        done=touch("{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma.done")
    threads:
        8
    log:
        "{bench_dir}/output_singlem_joint_{hash}/logs/singlem_joint_{hash}/{sample}.log"
    shell:
        "bash -c 'cd singlem_sylph_condense_regime && " \
        "git rev-parse HEAD | grep -q ^{wildcards.hash}' && " \
        "pixi run -e singlem-regime3 " \
        "singlem renew --threads {threads} --input-archive-otu-table {input.report} " \
        "--archive-otu-table {output.report} --metapackage {input.db} &> {log}"

rule singlem_joint_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_joint_{hash}/logs/singlem_joint_{hash}/{sample}.log"
    shell:
        "bash -c 'cd singlem_sylph_condense_regime && " \
        "git rev-parse HEAD | grep -q ^{wildcards.hash}' && " \
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --joint " \
        "--sylph-profile {input.sylph} " \
        "--metapackage {input.db} &> {log}"

rule singlem_nnls_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_nnls_{hash}/singlem_nnls_{hash}/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_nnls_{hash}/singlem_nnls_{hash}/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_nnls_{hash}/logs/singlem_nnls_{hash}/{sample}.log"
    shell:
        "bash -c 'cd singlem_sylph_condense_regime && " \
        "git rev-parse HEAD | grep -q ^{wildcards.hash}' && " \
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --joint " \
        "--sylph-profile {input.sylph} " \
        "--joint-sylph-weight 0 --joint-absence-weight 0 " \
        "--metapackage {input.db} &> {log}"

rule test_bench9:
    input:
        "9_related/output_sylph/sylph/dominance50/sample5.tax"

rule singlem_truecov_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax",
        done="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_truecov_{hash}/singlem_truecov_{hash}/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_truecov_{hash}/singlem_truecov_{hash}/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_truecov_{hash}/logs/singlem_truecov_{hash}/{sample}.log"
    shell:
        "bash -c 'cd singlem_sylph_condense_regime && " \
        "git rev-parse HEAD | grep -q ^{wildcards.hash}' && " \
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --joint " \
        "--sylph-profile {input.sylph} " \
        "--alpha 1.0 " \
        "--metapackage {input.db} &> {log}"

rule singlem_inject_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint_{hash}/singlem_joint_{hash}/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_inject_{hash}/singlem_inject_{hash}/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_inject_{hash}/singlem_inject_{hash}/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_inject_{hash}/logs/singlem_inject_{hash}/{sample}.log"
    shell:
        "bash -c 'cd singlem_sylph_condense_regime && " \
        "git rev-parse HEAD | grep -q ^{wildcards.hash}' && " \
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} " \
        "--sylph-profile {input.sylph} " \
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
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.paired.sylsp.done")
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
        report = "{bench_dir}/output_sylph/sylph/{sample}.tsv{eff}",
        gtdb_bac_tax = "bac120_taxonomy_r207.tsv",
        gtdb_ar_tax = "ar53_taxonomy_r207.tsv",
    wildcard_constraints:
        eff="(\.eff)?"
    output:
        profile = "{bench_dir}/output_sylph/sylph/{sample}.profile{eff}",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.profile{eff}.done")
    shell:
        "pixi run -e singlem " \
        "python3 bin/sylph_to_condensed.py --sylph-genome {input.report} " \
        "--sample {wildcards.sample} " \
        "--bac-tax {input.gtdb_bac_tax} " \
        "--arc-tax {input.gtdb_ar_tax} > {output.profile}"

rule sylph_profile_eff:
    input:
        sp = "{bench_dir}/output_sylph/sylph/{sample}.paired.sylsp",
        db = sylph_package
    output:
        report="{bench_dir}/output_sylph/sylph/{sample}.tsv.eff",
        done=touch("{bench_dir}/output_sylph/sylph/{sample}.tsv.eff.done")
    threads: 8
    resources:
        mem_mb=32000
    log:
        "{bench_dir}/output_sylph/logs/sylph/{sample}.eff.log"
    shell:
        "pixi run -e sylph " \
        "sylph profile {input.db} {input.sp} -t {threads} " \
        "-o {output.report} 2> {log}"

rule sylph_tax:
    input:
        profile = "{bench_dir}/output_sylph/sylph/{sample}.profile{eff}"
    output:
        tax = "{bench_dir}/output_sylph/sylph/{sample}.tax{eff}"
    wildcard_constraints:
        eff="(\.eff)?"
    shell:
        "sed '1s/coverage/Eff_cov/' {input.profile} > {output.tax}"

###############################################################################################
###############################################################################################
###############################################################################################
######### metax

rule metax_cami3_profile:
    input:
        r1 = "{bench_dir}/local_reads/{sample}.1.fq.gz",
        r2 = "{bench_dir}/local_reads/{sample}.2.fq.gz",
        db = metax_cami3_package,
        dump_dir = metax_cami3_taxdmp
    output:
        report="{bench_dir}/output_metax_cami3/metax_cami3/{sample}.profile.txt",
        done=touch("{bench_dir}/output_metax_cami3/metax_cami3/{sample}.done")
    threads: 8
    resources:
        mem_mb=128000
    log:
        "{bench_dir}/output_metax_cami3/logs/metax_cami3/{sample}.profile.log"
    params:
        prefix = lambda wildcards: f"{wildcards.bench_dir}/output_metax_cami3/metax_cami3/{wildcards.sample}"
    shell:
        "pixi run -e metax "
        "metax profile --dmp-dir {input.dump_dir} "
        "--db {input.db} "
        "-i {input.r1},{input.r2} "
        "-p "
        "-o {params.prefix} "
        "-t {threads} 2> {log}"

rule metax_cami3_report_to_condensed:
    input:
        report = "{bench_dir}/output_metax_cami3/metax_cami3/{sample}.profile.txt",
        dump_dir = metax_cami3_taxdmp
    output:
        profile = "{bench_dir}/output_metax_cami3/metax_cami3/{sample}.profile",
        done=touch("{bench_dir}/output_metax_cami3/metax_cami3/{sample}.profile.done")
    shell:
        "echo -e 'coverage\ttaxonomy\tsample' > {output.profile} && "
        "cut -f2,5 {input.report} | "
        "pixi run -e taxonkit taxonkit reformat2 -I 1 -f 'd__{{domain|acellular root|superkingdom}};p__{{phylum}};c__{{class}};o__{{order}};f__{{family}};g__{{genus}};s__{{species}}' | "
        "sed 's/$/\t{wildcards.sample}/' | "
        "cut -f2- >> {output.profile}"
