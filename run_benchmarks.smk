from os.path import join, dirname, abspath
import polars as pl

datasets_bench5 = [f'marine{i}' for i in range(1)]

singlem_metapackage = "tool_reference_data/S4.1.0.GTDB_r207.metapackage_20240502.smpkg"
sylph_package = "tool_reference_data/gtdb_database.syldb"

#####################################################################

# ## bench 5 novelty

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

rule download_shakya:
    input:
        [f'7_shakya_synthetic/local_reads/{sample}.1.fq.gz' for sample in ['SRR606249', 'SRR606245']],
        [f'7_shakya_synthetic/local_reads/{sample}.2.fq.gz' for sample in ['SRR606249', 'SRR606245']],

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

# ## 6 zymo

def parse_samples(report_file):
    return pl.read_csv(report_file, separator = "\t", columns = "run_accession").to_series(0)

zymo_samples = parse_samples("6_zymo_synthetic/filereport_read_run_ERP121404_D6300.txt")

rule download_zymo:
    input:
        [f'6_zymo_synthetic/local_reads/{sample}.1.fq.gz' for sample in zymo_samples],
        [f'6_zymo_synthetic/local_reads/{sample}.2.fq.gz' for sample in zymo_samples],

rule download_zymo_reference_genomes:
    output:
        "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2.zip"
    log:
        "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2-download.log"
    shell:
        "bash -c " \
        "'mkdir -p 6_zymo_synthetic && " \
        "cd 6_zymo_synthetic && " \
        "wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip' &> {log}"

rule extract_zymo_reference_genomes:
    input:
        "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2.zip"
    output:
        directory("6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2")
    log:
        "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2-extract.log"
    shell:
        "bash -c " \
        "'cd 6_zymo_synthetic && " \
        "tar -xzf ZymoBIOMICS.STD.refseq.v2.zip' &> {log}"

# ## 8 cami2 strain madness

rule download_cami_strain:
    input:
        [f'8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/anonymous_reads.fq.gz' for sample_number in range(100)],
        "8_cami2_strain/strmgCAMI2_setup-extract.done",
        "8_cami2_strain/short_read/source_genomes/"

rule split_cami_strain:
    input:
        [f'8_cami2_strain/split_reads/strain{sample_number}.{dir}.fq.gz' for sample_number in range(100) for dir in [1,2]]

rule download_cami_strain_reads:
    output:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads.tar.gz"
    log:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads-download.log"
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain/short_read && " \
        "cd 8_cami2_strain/short_read && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule extract_cami_strain_reads:
    input:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads.tar.gz"
    output:
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/anonymous_reads.fq.gz"
    log:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads-extract.log"
    shell:
        "bash -c " \
        "'cd 8_cami2_strain/short_read && " \
        "tar -xzf strmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule split_cami_strain_reads:
    input:
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/anonymous_reads.fq.gz"
    output:
        r1="8_cami2_strain/split_reads/strain{sample_number}.1.fq.gz",
        r2="8_cami2_strain/split_reads/strain{sample_number}.2.fq.gz",
        done=touch("8_cami2_strain/split_reads/strain{sample_number}.done")
    log:
        "8_cami2_strain/split_reads/strain{sample_number}.log"
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain/split_reads && zcat {input[0]} | " \
        "bin/deinterleave_fastq.sh {output.r1} {output.r2} compress' &> {log}"

rule download_cami_strain_genomes:
    output:
        "8_cami2_strain/strmgCAMI2_genomes.tar.gz"
    log:
        "8_cami2_strain/strmgCAMI2_genomes-download.log"
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain && " \
        "cd 8_cami2_strain && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/strmgCAMI2_genomes.tar.gz' &> {log}"

rule extract_cami_strain_genomes:
    input:
        "8_cami2_strain/strmgCAMI2_genomes.tar.gz"
    output:
        directory("8_cami2_strain/short_read/source_genomes/")
    log:
        "8_cami2_strain/strmgCAMI2_genomes-extract.log"
    shell:
        "bash -c " \
        "'cd 8_cami2_strain && " \
        "tar -xzf strmgCAMI2_genomes.tar.gz' &> {log}"

rule download_cami_strain_setup:
    output:
        "8_cami2_strain/strmgCAMI2_setup.tar.gz",
    log:
        '8_cami2_strain/strmgCAMI2_setup-download.log'
    shell:
        "bash -c 'cd 8_cami2_strain && " \
        "rm -f strmgCAMI2_setup.tar.gz && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_setup.tar.gz' &> {log}"
 
rule extract_cami_strain_setup:
    input:
        "8_cami2_strain/strmgCAMI2_setup.tar.gz"
    output:
        touch("8_cami2_strain/strmgCAMI2_setup-extract.done")
    log:
        "8_cami2_strain/strmgCAMI2_setup-extract.log"
    shell:
        "bash -c " \
        "'cd 8_cami2_strain && " \
        "tar -xzf strmgCAMI2_setup.tar.gz' &> {log}"
       
rule download_cami2_profiles:
    output:
        '8_cami2_strain/taxonomic_profiling_cami2.tar.gz'
    log:
        '8_cami2_strain/cami2_profiles-download.log'
    shell:
        "bash -c 'cd 8_cami2_strain && " \
        "rm -f taxonomic_profiling_cami2.tar.gz && " \
        "wget https://zenodo.org/records/5006866/files/taxonomic_profiling_cami2.tar.gz?download=1 -O taxonomic_profiling_cami2.tar.gz' &> {log}"

rule extract_cami2_profiles:
    input:
        "8_cami2_strain/taxonomic_profiling_cami2.tar.gz"
    output:
        directory("8_cami2_strain/taxonomic_profiling_cami2/")
    log:
        "8_cami2_strain/taxonomic_profiling_cami2-extract.log"
    shell:
        "bash -c " \
        "'cd 8_cami2_strain && " \
        "tar -xzf taxonomic_profiling_cami2.tar.gz' &> {log}"

rule download_gtdbtk_r207_data:
    output:
        'tool_reference_data/gtdbtk_r207_v2_data.tar.gz'
    log:
        "tool_reference_data/gtdbtk_r207_data-download.log"
    shell:
        "bash -c "\
        "'cd tool_reference_data && "\
        "wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz' &> {log}"

rule extract_gtdbtk_r207_data:
    input:
        "tool_reference_data/gtdbtk_r207_v2_data.tar.gz"
    output:
        directory("tool_reference_data/release207_v2/")
    log:
        "tool_reference_data/gtdbtk_r207_v2_data-extract.log"
    shell:
        "bash -c " \
        "'cd tool_reference_data && " \
        "tar -xzf gtdbtk_r207_v2_data.tar.gz' &> {log}"

def gtdbtk_data_path(wildcards):
    return abspath("tool_reference_data/release207_v2")

rule gtdbtk_identify:
    input:
        fa="8_cami2_strain/short_read/source_genomes",
        data_path=gtdbtk_data_path
    output:
        output_dir=directory("8_cami2_strain/gtdbtk_r207/identify")
    log:
        "8_cami2_strain/gtdbtk_r207/identify.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        "pixi run -e gtdbtk " \
        "gtdbtk identify --genome_dir {input.fa} " \
        "--out_dir {output.output_dir} " \
        "--extension .fasta " \
        "&> {log}"

rule gtdbtk_align:
    input:
        id="8_cami2_strain/gtdbtk_r207/identify",
        data_path=gtdbtk_data_path
    output:
        output_dir=directory("8_cami2_strain/gtdbtk_r207/align")
    log:
        "8_cami2_strain/gtdbtk_r207/align.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        "pixi run -e gtdbtk " \
        "gtdbtk align --identify_dir {input.id} " \
        "--out_dir {output.output_dir} " \
        "&> {log}"

rule gtdbtk_classify:
    input:
        fa="8_cami2_strain/short_read/source_genomes",
        al="8_cami2_strain/gtdbtk_r207/align",
        data_path=gtdbtk_data_path
    output:
        output_dir=directory("8_cami2_strain/gtdbtk_r207/classify"),
        done=touch("8_cami2_strain/gtdbtk_r207/classify.done")
    log:
        "8_cami2_strain/gtdbtk_r207/classify.log"
    resources:
        mem_mb=64000
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        "pixi run -e gtdbtk " \
        "gtdbtk classify --genome_dir {input.fa} " \
        "--align_dir {input.al} " \
        "--extension .fasta " \
        #"--scratch_data " \
        "--out_dir {output.output_dir} " \
        "&> {log}"

rule gtdbtk_classify_wf:
    input:
        fa="8_cami2_strain/short_read/source_genomes",
        data_path=gtdbtk_data_path
    output:
        output_dir=directory("8_cami2_strain/gtdbtk_r207_wf")
    log:
        "8_cami2_strain/gtdb_r207_wf/classify_wf.log"
    shell:
        "GTDBTK_DATA_PATH={input.data_path} " \
        "pixi run -e gtdbtk " \
        "gtdbtk classify_wf --genome-dir {input.fa} " \
        "--out_dir {output.output_dir} " \
        "&> {log}"

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
tools_bench10 = ['singlem', 'sylph', 'singlem_dev', 'singlem_joint', 'singlem_inject',
                 'singlem_nnls', 'singlem_truecov']

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
        report="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma",
        done=touch("{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma.done")
    threads:
        8
    log:
        "{bench_dir}/output_singlem_joint/logs/singlem_joint/{sample}.log"
    wildcard_constraints:
        regime3="joint|inject"
    shell:
        "pixi run -e singlem-regime3 " \
        "singlem renew --threads {threads} --input-archive-otu-table {input.report} " \
        "--archive-otu-table {output.report} --metapackage {input.db} &> {log}"

rule singlem_joint_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_joint/singlem_joint/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_joint/logs/singlem_joint/{sample}.log"
    shell:
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --joint " \
        "--sylph-profile {input.sylph} " \
        "--metapackage {input.db} &> {log}"

rule singlem_nnls_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_nnls/singlem_nnls/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_nnls/singlem_nnls/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_nnls/logs/singlem_nnls/{sample}.log"
    shell:
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
        report="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax",
        done="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_truecov/singlem_truecov/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_truecov/singlem_truecov/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_truecov/logs/singlem_truecov/{sample}.log"
    shell:
        "pixi run -e singlem-regime3 " \
        "singlem condense --input-archive-otu-table {input.report} " \
        "-p {output.profile} --joint " \
        "--sylph-profile {input.sylph} " \
        "--alpha 1.0 " \
        "--metapackage {input.db} &> {log}"

rule singlem_inject_run_condense:
    input:
        report="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma",
        sylph="{bench_dir}/output_sylph/sylph/{sample}.tax.eff",
        done="{bench_dir}/output_singlem_joint/singlem_joint/{sample}.sma.done",
        db=singlem_metapackage
    output:
        profile="{bench_dir}/output_singlem_inject/singlem_inject/{sample}.profile",
        done=touch("{bench_dir}/output_singlem_inject/singlem_inject/{sample}.profile.done")
    log:
        "{bench_dir}/output_singlem_inject/logs/singlem_inject/{sample}.log"
    shell:
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
