import os.path
from os.path import join, abspath
import os
import polars as pl

output_directory = 'tool_reference_data'
# mkdir
os.makedirs(output_directory, exist_ok=True)

# https://zenodo.org/records/11107165/files/S4.1.0.GTDB_r207.metapackage_20240502.smpkg.zb.tar.gz?download=1
singlem_metapackage = join(output_directory, 'S4.1.0.GTDB_r207.metapackage_20240502.smpkg')
singlem_metapackage_tgz = singlem_metapackage + '.zb.tar.gz'

metaphlan_db = join(output_directory, 'metaphlan_bowtiedb')
metaphlan42_db = join(output_directory, 'metaphlan42_bowtiedb')
metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'

motus_db = join(output_directory, 'motus', 'db_mOTU')
motus_gtdb_tsv = join(output_directory, 'motus', 'mOTUs_3.0.0_GTDB_tax.tsv')

kraken_db = join(output_directory, "struo2-kraken-GTDB_release207")

# Must be at genome level otherwise sourmash tax annotate fails
sourmash_db_taxonomy = join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.sqldb')
sourmash_db_dna = join(output_directory, 'sourmash', 'gtdb-rs207.genomic-reps.dna.k31.zip')

kaiju_db_original_dir = join(output_directory, 'kaiju')
kaiju_db_original_progenomes_fmi = kaiju_db_original_dir + '/kaiju_db_progenomes.fmi'
kaiju_db_original_progenomes_nodes = kaiju_db_original_dir + '/nodes.dmp'
kaiju_db_original_progenomes_names = kaiju_db_original_dir + '/names.dmp'

map2b_checkout_dir = join(output_directory, 'MAP2B-checkout')
map2b_db = os.path.join(map2b_checkout_dir, 'database/GTDB')

# metabuli_db_dir = join(output_directory, 'metabuli')
# metabuli_db = join(metabuli_db_dir, 'gtdb')

weebill_syldb = os.path.join(output_directory, "weebill_r207.100.syldb")
weebill_db = os.path.join(output_directory, "weebill_r207.100.syl2db")
weebill_binary = 'weebill/target/release/weebill'

# map2b is excluded: its database is GTDB r202 (a release behind this benchmark's
# r207) and its DownloadDB.py figshare links currently return empty 202
# responses, so the download cannot complete. See the commented-out map2b rules
# below.s
## metabuli download is not scripted because it is via sharepoint, which gives an indirect link.
tools = ['singlem', 'metaphlan', 'motus', 'kraken', 'sourmash', 'kaiju', 'metaphlan42']

all_tools_input = [join(output_directory, f'{tool}.done') for tool in tools]

# Benchmarks 8 and 9's novel-in-r214 community members, fetched from NCBI rather
# than from the Zenodo shadow pool; see bench89_novel_genomes_download.
novel_r214_genome_accessions_file = 'novel_r214_genome_accessions.txt'
with open(novel_r214_genome_accessions_file) as f:
    novel_r214_genome_accessions = [line.strip() for line in f if line.strip()]
novel_r214_genome_fastas = expand("novel_r214_genomes/{accession}.fna",
                                  accession=novel_r214_genome_accessions)

cami2_marine_reads = [
    "3_cami2_marine/split_reads/marine{sample_number}.done".format(sample_number=sample_number)
    for sample_number in range(10)
]

def parse_zymo_report_samples(report_file):
    return pl.read_csv(report_file, separator = "\t", columns = "run_accession").to_series(0)

zymo_samples = parse_zymo_report_samples(workflow.source_path("6_zymo_synthetic/filereport_read_run_ERP121404_D6300.txt"))

zymo_reads = [
    f'6_zymo_synthetic/local_reads/{sample}.{dir}.fq.gz' for sample in zymo_samples
    for dir in [1, 2]
]

shakya_reads = [
    f'7_shakya_synthetic/local_reads/{sample}.{dir}.fq.gz'
    for sample in ['SRR606249', 'SRR606245']
    for dir in [1, 2]
]

cami_strain_reads = [
    f"8_cami2_strain/local_reads/strain{sample_number}.{dir}.fq.gz"
    for sample_number in range(100)
    for dir in [1,2]
]
cami_strain_coverage_definitions = [
    f"8_cami2_strain/coverage_definitions/strain{sample_number}.tsv"
    for sample_number in range(100)
]
metax_gut_reads = [
    f"18_metax_gut/split_reads/sample_{sample_number}.{dir}.fq.gz"
    for sample_number in range(10)
    for dir in [1, 2]
]

def parse_clinical_stool_samples(report_file):
    return pl.read_csv(report_file, separator = "\t", columns = "run").to_series(0)

clinical_stool_samples = parse_clinical_stool_samples(workflow.source_path("19_clinical_stool/PRJNA786578_annotate.tsv"))

clinical_stool_reads = [
    f"19_clinical_stool/local_reads/{sample}.{dir}.fq.gz"
    for sample in clinical_stool_samples
    for dir in [1, 2]
]

all_benchmarks_datasets = [
    join(output_directory, 'shadow-genomes.done'),
    cami2_marine_reads,
    "2_phylogenetic_novelty/genomes",
    "2_phylogenetic_novelty/genome_pairs",
    # benchmarks 8 and 9's novel-in-r214 members; see the rule at the bottom.
    novel_r214_genome_fastas,
    zymo_reads,
    "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2/Genomes",
    shakya_reads,
    cami_strain_reads,
    cami_strain_coverage_definitions,
    "8_cami2_strain/genomes",
    metax_gut_reads
]

rule all:
    input:
        all_tools_input,
        all_benchmarks_datasets

# tools
rule all_tools:
    input:
        all_tools_input

rule metaphlan:
    output:
        done=touch(join(output_directory, 'metaphlan.done')),
        metaphlan_db=directory(metaphlan_db)
    log:
        join(output_directory, 'metaphlan.log')
    localrule:
        True
    shell:
        'eval "$(pixi shell-hook -e metaphlan)" && '
        'metaphlan --install --bowtie2db {metaphlan_db} --index {metaphlan_index} &> {log}'

## metaphlan 4.2.2
rule metaphlan42:
    output:
        done=touch(join(output_directory, 'metaphlan42.done')),
        metaphlan_db=directory(metaphlan42_db)
    log:
        join(output_directory, 'metaphlan42.log')
    localrule:
        True
    shell:
        'pixi run --environment metaphlan42 '
        'metaphlan --install --db_dir {metaphlan42_db} --index {metaphlan_index} &> {log}'

rule kraken_download:
    output:
        done=touch(join(output_directory, 'kraken.done')),
    log:
        abspath(join(output_directory, 'kraken-download.log'))
    localrule:
        True
    shell:
        "mkdir -pv {kraken_db} && cd {kraken_db} && wget --no-directories -r -np -c --progress=dot:mega -e robots=off http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release207/kraken2/ &> {log}"

rule sourmash:
    input:
        join(output_directory, 'sourmash-dna.done'),
        join(output_directory, 'sourmash-tax.done'),
    output:
        done=touch(join(output_directory, 'sourmash.done')),
    localrule:
        True

rule sourmash_sigs:
    output:
        done=touch(join(output_directory, 'sourmash-dna.done')),
        sourmash_dna=sourmash_db_dna,
    log:
        join(output_directory, 'sourmash-sigs-download.log')
    localrule:
        True
    shell:
        'wget -O {output.sourmash_dna} https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip &> {log}'

rule sourmash_tax_download:
    output:
        done=touch(join(output_directory, 'sourmash-tax-download.done')),
        sourmash_taxonomy=join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.reps.csv.gz'),
    log:
        join(output_directory, 'sourmash-tax-download.log')
    localrule:
        True
    shell:
        'wget -O {output.sourmash_taxonomy} https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.reps.csv.gz &> {log}'

rule sourmash_tax_convert:
    input:
        done=join(output_directory, 'sourmash-tax-download.done'),
        sourmash_taxonomy=join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.reps.csv.gz'),
    output:
        done=touch(join(output_directory, 'sourmash-tax.done')),
        sourmash_taxonomy=sourmash_db_taxonomy,
    log:
        join(output_directory, 'sourmash-tax-convert.log')
    shell:
        'eval "$(pixi shell-hook -e sourmash)" && '
        'sourmash tax prepare -t {input.sourmash_taxonomy} -o {output.sourmash_taxonomy} -F sql &> {log}'

rule kaiju_download:
    output:
        done=touch(join(output_directory, 'kaiju-download.done'))
    log:
        abspath(join(output_directory, 'kaiju-download.log'))
    localrule:
        True
    shell:
        'mkdir -p {kaiju_db_original_dir} && '
        'cd {kaiju_db_original_dir} && '
        'wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_progenomes_2023-05-25.tgz '
        '-O kaiju_db_progenomes_2023-05-25.tgz &> {log}'

rule kaiju_extract:
    input:
        done=join(output_directory, 'kaiju-download.done'),
    output:
        done=touch(join(output_directory, 'kaiju.done')),
        fmi=kaiju_db_original_progenomes_fmi,
        nodes=kaiju_db_original_progenomes_nodes,
        names=kaiju_db_original_progenomes_names
    log:
        abspath(join(output_directory, 'kaiju-extract.log'))
    params:
        outdir = join(output_directory, 'kaiju')
    shell:
        'cd {params.outdir} && tar -xzf ../kaiju/kaiju_db_progenomes_2023-05-25.tgz &> {log}'

# map2b download disabled: its DownloadDB.py pulls the GTDB r202 database from
# figshare, whose links currently return empty 202 ("preparing") responses, so
# the download never completes. (Updating MAP2B v1.5 -> v1.8 refreshes the
# figshare IDs but the links still stall and the DB stays GTDB r202.)
# rule map2b_checkout:
#     output:
#         done=touch(join(output_directory, 'map2b-checkout.done')),
#     log:
#         join(output_directory, 'map2b-checkout.log')
#     shell:
#         'rm -rf {map2b_checkout_dir} && git clone --branch MAP2Bv1.5 https://github.com/sunzhengCDNM/MAP2B {map2b_checkout_dir} &> {log}'

# rule map2b:
#     input:
#         join(output_directory, 'map2b-checkout.done'),
#     output:
#         done=touch(join(output_directory, 'map2b.done')),
#         map2b_db=directory(map2b_db)
#     log:
#         join(output_directory, 'map2b.log')
#     shell:
#         'eval "$(pixi shell-hook -e map2b)" && '
#         'python3 {map2b_checkout_dir}/scripts/DownloadDB.py -l {map2b_checkout_dir}/config/GTDB.CjePI.database.list -d {map2b_checkout_dir}/database/GTDB &> {log}'

## metabuli database download doesn't work because it is via sharepoint, which gives an indirect link.
# rule metabuli_download:
#     output:
#         done=touch(join(output_directory, 'metabuli-download.done')),
#     log:
#         join(output_directory, 'metabuli-download.log')
#     shell:
#         "mkdir -p {metabuli_db_dir} && wget 'https://connectqutedu.sharepoint.com/:u:/s/metabuli_gtdb_207/EYk7N71mp-NAtET5_X_fBDABM6AC_DCbxGiDc2rdVVlNiw?e=3i6ak0' -O {metabuli_db_dir}/gtdb.tar.gz &> {log}"

# rule metabuli_extract:
#     input:
#         done=join(output_directory, 'metabuli-download.done'),
#     output:
#         done=touch(join(output_directory, 'metabuli.done')),
#         db=directory(metabuli_db)
#     log:
#         abspath(join(output_directory, 'metabuli-extract.log'))
#     shell:
#         'cd {metabuli_db_dir} && tar -xzf gtdb.tar.gz &> {log}'

rule motus_db:
    output:
        done=touch(join(output_directory, 'motus-download.done')),
    log:
        abspath(join(output_directory, 'motus-download.log'))
    localrule:
        True
    shell:
        'eval "$(pixi shell-hook -e motus)" && '
        'motus downloadDB &> {log}'

rule motus_db_move:
    input:
        done=join(output_directory, 'motus-download.done'),
    output:
        motus_db=directory(motus_db),
        done=touch(join(output_directory, 'motus-db.done')),
    log:
        abspath(join(output_directory, 'motus-db-move.log'))
    shell:
        # ls /mnt/hpccs01/work/microbiome/msingle/mess/124_singlem-benchmarking/.snakemake/conda/b06521a4ea0bdbb4dd2eabbe19701683_/lib/python3.9/site-packages/motus/db_mOTU/
        'eval "$(pixi shell-hook -e motus)" && '
        'bin/migrate_motus_db.py {output.motus_db}'

rule motus_gtdb:
    output:
        done=touch(join(output_directory, 'motus-gtdb.done')),
        motus_gtdb_tsv=motus_gtdb_tsv,
    log:
        join(output_directory, 'motus-gtdb.log')
    localrule:
        True
    shell:
        'wget https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv '
        '-O {motus_gtdb_tsv} &> {log}'

rule motus:
    input:
        join(output_directory, 'motus-db.done'),
        join(output_directory, 'motus-gtdb.done'),
    output:
        done=touch(join(output_directory, 'motus.done'))
    localrule:
        True

rule singlem_download:
    output:
        done=touch(join(output_directory, 'singlem-download.done')),
        singlem_metapackage_tgz=singlem_metapackage_tgz
    log:
        join(output_directory, 'singlem-download.log')
    localrule:
        True
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "wget 'https://zenodo.org/records/11107165/files/S4.1.0.GTDB_r207.metapackage_20240502.smpkg.zb.tar.gz?download=1' "
        "-O {output.singlem_metapackage_tgz} &> {log}"

rule singlem_extract:
    input:
        done=join(output_directory, 'singlem-download.done'),
        singlem_metapackage_tgz=singlem_metapackage_tgz,
    output:
        done=touch(join(output_directory, 'singlem.done')),
        singlem_metapackage=directory(singlem_metapackage)
    log:
        abspath(join(output_directory, 'singlem-extract.log'))
    params:
        singlem_metapackage_basename = os.path.basename(singlem_metapackage)
    shell:
        "bash -c 'cd {output_directory} && "
        "tar -xzf {params.singlem_metapackage_basename}.zb.tar.gz && "
        "mv -v {params.singlem_metapackage_basename}.zb/payload_directory ../{output.singlem_metapackage}' &> {log}"

rule weebill_db_convert:
    input:
        weebill_syl=weebill_syldb
    output:
        weebill_syl2db=weebill_db
    resources:
        mem_mb=64000,
        runtime=120
    log:
        join(output_directory, "singlem-regime3.weebill.log")
    shell:
        f"{weebill_binary} "
        "db-convert {input.weebill_syl} -o {output.weebill_syl2db} &> {log}"
        
rule metax_download:
    output:
        zip = join(output_directory, "taxdmp_2026-07-01.zip"),
        done = touch(join(output_directory, "metax-download.done"))
    log:
        join(output_directory, "metax-download.log")
    localrule:
        True
    shell:
        "wget 'https://s3.bi.denbi.de/cami3__human-gut-challenge/databases/taxdmp_2026-07-01.zip?AWSAccessKeyId=d4098cac781f4b079ab33853731fcc30&Signature=fl9Y128uTIjqbZLi2ZO5BqFCDHk%3D&Expires=1815734514' "
        "-O {output.zip} &> {log}"

rule metax_extract:
    input:
        join(output_directory, "taxdmp_2026-07-01.zip")
    output:
        touch(join(output_directory, "metax.done")),
        d1 = directory(join(output_directory, ".taxonkit")),
    log:
        join(output_directory, "metax-extract.log")
    shell:
        "mkdir -p {output.d1} && "
        "unzip {input} -d {output.d1} &> {log}"

#rule gtdb_download_bac120:
#    output:
#        done=touch(join(output_directory, 'gtdb-download-bac.done')),
#        tar = 'bac120_metadata_r207.tar.gz'
#    log:
#        join(output_directory, 'gtdb-bac.log')
#    shell:
#        'wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz -O {output.tar} &> {log}'

#rule gtdb_extract_bac120:
#    input:
#        done=join(output_directory, 'gtdb-download-bac.done'),
#        tar = 'bac120_metadata_r207.tar.gz'
#    output:
#        done=touch(join(output_directory, 'gtdb-bac.done')),
#    log:
#        join(output_directory, 'gtdb-extract.log')
#    shell:
#        'tar -xzf {input.tar} &> {log}'

#rule gtdb_download_ar53:
#    output:
#        done=touch(join(output_directory, 'gtdb-download-ar.done')),
#        tar = 'ar53_metadata_r207.tar.gz'
#    log:
#        join(output_directory, 'gtdb-ar.log')
#    shell:
#        'wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz -O {output.tar} &> {log}'

#rule gtdb_extract_ar53:
#    input:
#        done=join(output_directory, 'gtdb-download-ar.done'),
#        tar = 'ar53_metadata_r207.tar.gz'
#    output:
#        done=touch(join(output_directory, 'gtdb-ar.done')),
#    log:
#        join(output_directory, 'gtdb-extract-ar.log')
#    shell:
#        'tar -xzf {input.tar} &> {log}'

# benchmark datasets
rule shadow_genomes_download:
    output:
        done=touch(join(output_directory, 'shadow-genomes-download.done')),
    log:
        join(output_directory, 'shadow-genomes-download.log')
    localrule:
        True
    shell:
        "rm -f reference_genomes.tar.gz && "
        "bash -c 'cd {output_directory} && "
        "wget 'https://zenodo.org/records/11002998/files/reference_genomes.tar.gz?download=1' "
        "-O reference_genomes.tar.gz' &> {log}"
        
rule shadow_genomes_extract:
    input:
        done=join(output_directory, 'shadow-genomes-download.done'),
    output:
        done=touch(join(output_directory, 'shadow-extract.done')),
    log:
        join(output_directory, 'shadow-genomes.log')
    shell:
        "bash -c 'cd {output_directory} && tar xf reference_genomes.tar.gz' &> {log}"

rule reference_genomes_link:
    input:
        done=join(output_directory, 'shadow-extract.done'),
    output:
        done=touch(join(output_directory, 'shadow-genomes.done')),
    log:
        join(output_directory, 'reference-genomes-link.log')
    localrule:
        True
    shell:
        "bash -c 'ln -sf {output_directory}/reference_genomes && "
        "cd {output_directory} && "
        "ln -sf reference_genomes/shadow . && "
        "ln -sf reference_genomes/new_in_r214 .' &> {log}"

rule genome_pairs_fst_download:
    output:
        done=touch(join(output_directory, 'genome-pairs-fst-download.done')),
    log:
        join(output_directory, 'genome-pairs-fst-download.log')
    localrule:
        True
    shell:
        "bash -c 'cd 2_phylogenetic_novelty/genomes && "
        "datasets download genome accession --inputfile ../genome_accessions.txt' &> {log}"

rule genome_pairs_fst_extract:
    output: 
        done=touch(join(output_directory, 'genome-pairs-fst-extract.done'))
    log:
        join(output_directory, 'genome-pairs-fst-extract.log')
    shell:
        "bash -c 'cd 2_phylogenetic_novelty/genomes && "
        "{{ rm -r ncbi_dataset/data README.md md5sum.txt; unzip ncbi_dataset.zip; }} && "
        # Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
        "parallel --col-sep \"\t\" cp {{1}} {{2}} :::: ../genome_ncbi_names.tsv' &> {log}"

rule genome_pairs_snd_download:
    output:
        done=touch(join(output_directory, 'genome-pairs-snd-download.done')),
    log:
        join(output_directory, 'genome-pairs-snd-download.log')
    localrule:
        True
    shell:
        "bash -c 'cd 2_phylogenetic_novelty/genome_pairs && "
        "datasets download genome accession --inputfile ../genome_pairs_accessions.txt' &> {log}"

rule genome_pairs_snd_extract:
    output: 
        done=touch(join(output_directory, 'genome-pairs-snd-extract.done'))
    log:
        join(output_directory, 'genome-pairs-snd-extract.log')
    shell:
        "bash -c 'cd 2_phylogenetic_novelty/genome_pairs && "
        "{{ rm -r ncbi_dataset/data README.md md5sum.txt; unzip ncbi_dataset.zip; }} && "
        # Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
        "parallel --col-sep \"\t\" cp {{1}} {{2}} :::: ../genome_pairs_ncbi_names.tsv' &> {log}"

rule genome_pairs_download:
    input:
        fst = join(output_directory, 'genome-pairs-fst-extract.done'),
        snd = join(output_directory, 'genome-pairs-snd-extract.done')
    output:
        done=touch(join(output_directory, 'genome-pairs-download.done')),
    localrule:
        True

rule download_cami_marine_reads:
    output:
        "3_cami2_marine/simulation_short_read/marmgCAMI2_sample_{sample_number}_reads.tar.gz"
    log:
        "3_cami2_marine/simulation_short_read/marmgCAMI2_sample_{sample_number}_reads-download.log"
    localrule:
        True
    shell:
        "bash -c 'mkdir -p 3_cami2_marine/simulation_short_read && "
        "cd 3_cami2_marine/simulation_short_read && "
        "rm -f marmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz && "
        "wget https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule extract_cami_marine_reads:
    input:
        "3_cami2_marine/simulation_short_read/marmgCAMI2_sample_{sample_number}_reads.tar.gz"
    output:
        "3_cami2_marine/simulation_short_read/simulation_short_read/2018.08.15_09.49.32_sample_{sample_number}/reads/anonymous_reads.fq.gz"
    log:
        "3_cami2_marine/simulation_short_read/marmgCAMI2_sample_{sample_number}_reads-extract.log"
    shell:
        "bash -c 'cd 3_cami2_marine/simulation_short_read && "
        "tar -xzf marmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule split_cami_marine_reads:
    input:
        "3_cami2_marine/simulation_short_read/simulation_short_read/2018.08.15_09.49.32_sample_{sample_number}/reads/anonymous_reads.fq.gz"
    output:
        r1="3_cami2_marine/split_reads/marine{sample_number}.1.fq.gz",
        r2="3_cami2_marine/split_reads/marine{sample_number}.2.fq.gz",
        done=touch("3_cami2_marine/split_reads/marine{sample_number}.done")
    log:
        "3_cami2_marine/split_reads/marine{sample_number}.log"
    shell:
        "bash -c 'mkdir -p 3_cami2_marine/split_reads && "
        "zcat {input[0]} |./bin/deinterleave_fastq.sh {output.r1} {output.r2} compress' &> {log}"

rule bench2_genomes_download:
    output:
        touch("2_phylogenetic_novelty/reference_genomes-download.done"),
    log:
        "2_phylogenetic_novelty/reference_genomes-download.log"
    localrule:
        True
    shell:
        "cd 2_phylogenetic_novelty && "
        "wget 'https://zenodo.org/records/12525852/files/bench2_genomes.tar.gz?download=1' "
        "-O bench2_genomes.tar.gz &> ../{log}"

rule bench2_genomes_extract:
    input:
        "2_phylogenetic_novelty/reference_genomes-download.done"
    output:
        touch("2_phylogenetic_novelty/reference_genomes-extract.done"),
        d1 = directory("2_phylogenetic_novelty/genomes"),
        d2 = directory("2_phylogenetic_novelty/genome_pairs"),
    log:
        "2_phylogenetic_novelty/reference_genomes-extract.log"
    shell:
        """
        cd 2_phylogenetic_novelty && tar -xzf bench2_genomes.tar.gz &> ../{log}
        """

# Benchmarks 8 and 9 both mix genomes from the Zenodo shadow pool (known in GTDB
# r207) with genomes whose species is new in r214. The novel ones are not in that
# tarball, and unlike benchmark 2's they were never packaged, so fetch them from NCBI
# by accession. Which accession is listed matters: it is the one whose contig names
# the original reads were simulated from, which is decidable from the read names --
# RefSeq (GCF_) copies prefix their contigs NZ_, GenBank (GCA_) ones do not. Five of
# the seven are RefSeq; GCA_934272745.1 has no RefSeq copy, and GCA_021151785.1's
# reads carry bare JAGKWP0100000NN.1 contig names, so it was the GenBank copy.
rule bench89_novel_genomes_download:
    input:
        accessions = novel_r214_genome_accessions_file,
    output:
        fastas = novel_r214_genome_fastas,
    log:
        "novel_r214_genomes-download.log"
    shell:
        # datasets nests each genome as
        # ncbi_dataset/data/<ACC>/<ACC>_<asm>_genomic.fna; flatten to <ACC>.fna,
        # which is the name community.tsv refers to.
        """
        bash -c 'mkdir -p novel_r214_genomes && cd novel_r214_genomes && \
        datasets download genome accession --inputfile ../{input.accessions} && \
        unzip -o ncbi_dataset.zip && \
        while read acc; do cp ncbi_dataset/data/$acc/*.fna $acc.fna; done < ../{input.accessions}' &> {log}
        """

rule download_zymo_reads:
    input:
        zymo_reads

rule download_zymo_reference_genomes:
    output:
        done = touch("6_zymo_synthetic/zymo_refseq.v2-download.done"),
        zip = "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2.zip"
    log:
        "6_zymo_synthetic/reference_genomes-download.log"
    shell:
        "mkdir -p 6_zymo_synthetic && "
        "wget 'https://zenodo.org/records/3935737/files/ZymoBIOMICS.STD.refseq.v2.zip?download=1' "
        "-O {output.zip} &> {log}"

rule extract_zymo_reference_genomes:
    input:
        "6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2.zip"
    output:
        touch("6_zymo_synthetic/reference_genomes-extract.done"),
        d1 = directory("6_zymo_synthetic/ZymoBIOMICS.STD.refseq.v2/Genomes"),
    log:
        "6_zymo_synthetic/reference_genomes-extract.log"
    shell:
        """
        unzip {input} -d 6_zymo_synthetic/ &> {log}
        """

rule download_shakya:
    input:
        shakya_reads

rule download_fastq:
    output:
        r1="<bench_dir>/local_reads/{sample}.1.fq.gz",
        r2="<bench_dir>/local_reads/{sample}.2.fq.gz",
        done=touch("<bench_dir>/local_reads/{sample}.done")
    log:
        "<bench_dir>/local_reads/{sample}.log"
    threads:
        1
    pathvars:
        bench_dir = "7_shakya_synthetic"
    params:
        reads_dir = lambda wildcards, input, output: os.path.dirname(output.r1)
    shell:
        "pixi run -e kingfisher "
        "kingfisher get -r {wildcards.sample} "
        "--output_directory {params.reads_dir} "
        "-m ena-ftp prefetch -f fastq.gz "
        "&& mv {params.reads_dir}/{wildcards.sample}_1.fastq.gz {output.r1} "
        "&& mv {params.reads_dir}/{wildcards.sample}_2.fastq.gz {output.r2} "
        "&> {log}"

rule download_cami_strain:
    input:
        cami_strain_reads,
        cami_strain_coverage_definitions,
        "8_cami2_strain/genomes/"

rule download_cami_strain_reads:
    output:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads.tar.gz"
    log:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads-download.log"
    localrule: True
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain/short_read && " \
        "cd 8_cami2_strain/short_read && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule extract_cami_strain_reads:
    input:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads.tar.gz"
    output:
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/anonymous_reads.fq.gz",
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/reads_mapping.tsv.gz"
    log:
        "8_cami2_strain/short_read/strmgCAMI2_sample_{sample_number}_reads-extract.log"
    localrule: True
    shell:
        "bash -c " \
        "'cd 8_cami2_strain/short_read && " \
        "tar -xzf strmgCAMI2_sample_{wildcards.sample_number}_reads.tar.gz' &> {log}"

rule split_cami_strain_reads:
    input:
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/anonymous_reads.fq.gz"
    output:
        r1="8_cami2_strain/local_reads/strain{sample_number}.1.fq.gz",
        r2="8_cami2_strain/local_reads/strain{sample_number}.2.fq.gz",
        done=touch("8_cami2_strain/local_reads/strain{sample_number}.done")
    log:
        "8_cami2_strain/local_reads/strain{sample_number}.log"
    localrule: True
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain/local_reads && zcat {input[0]} | " \
        "bin/deinterleave_fastq.sh {output.r1} {output.r2} compress' &> {log}"

rule extract_cami_strain_read_mapping_counts:
    input:
        "8_cami2_strain/short_read/short_read/2018.09.07_11.43.52_sample_{sample_number}/reads/reads_mapping.tsv.gz"
    output:
        "8_cami2_strain/read_mapping_counts/strain{sample_number}.tsv"
    shell:
        "mkdir -p 8_cami2_strain/read_mapping_counts && "
        "pigz -cd {input[0]} |cut -f2 |sort |uniq -c "
        "|sed 's/^  *//g' > {output}"

rule download_cami_strain_genomes:
    output:
        "8_cami2_strain/strmgCAMI2_genomes.tar.gz"
    log:
        "8_cami2_strain/strmgCAMI2_genomes-download.log"
    localrule: True
    shell:
        "bash -c " \
        "'mkdir -p 8_cami2_strain && " \
        "cd 8_cami2_strain && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/strmgCAMI2_genomes.tar.gz' &> {log}"

rule extract_cami_strain_genomes:
    input:
        "8_cami2_strain/strmgCAMI2_genomes.tar.gz"
    output:
        directory("8_cami2_strain/genomes/")
    log:
        "8_cami2_strain/strmgCAMI2_genomes-extract.log"
    localrule: True
    shell:
        "mkdir -p 8_cami2_strain/genomes_extract && " 
        "tar -xzf {input[0]} -C 8_cami2_strain/genomes_extract && "
        "mv 8_cami2_strain/genomes_extract/short_read/source_genomes {output[0]} && "
        "rmdir 8_cami2_strain/genomes_extract/short_read && "
        "rmdir 8_cami2_strain/genomes_extract &> {log}"

rule extract_cami_strain_genome_lengths:
    input:
        "8_cami2_strain/genomes/"
    output:
        directory("8_cami2_strain/genome_lengths")
    log:
        "8_cami2_strain/genome_lengths.log"
    shell:
        "mkdir -p {output} && "
        "parallel 'cut -f1,2 {{}} > {output}/{{/.}}.tsv' ::: {input}/*.fai"

# setup has genome_to_id.tsv and coverage_newXX.tsv
rule download_cami_strain_setup:
    output:
        "8_cami2_strain/strmgCAMI2_setup.tar.gz",
    log:
        '8_cami2_strain/strmgCAMI2_setup-download.log'
    localrule: True
    shell:
        "bash -c 'cd 8_cami2_strain && " \
        "rm -f strmgCAMI2_setup.tar.gz && " \
        "wget https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_setup.tar.gz' &> {log}"
 
rule extract_cami_strain_setup:
    input:
        "8_cami2_strain/strmgCAMI2_setup.tar.gz"
    output:
        directory("8_cami2_strain/setup_extract/short_read"),
        touch("8_cami2_strain/strmgCAMI2_setup-extract.done")
    log:
        "8_cami2_strain/strmgCAMI2_setup-extract.log"
    localrule: True
    shell:
        "mkdir -p 8_cami2_strain/setup_extract && " \
        "tar -xzf {input[0]} -C 8_cami2_strain/setup_extract' &> {log}"

# # regenerate coverage definitions instead of using these ones in setup
# rule copy_cami_strain_coverage_definitions:
#     input:
#         setup_dir="8_cami2_strain/setup_extract/short_read"
#     output:
#         "8_cami2_strain/coverage_definitions_2/strain{sample_number}.tsv"
#     localrule: True
#     shell:
#         "mkdir -p 8_cami2_strain/coverage_definitions && "
#         "mv {input.setup_dir}/coverage_new{wildcards.sample_number}.tsv "
#         "{output[0]}"

rule generate_cami_strain_coverage_definitions:
    input:
        setup_dir="8_cami2_strain/setup_extract/short_read",
        genome_lengths_directory="8_cami2_strain/genome_lengths",
        read_mapping_counts="8_cami2_strain/read_mapping_counts/{sample}.tsv"
    output:
        "8_cami2_strain/coverage_definitions/{sample}.tsv"
    log:
        "8_cami2_strain/logs/generate_coverage_defintions-{sample}.log"
    localrule: True
    shell:
        "python3 {workflow.basedir}/bin/cami_coverage_definitions.py "
        "--genome-to-id {input.setup_dir}/genome_to_id.tsv " 
        "--genome-lengths-directory {input.genome_lengths_directory} " \
        "--read-mapping-counts {input.read_mapping_counts} " \
        "> {output[0]} 2> {log}"

## profiles not useful
# rule download_cami2_profiles:
#     output:
#         '8_cami2_strain/taxonomic_profiling_cami2.tar.gz'
#     log:
#         '8_cami2_strain/cami2_profiles-download.log'
#     shell:
#         "bash -c 'cd 8_cami2_strain && " \
#         "rm -f taxonomic_profiling_cami2.tar.gz && " \
#         "wget https://zenodo.org/records/5006866/files/taxonomic_profiling_cami2.tar.gz?download=1 -O taxonomic_profiling_cami2.tar.gz' &> {log}"

# rule extract_cami2_profiles:
#     input:
#         "8_cami2_strain/taxonomic_profiling_cami2.tar.gz"
#     output:
#         directory("8_cami2_strain/taxonomic_profiling_cami2/")
#     log:
#         "8_cami2_strain/taxonomic_profiling_cami2-extract.log"
#     shell:
#         "bash -c " \
#         "'cd 8_cami2_strain && " \
#         "tar -xzf taxonomic_profiling_cami2.tar.gz' &> {log}"

rule download_gtdbtk_r207_data:
    output:
        join(output_directory, 'gtdbtk_r207_v2_data.tar.gz')
    log:
        join(output_directory, "gtdbtk_r207_data-download.log")
    shell:
        "bash -c "
        f"'cd {output_directory} && "
        "wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz' &> {log}"

rule extract_gtdbtk_r207_data:
    input:
        join(output_directory, "gtdbtk_r207_v2_data.tar.gz")
    output:
        directory(join(output_directory, "release207_v2/"))
    log:
        join(output_directory, "gtdbtk_r207_v2_data-extract.log")
    shell:
        "bash -c " \
        f"'cd {output_directory} && " \
        "tar -xzf gtdbtk_r207_v2_data.tar.gz' &> {log}"

rule download_metax_profiles:
    output:
        done = touch("18_metax_gut/metax_profiles-download.done"),
        zip = "18_metax_gut/metax_benchmark_profiles.zip"
    log:
        "18_metax_gut/metax_profiles-download.log"
    shell:
        "wget 'https://zenodo.org/records/20128301/files/metax_benchmark_profiles.zip?download=1' "
        "-O {output.zip} &> {log}"

rule extract_metax_profiles:
    input:
        zip = "18_metax_gut/metax_benchmark_profiles.zip"
    output:
        done = touch("18_metax_gut/metax_profiles-extract.done"),
        dir = directory("18_metax_gut/metax_benchmark_profiles")
    log:
        "18_metax_gut/metax_profiles-extract.log"
    shell:
        "unzip {input.zip} -d 18_metax_gut &> {log}"

rule metax_profiles_to_coverages:
    input:
        dir = "18_metax_gut/metax_benchmark_profiles"
    output:
        directory("18_metax_gut/coverages")
    log:
        "18_metax_gut/metax_profiles_to_condensed.log"
    shell:
        "mkdir -p {output} && "
        "pixi run -e cami-opal-lib "
        "python3 bin/biobox_to_coverages.py "
        "--profile {input.dir}/gut/Gold_standard.profile "
        "--outdir {output} &> {log}"

rule metax_coverages_to_taxids:
    input:
        [f"18_metax_gut/coverages/sample_{sample_number}.tsv" for sample_number in range(10)]
    output:
        "18_metax_gut/taxids.txt"
    shell:
        'cat {input} | cut -f1 | sed "s/\..\+$//" > {output}'

rule metax_taxids_to_genome_accessions:
    input:
        "18_metax_gut/taxids.txt"
    output:
        "18_metax_gut/taxids_to_genome_accessions.txt"
    log:
        "18_metax_gut/taxids_to_genome_accessions.log"
    shell:
        'datasets summary genome taxon --inputfile {input} --assembly-level complete --as-json-lines | '
        'dataformat tsv genome --fields accession,assminfo-refseq-category,assminfo-name,organism-tax-id,assminfo-biosample-description-organism-tax-id 2> {log} > {output}'

# convert ncbi to gtdb
rule metax_download_ncbi_to_gtdb_bacteria_translation:
    output:
        "18_metax_gut/ncbi_vs_gtdb_r232_bacteria.xlsx"
    shell:
        "wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/ncbi_vs_gtdb_r232_bacteria.xlsx "
        "-O 18_metax_gut"

rule metax_download_ncbi_to_gtdb_archaea_translation:
    output:
        "18_metax_gut/ncbi_vs_gtdb_r232_archaea.xlsx"
    shell:
        "wget https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/ncbi_vs_gtdb_r232_archaea.xlsx "
        "-O 18_metax_gut"

rule metax_ncbi_to_gtdb_translation:
    input:
        bac = "18_metax_gut/ncbi_vs_gtdb_r232_bacteria.xlsx",
        ar = "18_metax_gut/ncbi_vs_gtdb_r232_archaea.xlsx"
    output:
        "18_metax_gut/ncbi_to_gtdb_r232_mapping.tsv"
    log:
        "18_metax_gut/ncbi_to_gtdb_translation.log"
    shell:
        "pixi run -e fastexcel "
        "python3 bin/ncbi_to_gtdb_translation.py --bac {input.bac} --ar {input.ar} "
        "--output-file {output} &> {log}"

rule metax_coverages_to_genomes_list:
    input:
        [f"18_metax_gut/coverages/sample_{sample_number}.tsv" for sample_number in range(10)]
    output:
        "18_metax_gut/genomes_to_ids.txt"
    shell:
        """
        : > {output} && 
        for f in {input}; do
          cut -f1 "$f" | sed 's/\(.*\)/genomes/\1.fasta\t\1/' >> {output}
        done
        """

rule all_metax_gut_reads:
    input:
        metax_gut_reads

rule download_metax_gut_reads:
    output:
        done = touch("18_metax_gut/metax_gut_reads-download.done"),
        dir = directory("18_metax_gut/gut_dataset")
    log:
        "18_metax_gut/metax_gut_reads-download.log"
    shell:
        "wget 'https://research.bifo.helmholtz-hzi.de/downloads/metax/benchmark_datasets/gut/' "
        "-r -np -R \"index.html*\" -P 18_metax_gut/gut_dataset &> {log}"

rule split_metax_gut_reads:
    input:
        "18_metax_gut/gut_dataset/research.bifo.helmholtz-hzi.de/downloads/metax/benchmark_datasets/gut/sample_{sample_number}/reads/anonymous_reads.fq.gz"
    output:
        r1="18_metax_gut/split_reads/sample_{sample_number}.1.fq.gz",
        r2="18_metax_gut/split_reads/sample_{sample_number}.2.fq.gz",
        done=touch("18_metax_gut/split_reads/sample_{sample_number}.done")
    log:
        "18_metax_gut/split_reads/sample_{sample_number}.log"
    shell:
        "mkdir -p 18_metax_gut/split_reads && "
        "zcat {input[0]} |./bin/deinterleave_fastq.sh {output.r1} {output.r2} compress &> {log}"

# # don't need? datasets takes taxon ids directly     
# rule download_cami2_accession2taxid:
#     output:
#         done = touch("18_metax_gut/ncbi_accession2taxid-download.done"),
#         tar = "18_metax_gut/ncbi_taxonomy_accession2taxid.tar"
#     log:
#         "18_metax_gut/ncbi_accession2taxid-download.log"
#     shell:
#         "wget 'https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/ncbi_taxonomy_accession2taxid.tar' "
#         "-P 18_metax_gut &> {log}"
#         #"wget 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz' "

# # don't need?
# rule extract_cami2_accession2taxid:
#     input:
#         "18_metax_gut/ncbi_taxonomy_accession2taxid.tar"
#     output:
#         touch("18_metax_gut/ncbi-accession2taxid-extract.done"),
#         d1 = directory("18_metax_gut/ncbi-accession2taxid"),
#     log:
#         "18_metax_gut/ncbi-accession2taxid-extract.log"
#     shell:
#         "tar -xf {input} -C 18_metax_gut &> {log}"

# # don't need?
# rule download_cami2_refseqdb:
#     output:
#         touch("18_metax_gut/ncbi_refseqdb-download.done"),
#     log:
#         "18_metax_gut/ncbi_refseqdb-download.log"
#     shell:
#         "wget 'https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_2_DATABASES/RefSeq_genomic_20190108.tar' "
#        "-P 18_metax_gut &> {log}"

# 19 clinical stool
rule download_clinical_stool_reads:
    input:
        clinical_stool_reads

# rule annotate_clinical_stool_samples:
# output:
#     "19_clinical_stool/PRJNA786578_annotate.tsv"
# log:
#     "19_clinical_stool/samples-annotate.log"
# shell:
#     "pixi run -e kingfisher " \
#     "kingfisher annotate -p PRJNA786578 " \
#      "--output-file {output} --output-format tsv &> {log}"

use rule download_fastq as download_clinical_stool_sample with:
    pathvars:
        bench_dir = "19_clinical_stool"

# # can't wget from this url
# rule download_clinical_stool_metadata:
#     output:
#         zip="19_clinical_stool/microorganisms-10-00441-s001.zip"
#     log:
#         "19_clinical_stool/metadata-download.log"
#     shell:
#         "wget -O {output.zip} https://pmc.ncbi.nlm.nih.gov/articles/instance/8880012/bin/microorganisms-10-00441-s001.zip &> {log}"

rule extract_clinical_stool_metadata:
    input:
        zip="19_clinical_stool/microorganisms-10-00441-s001.zip"
    output:
        dir=directory("19_clinical_stool/metadata")
    log:
        "19_clinical_stool/metadata-extract.log"
    shell:
        "unzip {input.zip} -d {output.dir} &> {log}"
