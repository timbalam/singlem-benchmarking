#! /usr/bin/env python

# %%
import argparse
import logging
import sys
import os, re

import polars as pl
import tempfile
import extern



# %%
def read_gtdbtk(output_directory, taxonomy_only=True, remove_empty_ranks=False):
    if not taxonomy_only:
        raise NotImplementedError("Only taxonomy is supported for now")
    taxonomies = {}
    bac_taxonomy_file = os.path.join(output_directory, 'gtdbtk.bac120.summary.tsv')
    logging.debug('Reading taxonomy from %s' % bac_taxonomy_file)
    d = pl.read_csv(bac_taxonomy_file, separator='\t')
    if remove_empty_ranks:
        empty_ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    for row in d.rows(named=True):
        tax = row['classification']
        if remove_empty_ranks:
            tax = ';'.join([x for x in tax.split(';') if x.strip() not in empty_ranks])
        taxonomies[row['user_genome']] = tax
    logging.debug("Read %d taxonomies from Bacteria" % len(taxonomies))

    # Archaea
    arc_taxonomy_file = os.path.join(output_directory, 'gtdbtk.ar53.summary.tsv')
    logging.debug('Reading taxonomy from %s' % arc_taxonomy_file)
    d = pl.read_csv(arc_taxonomy_file, separator='\t')
    num_archaea = 0
    for row in d.rows(named=True):
        tax = row['classification']
        if remove_empty_ranks:
            tax = ';'.join([x for x in tax.split(';') if x.strip() not in empty_ranks])
        taxonomies[row['user_genome']] = tax
        num_archaea += 1
    logging.debug("Read %d new archaeal taxonomies, so %d total" % (num_archaea, len(taxonomies)))
    return taxonomies


# %%markdown
if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--coverage-file', required=True, help='Path to coverage file')
    parent_parser.add_argument('--known-genome-list', required=True, help='Path to genome list')
    # (env)cl5n007:20221031:~/m/msingle/mess/115_camisim_ish_benchmarking$ \ls -f ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes |grep GC |sed 's=\(.*\).fna=\1\t/work/microbiome/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes\1.fna=' >shadow_genome_paths.csv
    parent_parser.add_argument('--novel-genome-gtdbtk-output', required=True, help='Path to where new genomes have been run through gtdbtk')
    parent_parser.add_argument('--novel-genome-list', required=True, help='Path to genome list')
    parent_parser.add_argument('--novelty-ratio', required=True, help='Ratio of rank to parent rank novelty')
    parent_parser.add_argument('--gtdb-bac-metadata', required=True, help='Path to GTDB metadata file, for genome length')
    parent_parser.add_argument('--gtdb-ar-metadata', required=True, help='Path to GTDB metadata file, for genome length')
    parent_parser.add_argument('--output-condensed', required=True, help='Path to output file in singlem condensed format')
    parent_parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parent_parser.add_argument('--art', required=True, help='Path to ART binary (art_illumina)')

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


    # %%
    # class Args:
    #     read1 = 'r1.fq.gz'
    #     read2 = 'r2.fq.gz'
    #     coverage_file = '../4_complex_and_novel/coverage_definitions/coverage0.tsv'
    #     known_genome_list = '../1_novel_strains/shadow_genome_paths.csv'
    #     novel_genome_gtdbtk_output = '../4_complex_and_novel/gtdbtk_batchfile.random1000.gtdbtk_r207'
    #     novel_genome_list = '../4_complex_and_novel/gtdbtk_batchfile.random1000.csv'
    #     gtdb_bac_metadata = '../bac120_metadata_r207.tsv'
    #     gtdb_ar_metadata = "../ar53_metadata_r207.tsv"
    #     output_condensed = 'output_condensed.tsv'
    #     threads = 1
    #     art = 'art_illumina'
    #     novelty_ratio = 0.5
    # args = Args()

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # %%

    output1 = os.path.abspath(args.read1)
    output2 = os.path.abspath(args.read2)
    output_condensed = os.path.abspath(args.output_condensed)


# %%

    # Read coverages
    coverages = pl.read_csv(args.coverage_file,
                            separator = '\t', has_header = False,
                            new_columns = ['otu', 'coverage'])
    coverages = coverages.filter(pl.col('coverage') > 0)
    logging.info(f"Read {len(coverages)} coverages > 0.")


# %%
    # Remove RNODE ones which are plasmids etc.
    coverages = coverages.filter(~pl.col('otu').str.contains('RNODE'))
    logging.info(f"After removing plasmids etc, {len(coverages)} coverages > 0 remain.")


# %%

    genomes = (
        pl.read_csv(args.known_genome_list,
                    separator = '\t', has_header = False,
                    new_columns = ['genome', 'path'])
            .with_columns(
                # Make paths relative to input file
                pl.col('path').map_elements(
                    lambda x: os.path.normpath(os.path.join(
                        os.getcwd(),
                        os.path.dirname(args.known_genome_list),
                        x
                    )),
                    return_dtype = pl.String()
                )
            )
    )
    logging.info(f"Read {len(genomes)} genome fasta paths.")


# %%
    bac = pl.read_csv(args.gtdb_bac_metadata, separator = '\t',
                      infer_schema_length = 100000,
                      ignore_errors = True)
    ar = pl.read_csv(args.gtdb_ar_metadata, separator = '\t',
                     infer_schema_length = 100000,
                     ignore_errors = True)
    metadata = pl.concat([
        bac.select('accession', 'genome_size', 'gtdb_taxonomy'),
        ar.select('accession', 'genome_size', 'gtdb_taxonomy'),
    ])
    logging.info(f"Read {len(metadata)} GTDB metadata entries.")


# %%

    # get rid of GB_, RS_
    metadata = metadata.with_columns(pl.col('accession').str.slice(3).alias('genome'))

# %%

    # Shuffle genomes order so we get randomness
#    metadata = metadata.sample(fraction=1)

    known_info = (
        genomes.join(metadata, on = 'genome', how = 'inner')
            .select('path', 'genome', pl.col('gtdb_taxonomy').alias('taxonomy'))
    )

    r207_taxonomy = read_gtdbtk(args.novel_genome_gtdbtk_output, remove_empty_ranks=True)

# %%

    # Read list of genome paths
    novel_genome_list = (
        pl.read_csv(args.novel_genome_list, separator = '\t',
                    has_header = False,
                    new_columns = ['path', 'genome'])
            .with_columns(
                # Make paths relative to input file
                pl.col('path').map_elements(
                    lambda x: os.path.normpath(os.path.join(
                        os.getcwd(),
                        os.path.dirname(args.novel_genome_list),
                        x
                    )),
                    return_dtype = pl.String()
                )
            )
    )


# %%

    # Merge with GTDBTK output
    novel_info = novel_genome_list.with_columns(
        pl.col('genome').replace(r207_taxonomy).alias('taxonomy')
    ).with_columns(
        pl.col('taxonomy')
        .str.extract("(^d|;[pcofgs])__[^;]+$").alias("known_at")
    ).with_columns(
        pl.col('known_at')
        .replace_strict({";s": 0, ";g": 1, ";f": 2, ";o": 3, ";c": 4, ";p": 5, "d": 6})
        .alias('known_at')
    ).filter(pl.col('known_at') > 0)


# %%

    logging.info(f"Read {len(novel_info)} novel genomes.")

    #ranks = "rdpcofgs"
    sum_weights = sum(float(args.novelty_ratio) ** n for n in range(6))
    n_rank = (
        novel_info.group_by('known_at')
        .len("max_new")
        .with_columns(
            r = pl.lit(float(args.novelty_ratio)),
            n = pl.lit(len(coverages))
        )
        .with_columns(
            weight_new = pl.col("r") ** (pl.col("known_at") + 1)
        )
        .with_columns(
            frac_new = pl.col("weight_new") / sum_weights
        )
        .with_columns(
            n_new_r = (pl.col("n") * pl.col("frac_new")).round()
        )
        .with_columns(
            n_new = pl.min_horizontal(
                pl.col("n_new_r"),
                pl.col("max_new")
            )
        )
    )
    
    n_known = len(coverages) - n_rank["n_new"].sum()

# %%
    logging.info(f"Choosing {n_rank["n_new"].sum()} novel genomes and {n_known} known genomes.")
    novel_info_new = novel_info.join(
        n_rank.select("known_at", "n_new"),
        on = "known_at",
        how = "left",
        validate = "m:1"
    )
    chosen_df = pl.concat([
        known_info.sample(n_known),
        novel_info_new
        .with_columns(
            shuff_loc = pl.int_range(pl.len()).shuffle().over(pl.col("known_at"))
        )
        .filter(pl.col("shuff_loc") < pl.col("n_new"))
        .drop("shuff_loc", "known_at", "n_new")
    ])
    
    # Add coverage column
    chosen_df = chosen_df.with_columns(
        pl.lit(coverages['coverage']).alias('coverage')
    )
    

# %%
    read_length = 150

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = []
        with open(output_condensed, 'w') as f:
            f.write("sample\tcoverage\ttaxonomy\n")
            tax_to_coverage = {}

            for i, (fasta, genome_id, taxonomy, coverage) in enumerate(chosen_df.rows()):
                if taxonomy not in tax_to_coverage:
                    tax_to_coverage[taxonomy] = 0
                tax_to_coverage[taxonomy] += coverage

                sim_commands.append(
                    f"{args.art} -ss HSXt -i {fasta} -p -l {read_length} -f {coverage} -m 400 -s 10 -o simulated_reads/{i}. &>/dev/null"
                )

            for tax, cov in tax_to_coverage.items():
                f.write(f"{os.path.basename(args.coverage_file)}\t{cov}\t{tax}\n")
                
        logging.info(f"Simulating {len(sim_commands)} genomes ..")
        extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output1))
        extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output2))
    


