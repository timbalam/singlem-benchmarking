rule simulate_art_community_and_reads:
    input:
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


rule simulate_art_reads:
    output:
        r1=config["simulated"]

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
    