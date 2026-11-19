# Kraken2 database from GTDB R226 representative genomes

Builds a workable **Kraken2** database (plus the `acc2taxid.tsv` and `ktaxonomy.tsv`
sidecars used by winnow-tax) directly from the local GTDB representative-genome tarball
`gtdb_genomes_reps_r226.tar.gz` — no per-genome `kraken2-build --add-to-library` calls.

---

## What the build does

`build_kraken2_gtdb_r226_reps.sh` runs, in order:

1. **Taxdump + acc2taxid** — `make_gtdb_taxdump.py` reads the GTDB taxonomy
   (`bac120_taxonomy_r226.tsv.gz` + `ar53_taxonomy_r226.tsv.gz`) and the reps tarball and
   writes an NCBI-style `nodes.dmp` / `names.dmp` taxdump plus **`acc2taxid.tsv`**
   (accession → taxid).
2. **Library** — `materialize_kraken_library.py` streams the reps genomes into chunked
   `library/*.fna` files with embedded `kraken:taxid|...` headers.
3. **Hash** — `kraken2-build --build` (optionally `--no-masking`), then `--clean`.
4. **ktaxonomy.tsv** — `build_ktaxonomy.py` derives the pipe-delimited
   `taxid | parent | rank_code | depth | name` file from the taxdump (read from the work
   `taxonomy/`, which `--clean` does not touch). **This is produced by the build script
   now** — it is *not* downloaded.

## Taxonomy modes (`TAXONOMY_MODE`)

- `gtdb` — full GTDB taxonomy if `bac120_taxonomy_r226.tsv.gz` and `ar53_taxonomy_r226.tsv.gz`
  are present (GTDB names/ranks).
- `assembly` — fallback using local `assembly_summary_genbank.txt` / `assembly_summary_refseq.txt`;
  yields a usable but **species-only** tree under root.
- `auto` (default) — prefer full GTDB taxonomy, else fall back to assembly summaries.

## Inputs

Expected in `./kraken2_gtdb_r226_reps_work/` by default:

- `gtdb_genomes_reps_r226.tar.gz`  *(required)*
- For full GTDB taxonomy: `bac120_taxonomy_r226.tsv.gz`, `ar53_taxonomy_r226.tsv.gz`
- For assembly fallback: `assembly_summary_genbank.txt` and/or `assembly_summary_refseq.txt`

## Build

From the directory where you want the database outputs written:

```bash
export THREADS=16
export NO_MASKING=1
export TAXONOMY_MODE=auto

bash /path/to/winnow-tax/kraken_gtdb_r226_build_scripts/build_kraken2_gtdb_r226_reps.sh
```

Set `KRAKEN2_BUILD=/path/to/kraken2-build` if `kraken2-build` is not on
`PATH`.

## Outputs

- Final Kraken2 DB: `./kraken2_gtdb_r226_reps_db/`
- Sidecars in the launch directory: **`./acc2taxid.tsv`**, **`./ktaxonomy.tsv`**
- Other work products (taxdump, library) under `./kraken2_gtdb_r226_reps_work/`

Override locations with env vars: `WORKDIR`, `TARBALL`, `BAC_TAX`, `AR_TAX`, `DB_DIR`,
`ACC2TAXID`, `KTAXONOMY`, `CHUNK_SIZE_BYTES`.

---

## Resource planning (read before launching)

A full GTDB Kraken build is **heavy**:

- **RAM:** plan for **> 1 TB**. A previous **1000 GB** attempt was **OOM-killed during the
  final `kraken2-build` `build_db` step** — request more than 1 TB for the build node.
- **Wall time:** ~**1.5 days**.
- **Disk:** the final database is large; keep ample free space for the work dir + DB.
- Tip: `CHUNK_SIZE_BYTES` tunes the library chunk size; `kraken2-build --clean` (run
  automatically) removes intermediate library files after a successful build.

> **Run-time memory (using the DB, not building it):** Kraken2 loads the whole hash into
> RAM, so winnow-tax runs also need a big-memory node (~700 GB–1 TB) — see the container
> README at the repository root.

## Notes

- Full GTDB mode gives GTDB names/ranks; the assembly-summary fallback is local-only and
  more robust here but produces a species-only tree.
- This DB is the one used (and validated via the winnow-tax spike / CAMI3 / pathogen
  benchmarks) for real-sample profiling. (An earlier status note paused further build
  work over the resource cost above; that resource guidance still applies — provision
  accordingly.)
