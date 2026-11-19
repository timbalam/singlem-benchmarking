#!/usr/bin/env bash
set -euo pipefail

R="${R:-226}"
THREADS="${THREADS:-16}"
NO_MASKING="${NO_MASKING:-1}"
CHUNK_SIZE_BYTES="${CHUNK_SIZE_BYTES:-20000000000}"
TAXONOMY_MODE="${TAXONOMY_MODE:-auto}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$PWD/kraken2_gtdb_r${R}_reps_work}"
TARBALL="${TARBALL:-$WORKDIR/gtdb_genomes_reps_r${R}.tar.gz}"
BAC_TAX="${BAC_TAX:-$WORKDIR/bac120_taxonomy_r${R}.tsv.gz}"
AR_TAX="${AR_TAX:-$WORKDIR/ar53_taxonomy_r${R}.tsv.gz}"
ASSEMBLY_SUMMARY_GENBANK="${ASSEMBLY_SUMMARY_GENBANK:-$WORKDIR/assembly_summary_genbank.txt}"
ASSEMBLY_SUMMARY_REFSEQ="${ASSEMBLY_SUMMARY_REFSEQ:-$WORKDIR/assembly_summary_refseq.txt}"
TAXDUMP_DIR="${TAXDUMP_DIR:-$WORKDIR/taxonomy}"
ACC2TAXID="${ACC2TAXID:-$PWD/acc2taxid.tsv}"
KTAXONOMY="${KTAXONOMY:-$PWD/ktaxonomy.tsv}"
DB_DIR="${DB_DIR:-$PWD/kraken2_gtdb_r${R}_reps_db}"
KRAKEN2_BUILD="${KRAKEN2_BUILD:-kraken2-build}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"

mkdir -p "$WORKDIR"

if [[ ! -f "$TARBALL" ]]; then
  echo "Missing $TARBALL"
  echo "Set TARBALL=... or place gtdb_genomes_reps_r${R}.tar.gz in $WORKDIR"
  exit 1
fi

taxonomy_mode_resolved="$TAXONOMY_MODE"
if [[ "$taxonomy_mode_resolved" == "auto" ]]; then
  if [[ -f "$BAC_TAX" && -f "$AR_TAX" ]]; then
    taxonomy_mode_resolved="gtdb"
  elif [[ -f "$ASSEMBLY_SUMMARY_GENBANK" || -f "$ASSEMBLY_SUMMARY_REFSEQ" ]]; then
    taxonomy_mode_resolved="assembly"
  else
    echo "Could not resolve TAXONOMY_MODE automatically."
    echo "Provide both $BAC_TAX and $AR_TAX for full GTDB taxonomy, or provide assembly summary files for fallback mode."
    exit 1
  fi
fi

echo "[config] R=$R THREADS=$THREADS NO_MASKING=$NO_MASKING CHUNK_SIZE_BYTES=$CHUNK_SIZE_BYTES TAXONOMY_MODE=$taxonomy_mode_resolved"
echo "[config] WORKDIR=$WORKDIR"
echo "[config] TARBALL=$TARBALL"
echo "[config] DB_DIR=$DB_DIR"
echo "[config] ACC2TAXID=$ACC2TAXID"
echo "[config] KTAXONOMY=$KTAXONOMY"

rm -rf "$DB_DIR/taxonomy" "$DB_DIR/library"
mkdir -p "$DB_DIR/taxonomy" "$DB_DIR/library/added"

if [[ "$taxonomy_mode_resolved" == "gtdb" ]]; then
  if [[ ! -f "$BAC_TAX" || ! -f "$AR_TAX" ]]; then
    echo "GTDB taxonomy mode requested, but missing $BAC_TAX or $AR_TAX"
    exit 1
  fi
  echo "[taxdump] building full GTDB taxdump"
  python3 -u "$HERE/make_gtdb_taxdump.py" \
    --bac120-taxonomy "$BAC_TAX" \
    --ar53-taxonomy "$AR_TAX" \
    --tarball "$TARBALL" \
    --outdir "$TAXDUMP_DIR" \
    --acc2taxid "$ACC2TAXID"
elif [[ "$taxonomy_mode_resolved" == "assembly" ]]; then
  summary_args=()
  if [[ -f "$ASSEMBLY_SUMMARY_GENBANK" ]]; then
    summary_args+=(--assembly-summary "$ASSEMBLY_SUMMARY_GENBANK")
  fi
  if [[ -f "$ASSEMBLY_SUMMARY_REFSEQ" ]]; then
    summary_args+=(--assembly-summary "$ASSEMBLY_SUMMARY_REFSEQ")
  fi
  if [[ ${#summary_args[@]} -eq 0 ]]; then
    echo "Assembly-summary fallback requested, but no assembly summary files were found."
    exit 1
  fi
  echo "[taxdump] building assembly-summary fallback taxonomy"
  python3 -u "$HERE/make_species_taxdump_from_assembly_summary.py" \
    --tarball "$TARBALL" \
    "${summary_args[@]}" \
    --outdir "$TAXDUMP_DIR" \
    --acc2taxid "$ACC2TAXID"
else
  echo "Unknown TAXONOMY_MODE: $taxonomy_mode_resolved"
  exit 1
fi

cp -f "$TAXDUMP_DIR/nodes.dmp" "$DB_DIR/taxonomy/nodes.dmp"
cp -f "$TAXDUMP_DIR/names.dmp" "$DB_DIR/taxonomy/names.dmp"

echo "[library] materializing Kraken2 library chunks"
python3 -u "$HERE/materialize_kraken_library.py" \
  --tarball "$TARBALL" \
  --acc2taxid "$ACC2TAXID" \
  --outdir "$DB_DIR/library/added" \
  --chunk-size-bytes "$CHUNK_SIZE_BYTES"

build_args=(--db "$DB_DIR" --threads "$THREADS" --build)
if [[ "$NO_MASKING" == "1" ]]; then
  build_args+=(--no-masking)
fi

echo "[db] building k-mer database"
"$KRAKEN2_BUILD" "${build_args[@]}"

echo "[db] cleaning intermediate files"
"$KRAKEN2_BUILD" --db "$DB_DIR" --clean

# ktaxonomy.tsv for the winnow-tax mono-lineage rescue step: pipe-delimited
# taxid|parent|rank_code|depth|name, derived from the taxdump. Read from
# $TAXDUMP_DIR (the work taxonomy/, which kraken2-build --clean does NOT touch).
echo "[ktaxonomy] building $KTAXONOMY from taxdump $TAXDUMP_DIR"
python3 -u "$HERE/build_ktaxonomy.py" "$TAXDUMP_DIR" "$KTAXONOMY"

echo "[done] Kraken2 GTDB r${R} DB: $DB_DIR"
echo "[done] sidecars: acc2taxid=$ACC2TAXID  ktaxonomy=$KTAXONOMY"
