#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, Optional, Set

from gtdb_utils import iter_rep_accessions, normalize_genome_accession, strip_accession_version


def parse_summary_header(handle) -> Optional[Dict[str, int]]:
    for raw_line in handle:
        if raw_line.startswith("#assembly_accession\t"):
            columns = raw_line.rstrip("\n").lstrip("#").split("\t")
            return {name: idx for idx, name in enumerate(columns)}
    return None


def write_taxdump(nodes_dmp: Path, names_dmp: Path, names_by_taxid: Dict[int, str]) -> None:
    nodes_dmp.parent.mkdir(parents=True, exist_ok=True)
    with nodes_dmp.open("w", encoding="utf-8") as out:
        out.write("1\t|\t1\t|\tno rank\t|\t\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")
        for taxid in sorted(names_by_taxid):
            out.write(f"{taxid}\t|\t1\t|\tspecies\t|\t\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n")

    with names_dmp.open("w", encoding="utf-8") as out:
        out.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        for taxid in sorted(names_by_taxid):
            out.write(f"{taxid}\t|\t{names_by_taxid[taxid]}\t|\t\t|\tscientific name\t|\n")


def iter_summary_rows(summary_path: Path) -> Iterable[Dict[str, str]]:
    with summary_path.open("r", encoding="utf-8", errors="replace") as handle:
        header = parse_summary_header(handle)
        if header is None:
            raise SystemExit(f"Could not find #assembly_accession header in {summary_path}")
        required = ["assembly_accession", "taxid", "species_taxid", "organism_name"]
        for field in required:
            if field not in header:
                raise SystemExit(f"Missing column {field} in {summary_path}")
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            yield {
                "assembly_accession": row[header["assembly_accession"]],
                "taxid": row[header["taxid"]],
                "species_taxid": row[header["species_taxid"]],
                "organism_name": row[header["organism_name"]],
            }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a minimal NCBI-style taxdump using assembly_summary accessions and species_taxid values."
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--genomes-root", type=Path, help="Root directory containing GTDB rep genomes (*.fna[.gz])")
    source_group.add_argument("--tarball", type=Path, help="GTDB representative-genomes tarball")
    parser.add_argument(
        "--assembly-summary",
        dest="assembly_summaries",
        action="append",
        required=True,
        type=Path,
        help="assembly_summary file; provide once for GenBank and once for RefSeq when available",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for nodes.dmp/names.dmp")
    parser.add_argument("--acc2taxid", required=True, type=Path, help="Output TSV mapping accession<TAB>species_taxid")
    args = parser.parse_args()

    rep_accessions = iter_rep_accessions(genomes_root=args.genomes_root, tarball=args.tarball)
    if not rep_accessions:
        source = args.genomes_root if args.genomes_root is not None else args.tarball
        raise SystemExit(f"No representative genomes found from: {source}")

    acc_to_taxid = {}  # type: Dict[str, int]
    names_by_taxid = {}  # type: Dict[int, str]
    seen = 0
    kept = 0

    for summary_path in args.assembly_summaries:
        for row in iter_summary_rows(summary_path):
            seen += 1
            acc = normalize_genome_accession(row["assembly_accession"])
            acc_nover = strip_accession_version(acc)
            if acc not in rep_accessions and acc_nover not in rep_accessions:
                continue
            species_taxid = row["species_taxid"].strip()
            fallback_taxid = row["taxid"].strip()
            if species_taxid.isdigit() and int(species_taxid) > 0:
                taxid = int(species_taxid)
            elif fallback_taxid.isdigit() and int(fallback_taxid) > 0:
                taxid = int(fallback_taxid)
            else:
                continue
            name = row["organism_name"].strip() or acc
            acc_to_taxid[acc] = taxid
            acc_to_taxid[acc_nover] = taxid
            names_by_taxid.setdefault(taxid, name)
            kept += 1

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_taxdump(args.outdir / "nodes.dmp", args.outdir / "names.dmp", names_by_taxid)

    args.acc2taxid.parent.mkdir(parents=True, exist_ok=True)
    with args.acc2taxid.open("w", encoding="utf-8") as out:
        for acc in sorted(acc_to_taxid):
            out.write(f"{acc}\t{acc_to_taxid[acc]}\n")

    missing = 0
    for acc in rep_accessions:
        if acc in acc_to_taxid:
            continue
        missing += 1
    print(
        f"[make_species_taxdump_from_assembly_summary] rep_accessions={len(rep_accessions)} kept={kept}/{seen} "
        f"species_nodes={len(names_by_taxid)} missing_rep_accessions={missing}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
