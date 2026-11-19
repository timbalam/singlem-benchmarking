#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from gtdb_utils import iter_rep_accessions, normalize_genome_accession, strip_accession_version


RANK_BY_PREFIX = {
    "d__": "superkingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}


def find_lineage_column(fields):  # type: (List[str]) -> Optional[str]
    for value in fields:
        if "d__" in value and ";" in value:
            return value
    return None


def parse_lineage(lineage):  # type: (str) -> List[Tuple[str, str]]
    result = []  # type: List[Tuple[str, str]]
    for token in lineage.strip().split(";"):
        token = token.strip()
        if len(token) < 4:
            continue
        prefix = token[:3]
        if prefix not in RANK_BY_PREFIX:
            continue
        name = token[3:].strip()
        if not name:
            continue
        result.append((RANK_BY_PREFIX[prefix], name))
    return result


def write_taxdump(nodes_dmp, names_dmp, parent_by_taxid, rank_by_taxid, name_by_taxid):
    # type: (Path, Path, Dict[int, int], Dict[int, str], Dict[int, str]) -> None
    nodes_dmp.parent.mkdir(parents=True, exist_ok=True)

    # nodes.dmp format (NCBI style):
    # tax_id | parent_tax_id | rank | embl_code | division_id | inherited_div_flag | genetic_code_id |
    # inherited_GC_flag | mitochondrial_genetic_code_id | inherited_MGC_flag | GenBank_hidden_flag |
    # hidden_subtree_root_flag | comments |
    with nodes_dmp.open("w", encoding="utf-8") as out:
        for taxid in sorted(parent_by_taxid):
            parent = parent_by_taxid[taxid]
            rank = rank_by_taxid[taxid]
            out.write(
                f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            )

    # names.dmp format:
    # tax_id | name_txt | unique_name | name_class |
    with names_dmp.open("w", encoding="utf-8") as out:
        for taxid in sorted(name_by_taxid):
            name = name_by_taxid[taxid]
            out.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create an NCBI-style taxdump (nodes.dmp/names.dmp) from GTDB taxonomy TSVs, restricted to a representative-genome set."
    )
    parser.add_argument("--bac120-taxonomy", required=True, type=Path, help="Path to bac120_taxonomy_r226.tsv.gz")
    parser.add_argument("--ar53-taxonomy", required=True, type=Path, help="Path to ar53_taxonomy_r226.tsv.gz")
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--genomes-root", type=Path, help="Root directory containing GTDB rep genomes (*.fna[.gz])")
    source_group.add_argument("--tarball", type=Path, help="GTDB representative-genomes tarball")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory for nodes.dmp/names.dmp")
    parser.add_argument("--acc2taxid", required=True, type=Path, help="Output TSV mapping (accession<TAB>species_taxid)")
    args = parser.parse_args()

    rep_accessions = iter_rep_accessions(genomes_root=args.genomes_root, tarball=args.tarball)
    if not rep_accessions:
        source = args.genomes_root if args.genomes_root is not None else args.tarball
        raise SystemExit(f"No representative genomes found from: {source}")

    next_taxid = 2
    parent_by_taxid = {1: 1}  # type: Dict[int, int]
    rank_by_taxid = {1: "no rank"}  # type: Dict[int, str]
    name_by_taxid = {1: "root"}  # type: Dict[int, str]

    node_taxid_by_key = {}  # type: Dict[Tuple[int, str, str], int]
    acc_to_species_taxid = {}  # type: Dict[str, int]

    def get_or_create_node(parent_taxid: int, rank: str, name: str) -> int:
        nonlocal next_taxid
        key = (parent_taxid, rank, name)
        existing = node_taxid_by_key.get(key)
        if existing is not None:
            return existing
        taxid = next_taxid
        next_taxid += 1
        node_taxid_by_key[key] = taxid
        parent_by_taxid[taxid] = parent_taxid
        rank_by_taxid[taxid] = rank
        name_by_taxid[taxid] = name
        return taxid

    def ingest_taxonomy_tsv(path):  # type: (Path) -> Tuple[int, int]
        kept = 0
        seen = 0
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip("\n")
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 2:
                    continue
                raw_acc = fields[0]
                acc = normalize_genome_accession(raw_acc)
                acc_nover = strip_accession_version(acc)
                seen += 1
                if acc not in rep_accessions and acc_nover not in rep_accessions:
                    continue
                lineage_value = find_lineage_column(fields)
                if not lineage_value:
                    continue
                lineage = parse_lineage(lineage_value)
                if not lineage:
                    continue
                current = 1
                species_taxid = None
                for rank, name in lineage:
                    current = get_or_create_node(current, rank, name)
                    if rank == "species":
                        species_taxid = current
                if species_taxid is None:
                    # Fallback: use the lowest rank present.
                    species_taxid = current
                acc_to_species_taxid[acc] = species_taxid
                acc_to_species_taxid[acc_nover] = species_taxid
                kept += 1
        return kept, seen

    kept_bac, seen_bac = ingest_taxonomy_tsv(args.bac120_taxonomy)
    kept_ar, seen_ar = ingest_taxonomy_tsv(args.ar53_taxonomy)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_taxdump(
        nodes_dmp=args.outdir / "nodes.dmp",
        names_dmp=args.outdir / "names.dmp",
        parent_by_taxid=parent_by_taxid,
        rank_by_taxid=rank_by_taxid,
        name_by_taxid=name_by_taxid,
    )

    args.acc2taxid.parent.mkdir(parents=True, exist_ok=True)
    with args.acc2taxid.open("w", encoding="utf-8") as out:
        for acc in sorted(acc_to_species_taxid):
            out.write(f"{acc}\t{acc_to_species_taxid[acc]}\n")

    missing = 0
    for acc in rep_accessions:
        if acc in acc_to_species_taxid:
            continue
        missing += 1
    print(f"[make_gtdb_taxdump] rep_accessions={len(rep_accessions)} bac_kept={kept_bac}/{seen_bac} ar_kept={kept_ar}/{seen_ar} taxids={len(parent_by_taxid)} missing_rep_accessions={missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
