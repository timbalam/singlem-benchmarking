#!/usr/bin/env python3
"""Build a Kraken2-style ktaxonomy.tsv from a Kraken2 taxonomy/ directory
(names.dmp + nodes.dmp).

Output columns (pipe-separated, tab-around the pipes, to match the sink DB
format expected by build_validation_case.py / kraken_taxonomy_utils.py):

    taxid \t | \t parent_taxid \t | \t rank_code \t | \t depth \t | \t name

Rank codes follow the same letter scheme used by the sink ktaxonomy:
  superkingdom/domain -> D, phylum -> P, class -> C, order -> O,
  family -> F, genus -> G, species -> S, no rank -> 'R' for root else '-',
  subspecies/strain -> S1.
"""

import argparse
from typing import Dict, List, Tuple


RANK_MAP = {
    "superkingdom": "D",
    "domain": "D",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
    "subspecies": "S1",
    "strain": "S1",
}


def parse_dmp_line(line: str) -> List[str]:
    return [c.strip() for c in line.rstrip("\n").rstrip("\t|").split("\t|\t")]


def build(taxdir: str, out_path: str) -> None:
    # names: taxid -> scientific name
    names: Dict[int, str] = {}
    with open(f"{taxdir}/names.dmp") as f:
        for line in f:
            parts = parse_dmp_line(line)
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[int(parts[0])] = parts[1]

    # nodes: taxid -> (parent, rank_string)
    nodes: Dict[int, Tuple[int, str]] = {}
    with open(f"{taxdir}/nodes.dmp") as f:
        for line in f:
            parts = parse_dmp_line(line)
            if len(parts) >= 3:
                tid = int(parts[0])
                par = int(parts[1])
                rk = parts[2]
                nodes[tid] = (par, rk)

    # depth = distance from root (taxid 1)
    depth_cache: Dict[int, int] = {1: 0}

    def depth_of(tid: int) -> int:
        if tid in depth_cache:
            return depth_cache[tid]
        seen = []
        cur = tid
        while cur not in depth_cache:
            seen.append(cur)
            par = nodes.get(cur, (1, ""))[0]
            if par == cur:
                depth_cache[cur] = 0
                break
            cur = par
        base = depth_cache[cur]
        for i, t in enumerate(reversed(seen)):
            depth_cache[t] = base + i + 1
        return depth_cache[tid]

    with open(out_path, "w") as out:
        for tid in sorted(nodes.keys()):
            par, rk = nodes[tid]
            name = names.get(tid, "")
            if tid == 1:
                code = "R"
            else:
                code = RANK_MAP.get(rk, "-")
            d = depth_of(tid)
            out.write(f"{tid}\t|\t{par}\t|\t{code}\t|\t{d}\t|\t{name}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build ktaxonomy.tsv from a Kraken2 taxonomy directory containing nodes.dmp and names.dmp."
    )
    parser.add_argument("taxdir", help="Kraken2 taxonomy directory")
    parser.add_argument("out_path", help="Output ktaxonomy.tsv path")
    args = parser.parse_args()
    build(args.taxdir, args.out_path)
    print(f"wrote: {args.out_path}")
