#!/usr/bin/env python3

import gzip
import re
import tarfile
from pathlib import Path
from typing import Generator, IO, Iterable, Optional, Set, Tuple


_GENOME_SUFFIX_RE = re.compile(r"_genomic\.fna(\.gz)?$")


def normalize_genome_accession(raw: str) -> str:
    value = raw.strip()
    if value.startswith("GB_") or value.startswith("RS_"):
        value = value[3:]
    return value


def strip_accession_version(accession: str) -> str:
    return re.sub(r"\.\d+$", "", accession)


def accession_from_genome_name(name: str) -> Optional[str]:
    basename = Path(name).name
    if not _GENOME_SUFFIX_RE.search(basename):
        return None
    return normalize_genome_accession(_GENOME_SUFFIX_RE.sub("", basename))


def accession_from_genome_path(path: Path) -> Optional[str]:
    return accession_from_genome_name(path.name)


def iter_rep_accessions_from_genomes_root(genomes_root: Path) -> Set[str]:
    reps = set()  # type: Set[str]
    for pattern in ("*_genomic.fna.gz", "*_genomic.fna"):
        for genome_path in genomes_root.rglob(pattern):
            acc = accession_from_genome_path(genome_path)
            if not acc:
                continue
            reps.add(acc)
            reps.add(strip_accession_version(acc))
    return reps


def iter_rep_accessions_from_tarball(tarball: Path) -> Set[str]:
    reps = set()  # type: Set[str]
    with tarfile.open(str(tarball), "r|*") as handle:
        for member in handle:
            if not member.isfile():
                continue
            acc = accession_from_genome_name(member.name)
            if not acc:
                continue
            reps.add(acc)
            reps.add(strip_accession_version(acc))
    return reps


def iter_rep_accessions(genomes_root: Optional[Path] = None, tarball: Optional[Path] = None) -> Set[str]:
    if genomes_root is not None:
        return iter_rep_accessions_from_genomes_root(genomes_root)
    if tarball is not None:
        return iter_rep_accessions_from_tarball(tarball)
    raise ValueError("Either genomes_root or tarball must be provided")


def open_text_from_member(member_handle: IO[bytes], member_name: str):
    if member_name.endswith(".gz"):
        gz = gzip.GzipFile(fileobj=member_handle)
        return gz
    return member_handle


def iter_genome_members_from_tarball(tarball: Path) -> Generator[Tuple[str, IO[bytes]], None, None]:
    with tarfile.open(str(tarball), "r|*") as handle:
        for member in handle:
            if not member.isfile():
                continue
            acc = accession_from_genome_name(member.name)
            if not acc:
                continue
            member_handle = handle.extractfile(member)
            if member_handle is None:
                continue
            yield member.name, member_handle


def iter_genome_paths_from_root(genomes_root: Path) -> Iterable[Path]:
    for pattern in ("*_genomic.fna.gz", "*_genomic.fna"):
        for genome_path in sorted(genomes_root.rglob(pattern)):
            yield genome_path
