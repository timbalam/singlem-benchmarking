#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

from gtdb_utils import (
    accession_from_genome_name,
    iter_genome_members_from_tarball,
    iter_genome_paths_from_root,
    normalize_genome_accession,
    open_text_from_member,
    strip_accession_version,
)


def load_acc2taxid(path: Path) -> Dict[str, str]:
    mapping = {}  # type: Dict[str, str]
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            acc, taxid = line.split("\t", 1)
            mapping[acc] = taxid
    return mapping


def parse_header(raw_line: bytes) -> Tuple[str, bytes]:
    line = raw_line.rstrip(b"\r\n")
    payload = line[1:]
    parts = payload.split(None, 1)
    seqid = parts[0].decode("utf-8", errors="replace")
    rest = b""
    if len(parts) > 1:
        rest = b" " + parts[1]
    return seqid, rest


class ChunkWriter:
    def __init__(self, outdir: Path, chunk_size_bytes: int):
        self.outdir = outdir
        self.chunk_size_bytes = chunk_size_bytes
        self.chunk_index = 0
        self.current_size = 0
        self.fasta_handle = None
        self.prelim_handle = None
        self.chunks_written = 0

    def _open_chunk(self) -> None:
        self.chunk_index += 1
        fasta_path = self.outdir / f"library_{self.chunk_index:05d}.fna"
        prelim_path = self.outdir / f"prelim_map_{self.chunk_index:05d}.txt"
        self.fasta_handle = fasta_path.open("wb")
        self.prelim_handle = prelim_path.open("w", encoding="utf-8")
        self.current_size = 0
        self.chunks_written += 1

    def ensure_chunk(self) -> None:
        if self.fasta_handle is None:
            self._open_chunk()

    def rotate_if_needed(self) -> None:
        if self.fasta_handle is None:
            self._open_chunk()
            return
        if self.current_size < self.chunk_size_bytes:
            return
        self.close()
        self._open_chunk()

    def write_header(self, seqid: str, taxid: str, rest: bytes) -> None:
        self.ensure_chunk()
        rewritten_seqid = f"{seqid}|kraken:taxid|{taxid}"
        header = b">" + rewritten_seqid.encode("utf-8") + rest + b"\n"
        self.fasta_handle.write(header)
        self.prelim_handle.write(f"TAXID\t{rewritten_seqid}\t{taxid}\n")
        self.current_size += len(header)

    def write_sequence_line(self, raw_line: bytes) -> None:
        self.ensure_chunk()
        self.fasta_handle.write(raw_line)
        self.current_size += len(raw_line)

    def close(self) -> None:
        if self.fasta_handle is not None:
            self.fasta_handle.close()
            self.fasta_handle = None
        if self.prelim_handle is not None:
            self.prelim_handle.close()
            self.prelim_handle = None


def resolve_taxid(mapping: Dict[str, str], accession: str) -> Optional[str]:
    for candidate in (accession, normalize_genome_accession(accession), strip_accession_version(accession)):
        taxid = mapping.get(candidate)
        if taxid is not None:
            return taxid
    return None


def write_genome_stream(writer: ChunkWriter, stream, accession: str, taxid: str) -> int:
    seq_count = 0
    for raw_line in stream:
        if raw_line.startswith(b">"):
            seq_count += 1
            seqid, rest = parse_header(raw_line)
            writer.write_header(seqid, taxid, rest)
        else:
            writer.write_sequence_line(raw_line)
    return seq_count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create Kraken2 library FASTA chunks and prelim_map files from GTDB representative genomes."
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--genomes-root", type=Path, help="Root directory containing GTDB rep genomes (*.fna[.gz])")
    source_group.add_argument("--tarball", type=Path, help="GTDB representative-genomes tarball")
    parser.add_argument("--acc2taxid", required=True, type=Path, help="TSV mapping accession<TAB>taxid")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory, typically DB/library/added")
    parser.add_argument(
        "--chunk-size-bytes",
        type=int,
        default=20000000000,
        help="Rotate to a new library chunk after approximately this many output bytes (default: 20,000,000,000)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        help="Emit a progress line after this many genomes are written (default: 1000)",
    )
    args = parser.parse_args()

    mapping = load_acc2taxid(args.acc2taxid)
    if not mapping:
        raise SystemExit(f"Empty mapping file: {args.acc2taxid}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    writer = ChunkWriter(args.outdir, args.chunk_size_bytes)

    genomes_total = 0
    genomes_written = 0
    sequences_written = 0
    missing_taxid = 0

    if args.tarball is not None:
        iterator = iter_genome_members_from_tarball(args.tarball)
        for member_name, member_handle in iterator:
            genomes_total += 1
            acc = accession_from_genome_name(member_name)
            if not acc:
                continue
            taxid = resolve_taxid(mapping, acc)
            if taxid is None:
                missing_taxid += 1
                continue
            writer.rotate_if_needed()
            stream = open_text_from_member(member_handle, member_name)
            try:
                sequences_written += write_genome_stream(writer, stream, acc, taxid)
            finally:
                stream.close()
                member_handle.close()
            genomes_written += 1
            if args.progress_every > 0 and genomes_written % args.progress_every == 0:
                print(
                    f"[materialize_kraken_library] progress genomes_written={genomes_written} "
                    f"sequences_written={sequences_written} chunks={writer.chunks_written}",
                    flush=True,
                )
    else:
        for genome_path in iter_genome_paths_from_root(args.genomes_root):
            genomes_total += 1
            acc = accession_from_genome_name(genome_path.name)
            if not acc:
                continue
            taxid = resolve_taxid(mapping, acc)
            if taxid is None:
                missing_taxid += 1
                continue
            writer.rotate_if_needed()
            if genome_path.suffix == ".gz":
                import gzip

                handle = gzip.open(str(genome_path), "rb")
            else:
                handle = genome_path.open("rb")
            try:
                sequences_written += write_genome_stream(writer, handle, acc, taxid)
            finally:
                handle.close()
            genomes_written += 1
            if args.progress_every > 0 and genomes_written % args.progress_every == 0:
                print(
                    f"[materialize_kraken_library] progress genomes_written={genomes_written} "
                    f"sequences_written={sequences_written} chunks={writer.chunks_written}",
                    flush=True,
                )

    writer.close()
    print(
        f"[materialize_kraken_library] genomes_total={genomes_total} genomes_written={genomes_written} "
        f"sequences_written={sequences_written} missing_taxid={missing_taxid} chunks={writer.chunks_written}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
