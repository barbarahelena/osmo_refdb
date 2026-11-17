#!/usr/bin/env python3
"""
fetch_zenodo_member.py -- download ONE file out of a large zip on Zenodo,
without downloading the zip.

Cayman's prevalence-filtered GMGC catalogs are published as one 19.3 GB zip
(https://zenodo.org/records/10473258, gene_catalogues.zip) holding one
gzipped nucleotide FASTA per habitat, e.g.
  gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz  (2.99 GB)
The members are deflated inside the zip, but their content is already gzip,
so they are about as large as the .gz itself. This reads the zip's directory
with HTTP range requests, downloads the member's byte range resumably (to
<out>.part), inflates it locally and checks the member's CRC32.

Usage:
  python fetch_zenodo_member.py --list
  python fetch_zenodo_member.py --member gene_catalogues/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz \\
      --out /path/to/GMGC10.human-gut.95nr.0.5.percent.prevalence.fna.gz
  (--max-bytes N downloads only the first N compressed bytes, for format
  tests: the output is a truncated .gz and there is no CRC check)

Needs internet (a login or transfer node). Python >= 3.8, standard library only.
GMGC: Coelho et al. 2022, Nature 601; the prevalence-filtered sets: Ducarmon et
al., Nature Microbiology 2026 (Cayman). CC-BY-4.0.
"""
from __future__ import annotations

import argparse
import http.client
import io
import struct
import sys
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from pathlib import Path

DEFAULT_URL = "https://zenodo.org/api/records/10473258/files/gene_catalogues.zip/content"
CHUNK = 16 << 20


def request(url: str, start: int, end: int, tries: int = 6) -> tuple[bytes, str]:
    """bytes [start, end] of url (inclusive), and the final URL after redirects."""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.status != 206:
                    sys.exit(f"ERROR: the server ignored the Range header (HTTP {r.status}); "
                             "cannot fetch a single member")
                return r.read(), r.geturl()
        except (urllib.error.URLError, TimeoutError, ConnectionError,
                http.client.HTTPException) as exc:     # incl. IncompleteRead: cut-short chunk
            if attempt == tries - 1:
                raise
            wait = 2 ** attempt
            print(f"  retry in {wait}s ({exc})", file=sys.stderr)
            time.sleep(wait)
    raise AssertionError


class Remote(io.RawIOBase):
    """Read-only, seekable file over HTTP range requests (for zipfile)."""

    def __init__(self, url: str):
        req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            if r.status != 206:
                sys.exit(f"ERROR: the server ignored the Range header (HTTP {r.status}); "
                         "cannot fetch a single member")
            self.size = int(r.headers["Content-Range"].split("/")[1])
            self.url = r.geturl()          # the redirect target answers the later ranges
        self.pos = 0

    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def read(self, n=-1):
        n = self.size - self.pos if n < 0 else min(n, self.size - self.pos)
        if n <= 0:
            return b""
        data, _ = request(self.url, self.pos, self.pos + n - 1)
        self.pos += len(data)
        return data

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def crc32_of(path: Path, nbytes: int) -> int:
    crc = 0
    with open(path, "rb") as fh:
        left = nbytes
        while left:
            block = fh.read(min(left, 1 << 24))
            if not block:
                break
            crc = zlib.crc32(block, crc)
            left -= len(block)
    return crc & 0xFFFFFFFF


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", default=DEFAULT_URL)
    ap.add_argument("--list", action="store_true", help="list the zip's members and exit")
    ap.add_argument("--member", help="path of the member inside the zip")
    ap.add_argument("--out", type=Path, help="where to write it")
    ap.add_argument("--max-bytes", type=int, default=None,
                    help="only the first N bytes (format tests; no CRC check)")
    args = ap.parse_args()

    remote = Remote(args.url)
    zf = zipfile.ZipFile(io.BufferedReader(remote, 1 << 20))
    if args.list or not args.member:
        for i in zf.infolist():
            print(f"{i.file_size / 1e9:8.3f} GB  {'stored' if i.compress_type == 0 else 'deflated'}  {i.filename}")
        return
    if not args.out:
        sys.exit("ERROR: --out is required with --member")
    info = zf.getinfo(args.member)
    if info.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
        sys.exit(f"ERROR: {args.member}: unsupported zip compression method {info.compress_type}")

    # the member's data starts after its local header
    remote.seek(info.header_offset)
    sig, *_, nlen, xlen = struct.unpack("<4s5H3L2H", remote.read(30))
    if sig != b"PK\x03\x04":
        sys.exit("ERROR: bad local file header")
    start = info.header_offset + 30 + nlen + xlen
    total = info.compress_size if args.max_bytes is None else min(args.max_bytes, info.compress_size)

    part = args.out.with_name(args.out.name + ".part")
    part.parent.mkdir(parents=True, exist_ok=True)
    have = part.stat().st_size if part.exists() else 0
    if have > total:
        part.unlink()
        have = 0
    print(f"{args.member}: {total / 1e9:.3f} GB, resuming at {have / 1e9:.3f} GB", flush=True)
    t0 = time.time()
    with open(part, "ab") as fh:
        while have < total:
            end = min(have + CHUNK, total) - 1
            data, _ = request(remote.url, start + have, start + end)
            fh.write(data)
            have += len(data)
            print(f"  {have / 1e9:7.3f} / {total / 1e9:.3f} GB", flush=True)

    if info.compress_type == zipfile.ZIP_STORED:
        if args.max_bytes is None:
            got = crc32_of(part, total)
            if got != info.CRC:
                sys.exit(f"ERROR: CRC32 mismatch for {part} ({got:08x} != {info.CRC:08x}); "
                         "delete it and retry")
            print(f"CRC32 ok ({info.CRC:08x})")
        part.replace(args.out)
    else:
        # raw deflate -> the member's own bytes (here a .gz), streamed
        inflater = zlib.decompressobj(-zlib.MAX_WBITS)
        crc = n_out = 0
        tmp = args.out.with_name(args.out.name + ".inflating")
        with open(part, "rb") as src, open(tmp, "wb") as dst:
            while True:
                block = src.read(1 << 24)
                if not block:
                    break
                try:
                    out = inflater.decompress(block)
                except zlib.error as exc:
                    sys.exit(f"ERROR: inflating {part}: {exc}; delete it and retry")
                dst.write(out)
                crc = zlib.crc32(out, crc)
                n_out += len(out)
            out = inflater.flush()
            dst.write(out)
            crc = zlib.crc32(out, crc)
            n_out += len(out)
        if args.max_bytes is None:
            if (crc & 0xFFFFFFFF) != info.CRC or n_out != info.file_size:
                sys.exit(f"ERROR: CRC32/size mismatch for {args.member}; delete {part} and retry")
            print(f"CRC32 ok ({info.CRC:08x}), {n_out / 1e9:.3f} GB")
        tmp.replace(args.out)
        part.unlink()
    print(f"written {args.out} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
