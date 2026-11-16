#!/usr/bin/env python3
"""Hash the source files each Docker image is built from.

`just build` stores the hash as the image label `horn.source_sha256`; the
launcher recomputes it and warns when an image was built from other sources.
Usage: image_source.py IMAGE  (prints the hash)
"""
import hashlib
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
LABEL = "horn.source_sha256"
COMMON = ("pyproject.toml", "uv.lock")
# Paths copied into each image's production stage, as in its Dockerfile.
IMAGE_SOURCES = {
    "horn-analysis": ("packages/horn-core", "packages/horn-drivers", "packages/horn-analysis"),
    "horn-geometry": ("packages/horn-core", "packages/horn-geometry"),
    "horn-solver": ("packages",),
}


def image_source_files(image, root=ROOT):
    paths = (*COMMON, f"packages/{image}/Dockerfile", *IMAGE_SOURCES[image])
    listed = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z", "--", *paths],
        cwd=root).decode().split("\0")
    return sorted({name for name in listed if name and (root/name).is_file()})


def image_source_hash(image, root=ROOT):
    """One sha256 over the relative path and content of every source file."""
    combined = hashlib.sha256()
    for name in image_source_files(image, root):
        combined.update(name.encode() + b"\0" + hashlib.sha256((root/name).read_bytes()).digest())
    return combined.hexdigest()


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in IMAGE_SOURCES:
        sys.exit(f"usage: {sys.argv[0]} {{{','.join(IMAGE_SOURCES)}}}")
    print(image_source_hash(sys.argv[1]))
