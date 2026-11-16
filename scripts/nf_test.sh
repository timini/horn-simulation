#!/usr/bin/env bash
# Run the pinned nf-test release. Uses $NF_TEST when set, otherwise a cached
# copy under .tools/, downloading and checksum-verifying it on first use.
set -euo pipefail

VERSION=0.9.5
SHA256=b7679eb90cdc9642bfa89e9634db02ce6e699de53ad35ef0e2f4847634fc1641
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$ROOT/.tools/nf-test-$VERSION"

# Select Java 17 to 22 the same way the launcher does, so nf-test and the
# Nextflow it starts work without JAVA_HOME set in the shell.
if java_home=$(python3 -c 'import sys; sys.path.insert(0, sys.argv[1]); import run_pipeline; print(run_pipeline.java_environment().get("NXF_JAVA_HOME", ""))' "$ROOT/scripts" 2> /dev/null) \
        && [[ -n "$java_home" ]]; then
    export JAVA_HOME="$java_home" NXF_JAVA_HOME="$java_home"
    unset JAVA_CMD
fi

if [[ -n "${NF_TEST:-}" ]]; then
    exec "$NF_TEST" "$@"
fi

if [[ ! -x "$DEST/nf-test" ]]; then
    mkdir -p "$DEST"
    archive="$DEST/nf-test-$VERSION.tar.gz"
    curl -fsSL -o "$archive" \
        "https://github.com/askimed/nf-test/releases/download/v$VERSION/nf-test-$VERSION.tar.gz"
    if command -v sha256sum > /dev/null; then
        actual=$(sha256sum "$archive" | cut -d' ' -f1)
    else
        actual=$(shasum -a 256 "$archive" | cut -d' ' -f1)
    fi
    if [[ "$actual" != "$SHA256" ]]; then
        echo "nf-test $VERSION checksum mismatch: $actual" >&2
        rm -f "$archive"
        exit 1
    fi
    tar -xzf "$archive" -C "$DEST"
    rm -f "$archive"
    chmod +x "$DEST/nf-test"
fi

reported=$("$DEST/nf-test" version 2>&1 | grep -Eo '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)
if [[ "$reported" != "$VERSION" ]]; then
    echo "Expected nf-test $VERSION in $DEST, found '$reported'" >&2
    exit 1
fi
exec "$DEST/nf-test" "$@"
