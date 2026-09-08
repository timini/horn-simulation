#!/bin/bash
set -e
# The pinned image has OpenBLAS 0.3.20. Keep its x86 kernel consistent across
# hosted CPUs while guarding against intermittent P2 geometry/backend errors.
# Honour an explicit override for users validating another backend.
if [ "$(uname -m)" = "x86_64" ]; then
    export OPENBLAS_CORETYPE="${OPENBLAS_CORETYPE:-Nehalem}"
fi
source /usr/local/bin/dolfinx-complex-mode
exec "$@"
