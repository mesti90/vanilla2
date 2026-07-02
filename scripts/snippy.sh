#!/usr/bin/env bash
set -euo pipefail

OUTDIR="$1"
TMPDIR="$2"
REF="$3"
R1="$4"
R2="$5"
THREADS="$6"

rm -rf "$OUTDIR" "$TMPDIR"

snippy \
	--outdir "$OUTDIR" \
	--tmpdir "$TMPDIR" \
	--ref "$REF" \
	--R1 "$R1" \
	--R2 "$R2" \
	--cpus "$THREADS" \
	--force \
	--quiet