#!/usr/bin/env bash
# RT-Rep-Seq test run on the bundled simulated data (test/data, 12 barcodes, ~13,000 reads).
#
# It imitates a sequencing run: the six fastq.gz files are dropped one by one (every 20 s) into a
# MinKNOW-style folder while RT-Rep-Seq is watching it. Takes ~5-15 min, mostly Metabuli loading
# its database. Expected result: see test/README.md.
#
#   test/run_test.sh              live test (default)
#   test/run_test.sh --all-now    copy all files at once (single analysis round)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN="${RT_REP_SEQ_TEST_DIR:-$HOME/RT-Rep-Seq_test_run}"
FASTQ="$RUN/SIM_test/no_sample/20261118_1930_SIM_SIM00001_rtrepseqtest/fastq_pass"
DELAY=20; [[ "${1:-}" == "--all-now" ]] && DELAY=0
command -v RT-Rep-Seq >/dev/null || { echo "RT-Rep-Seq not found -- run ./install.sh first (and open a new terminal)"; exit 1; }

if [[ -d "$RUN" ]]; then echo "$RUN exists -- remove it first for a fresh test: rm -rf $RUN"; exit 1; fi
mkdir -p "$FASTQ"
files=("$HERE"/data/*.fastq.gz)
cp "${files[0]}" "$FASTQ/"
# feed the remaining files in the background, like MinKNOW writing new chunks
( for f in "${files[@]:1}"; do sleep "$DELAY"; cp "$f" "$FASTQ/.tmp_$(basename "$f")"; mv "$FASTQ/.tmp_$(basename "$f")" "$FASTQ/$(basename "$f")"; done ) &
cd "$FASTQ"
echo "Simulated MinKNOW folder: $FASTQ  (${#files[@]} files arrive every ${DELAY}s)"
exec RT-Rep-Seq --min-reads 500
