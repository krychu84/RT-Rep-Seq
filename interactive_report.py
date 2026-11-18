#!/usr/bin/env python3
"""interactive_report.py -- (re)build RT-Rep-Seq's interactive report.html for an
existing results/ folder (batches/batchNN/results, latest/, or a plain --results-dir).

Since v1.1.0 the pipeline writes this interactive report itself (step 5, report.html);
this script is for rebuilding it -- e.g. with a different slider start value -- or for
upgrading results folders made by v1.0.0. Read-only on everything except the output file.

The report (see report_html.py): read-count threshold slider; 2 x 96-well plate maps
(green >= threshold, red below); species table, LCP heatmap (re-clustered live) and LCP
histograms for the passing barcodes; TSV/PNG downloads; A4 print layout.

Usage:
  python interactive_report.py --results-dir /path/to/results \
      [--counts /path/to/demux_read_counts.tsv] [--threshold 10000] [--out report.html]

Barcodes are only shown if they were classified/profiled: results made by v1.0.0 in
one-shot mode only contain the barcodes that passed --min-reads (rerun the pipeline with
`--skip demux` to analyse the rest).
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
import report_html  # noqa: E402
from rt_rep_seq import REPO_URL, VERSION  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', required=True)
    ap.add_argument('--counts', default=None,
                     help='demux_read_counts.tsv (default: <results-dir>/demux_read_counts.tsv, '
                          'then ../work/demultiplexed_fastq/demux_read_counts.tsv)')
    ap.add_argument('--threshold', type=int, default=10000, help='slider start value (default: 10000)')
    ap.add_argument('--bin-bp', type=int, default=10, help='histogram bin width, bp (default: 10)')
    ap.add_argument('--min-read-length', type=int, default=500,
                     help='reads >= this length (bp) are used for the LCP analysis (default: 500)')
    ap.add_argument('--run-name', default=None)
    ap.add_argument('--gtdb-db', default=None, help='Metabuli DB used (for the Methods panel)')
    ap.add_argument('--metabuli-env', default=None, help='Metabuli conda env (for its version in Methods)')
    ap.add_argument('--out', default=None, help='output HTML (default: <results-dir>/report.html)')
    args = ap.parse_args()

    counts = args.counts
    if not counts:
        for c in (os.path.join(args.results_dir, 'demux_read_counts.tsv'),
                  os.path.join(args.results_dir, '..', 'work', 'demultiplexed_fastq', 'demux_read_counts.tsv')):
            if os.path.exists(c):
                counts = c
                break
    payload = report_html.build_payload(
        args.results_dir, counts_path=counts, default_threshold=args.threshold, binsize=args.bin_bp,
        min_len=args.min_read_length,
        meta={'run': args.run_name or os.path.basename(os.path.abspath(os.path.join(args.results_dir, '..'))),
              'repo_url': REPO_URL, 'version': VERSION,
              'methods': report_html.methods_meta(
                  onrep_config=os.path.join(HERE, 'assets', 'onrep_config.json'),
                  gtdb_db=args.gtdb_db, metabuli_env=args.metabuli_env,
                  tool_config=dict({'repo_url': REPO_URL}, **report_html.load_tool_config(HERE)), version=VERSION,
                  params={'min_reads': args.threshold, 'n_reference_barcodes': 192,
                          'min_read_length': args.min_read_length})},
        logo_ucph=os.path.join(HERE, 'assets', 'ucph_logo.png'),
        logo_onrepseq=os.path.join(HERE, 'assets', 'onrepseq_logo.png'))
    out = report_html.write_report(args.out or os.path.join(args.results_dir, 'report.html'), payload)
    print(f'[interactive_report] {len(payload["reads"])} barcodes -> {out}')


if __name__ == '__main__':
    main()
