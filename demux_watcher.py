#!/usr/bin/env python3
"""demux_watcher.py -- incrementally demultiplex new fastq(.gz) files as they appear.

Watches --raw-dir for new *.fastq(.gz) files (e.g. a live MinKNOW output folder) and,
for each one, demultiplexes it and APPENDS the classified reads to --demux-dir. Already
-processed files are never reprocessed, and existing barcode data is never wiped or
overwritten -- this folder only grows.

Name --demux-dir `demultiplexed_fastq` (its basename matters) so that rt_rep_seq.py's
`--skip demux` mode can read it directly as `<work-dir>/demultiplexed_fastq/`. Pairs with
analysis_watcher.py, which watches --demux-dir for growth and reruns classify/LCPs/
heatmap/report on top of it whenever it changes.

Prints the same startup banner and resource-budget line as rt_rep_seq.py, and is
pinned to the same CPU/RAM budget by default (half the machine) -- see --cpu-frac/
--mem-frac/--no-banner/--no-mem-cap.

Usage:
  python demux_watcher.py --raw-dir /path/to/live_fastq_pass \\
      --demux-dir /path/to/demultiplexed_fastq [--poll-interval 5] [--min-reads 10000]

Stop gracefully by creating a file named .demux_watcher_stop inside --demux-dir.
"""
import argparse
import glob
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import rt_rep_seq as rt  # noqa: E402  (shared Resources/banner/resource_snapshot helpers)

PROCESSED_LOG = '.processed_files.txt'
STOP_MARKER = '.demux_watcher_stop'


def load_processed(demux_dir):
    path = os.path.join(demux_dir, PROCESSED_LOG)
    if not os.path.exists(path):
        return set()
    return set(open(path).read().splitlines())


def mark_processed(demux_dir, filename):
    with open(os.path.join(demux_dir, PROCESSED_LOG), 'a') as f:
        f.write(filename + '\n')


def rebuild_counts(demux_dir):
    """Recount every {barcode}.fastq already in demux_dir -- lets the watcher resume
    cleanly after a restart without trusting a possibly-stale counts file."""
    counts = Counter()
    for fq_path in glob.glob(os.path.join(demux_dir, 'repBC*.fastq')):
        bc = os.path.splitext(os.path.basename(fq_path))[0]
        with open(fq_path) as f:
            counts[bc] = sum(1 for _ in f) // 4
    return counts


def ready_files(paths, state, settle):
    """Files that are safe to read: non-empty, and size + mtime unchanged since the previous
    poll AND not modified for >= `settle` seconds. A file still being written (MinKNOW, cp,
    rsync) keeps changing and is skipped until it settles -- reading it early would demux a
    truncated file and, since files are only ever processed once, lose the rest of it."""
    now, ready = time.time(), []
    for p in paths:
        try:
            st = os.stat(p)
        except FileNotFoundError:
            state.pop(p, None)
            continue
        sig = (st.st_size, st.st_mtime)
        prev = state.get(p)
        state[p] = sig
        if st.st_size > 0 and prev == sig and now - st.st_mtime >= settle:
            ready.append(p)
    return ready


def write_counts(demux_dir, counts, min_reads):
    path = os.path.join(demux_dir, 'demux_read_counts.tsv')
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        f.write('barcode\treads\tpassed\n')
        for bc, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            f.write(f"{bc}\t{n}\t{'yes' if n >= min_reads else 'no'}\n")
    os.replace(tmp, path)  # atomic -- analysis_watcher never sees a half-written file


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--raw-dir', required=True, help='directory receiving new *.fastq(.gz) files')
    ap.add_argument('--demux-dir', required=True,
                     help="persistent output dir -- name it 'demultiplexed_fastq' so "
                          'rt_rep_seq.py --skip demux can read it directly')
    ap.add_argument('--onrep-barcodes',
                     default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                           'assets', 'onrep_barcodes.fasta'))
    ap.add_argument('--onrep-config',
                     default=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                           'assets', 'onrep_config.json'))
    ap.add_argument('--min-reads', type=int, default=10000,
                     help='threshold recorded in demux_read_counts.tsv\'s "passed" column -- '
                          'bookkeeping only, this watcher never stops early (default: 10000)')
    ap.add_argument('--poll-interval', type=int, default=5,
                     help='seconds between checks for new raw files (default: 5)')
    ap.add_argument('--settle', type=int, default=10,
                     help='only read a file once it has not changed for this many seconds '
                          '(default: 10) -- guards against half-written files')
    ap.add_argument('--max-retries', type=int, default=3,
                     help='a file that fails to read or yields 0 reads is retried this many '
                          'times before being skipped for good (default: 3)')
    ap.add_argument('--cpu-frac', type=float, default=0.5,
                     help='share of CPU cores this process is pinned to (default: 0.5)')
    ap.add_argument('--mem-frac', type=float, default=0.5,
                     help='share of system RAM this process is capped to via a '
                          'systemd-run user scope (default: 0.5)')
    ap.add_argument('--no-mem-cap', action='store_true')
    ap.add_argument('--no-banner', action='store_true')
    args = ap.parse_args()

    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass

    rt.TOOL_RESOURCES = rt.Resources(args.cpu_frac, args.mem_frac)
    if not args.no_mem_cap:
        rt.TOOL_RESOURCES.reexec_with_mem_cap(sys.argv)
    rt.TOOL_RESOURCES.apply_cpu_cap()
    if not args.no_banner:
        rt.banner()
    cap_note = '' if (args.no_mem_cap or rt.TOOL_RESOURCES.mem_capped) else ' (RAM cap not enforced)'
    print(f'[demux-watcher] resource budget: CPU {len(rt.TOOL_RESOURCES.cpus)}/{os.cpu_count()} '
          f'cores (cap {args.cpu_frac * 100:.0f}%), RAM cap '
          f'{rt.human_bytes(rt.TOOL_RESOURCES.mem_limit)}{cap_note}', flush=True)

    from torchlex.config import load_config
    from torchlex.demultiplexer import Demultiplexer
    from torchlex.fastq_reader import FastqFile

    config = load_config(args.onrep_config)
    config['amplicones']['barcodes_path'] = args.onrep_barcodes
    demux = Demultiplexer(config)

    os.makedirs(args.demux_dir, exist_ok=True)
    processed = load_processed(args.demux_dir)
    counts = rebuild_counts(args.demux_dir)
    if counts:
        write_counts(args.demux_dir, counts, args.min_reads)

    print(f'[demux-watcher] observing {args.raw_dir} for new *.fastq(.gz) -- checking every '
          f'{args.poll_interval}s -- {len(processed)} file(s) already processed -- '
          f'output -> {args.demux_dir}', flush=True)

    file_state, failures = {}, Counter()
    rt.status_init(os.path.join(os.path.dirname(os.path.abspath(args.demux_dir)), '.status', 'demux.json'),
                   state='waiting', raw_dir=os.path.abspath(args.raw_dir), min_reads=args.min_reads, pid=os.getpid())

    def report_status(state, current=None, n_present=None):
        rt.status_update(force=True, state=state, current=current, files_present=n_present,
                         files_done=len(processed), reads=sum(counts.values()), barcodes=len(counts),
                         n_pass=sum(1 for c in counts.values() if c >= args.min_reads))
    while True:
        if os.path.exists(os.path.join(args.demux_dir, STOP_MARKER)):
            print('[demux-watcher] stop marker found -- exiting', flush=True)
            break

        candidates = sorted(glob.glob(os.path.join(args.raw_dir, '*.fastq')) +
                             glob.glob(os.path.join(args.raw_dir, '*.fastq.gz')))
        new_files = [p for p in candidates if os.path.basename(p) not in processed]
        pending_n = len(new_files)
        new_files = ready_files(new_files, file_state, args.settle)

        if not new_files:
            report_status('waiting for new files' if not pending_n else 'waiting for files to finish writing',
                          n_present=len(candidates))
            time.sleep(args.poll_interval)
            continue

        for path in new_files:
            name = os.path.basename(path)
            report_status('demultiplexing', current=name, n_present=len(candidates))
            try:
                fq = FastqFile(path)   # reads the whole file before anything is written out
            except Exception as e:     # truncated gzip etc. -- nothing appended yet, safe to retry
                failures[name] += 1
                print(f'[demux-watcher] {name}: could not read ({e.__class__.__name__}: {e}) -- '
                      f'attempt {failures[name]}/{args.max_retries}', flush=True)
                if failures[name] >= args.max_retries:
                    print(f'[demux-watcher] {name}: giving up, skipped', flush=True)
                    mark_processed(args.demux_dir, name)
                    processed.add(name)
                file_state.pop(path, None)
                continue
            if len(fq) == 0 and failures[name] + 1 < args.max_retries:
                failures[name] += 1
                print(f'[demux-watcher] {name}: 0 reads -- probably still being written, will retry '
                      f'(attempt {failures[name]}/{args.max_retries})', flush=True)
                file_state.pop(path, None)
                continue
            classified = demux.demultiplex_file(fq, args.demux_dir)  # appends {barcode}.fastq
            for row in classified.rows:
                counts[row.barcode] += 1
            write_counts(args.demux_dir, counts, args.min_reads)
            mark_processed(args.demux_dir, os.path.basename(path))
            processed.add(os.path.basename(path))
            n_ge = sum(1 for c in counts.values() if c >= args.min_reads)
            print(f'[demux-watcher] {os.path.basename(path)}: {len(fq)} reads, '
                  f'{len(classified.rows)} classified | {len(counts)} barcodes seen so far, '
                  f'{n_ge} >= {args.min_reads} reads  ({rt.resource_snapshot()})', flush=True)


if __name__ == '__main__':
    main()
