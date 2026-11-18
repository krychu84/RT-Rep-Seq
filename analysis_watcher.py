#!/usr/bin/env python3
"""analysis_watcher.py -- rerun classify/LCPs/heatmap/report whenever a persistent
demultiplexed_fastq/ folder (as filled in by demux_watcher.py) grows, saving each
completed cycle as its own never-overwritten batchNN folder, and optionally also
mirroring the newest cycle's results into a single fixed, overwritten --latest-dir.

Never runs two cycles concurrently: this loop is single-threaded, so it simply doesn't
check for new data again until the previous cycle's rt_rep_seq run_pipeline() call
returns. If more data arrived while a cycle was running, that's noticed immediately on
the next check (no extra sleep), so cycles run back-to-back for as long as data keeps
growing, and idle (just a cheap poll) once it stops.

Each cycle re-tags and reclassifies EVERY barcode currently in --demux-dir (not just
what's new) -- that folder only grows, so each batch is a fresh, complete analysis of
everything demultiplexed so far, not a delta.

Usage:
  python analysis_watcher.py --demux-dir /path/to/demultiplexed_fastq \\
      --batches-dir /path/to/batches --gtdb-db /path/to/gtdb_db --metabuli-env metabuli \\
      [--latest-dir /path/to/results] [--poll-interval 10] [--min-reads 1]

Stop gracefully by creating a file named .analysis_watcher_stop inside --demux-dir
(checked between cycles, never mid-cycle).
"""
import argparse
import glob
import os
import re
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import rt_rep_seq as rt  # noqa: E402

STOP_MARKER = '.analysis_watcher_stop'
BATCH_RE = re.compile(r'^batch(\d+)$')


def demux_state(demux_dir):
    """(barcode, byte size) snapshot of every per-barcode fastq -- growth check that's
    cheap even once the folder holds hundreds of MB, since it only stats file sizes."""
    return tuple(sorted(
        (os.path.basename(p), os.path.getsize(p))
        for p in glob.glob(os.path.join(demux_dir, 'repBC*.fastq'))))


def next_batch_dir(batches_dir):
    nums = [int(m.group(1)) for p in glob.glob(os.path.join(batches_dir, 'batch*'))
            for m in [BATCH_RE.match(os.path.basename(p))] if m]
    n = (max(nums) + 1) if nums else 1
    return os.path.join(batches_dir, f'batch{n:02d}')


def build_cycle_args(args, cycle_work_dir, cycle_results_dir):
    """A rt_rep_seq argparse.Namespace for one classify+lcps+heatmap+report cycle that
    reads demux data straight from args.demux_dir via --skip demux."""
    argv = [
        'analysis_watcher',
        '--input-dir', args.demux_dir,  # unused (demux skipped) but argparse requires it
        '--gtdb-db', args.gtdb_db,
        '--work-dir', cycle_work_dir,
        '--results-dir', cycle_results_dir,
        '--min-reads', str(args.qc_min_reads),
        '--analysis-min-reads', str(args.min_reads),
        '--skip', 'demux',
    ]
    if args.metabuli_env:
        argv += ['--metabuli-env', args.metabuli_env]
    if args.metabuli_bin:
        argv += ['--metabuli-bin', args.metabuli_bin]
    if args.threads:
        argv += ['--threads', str(args.threads)]
    if args.run_name:
        argv += ['--run-name', args.run_name]
    old_argv = sys.argv
    sys.argv = argv
    try:
        return rt.parse_args()
    finally:
        sys.argv = old_argv


def run_one_cycle(args, batch_dir):
    cycle_work_dir = os.path.join(batch_dir, '_work')
    cycle_results_dir = os.path.join(batch_dir, 'results')
    os.makedirs(cycle_work_dir, exist_ok=True)
    # symlink (not copy) the growing, potentially large demux folder into the place
    # rt_rep_seq's hardcoded `<work-dir>/demultiplexed_fastq/...` paths expect it
    os.symlink(os.path.abspath(args.demux_dir), os.path.join(cycle_work_dir, 'demultiplexed_fastq'))

    cycle_args = build_cycle_args(args, cycle_work_dir, cycle_results_dir)
    rt.run_pipeline(cycle_args)

    # scratch (tagged fasta, raw Metabuli output) isn't a deliverable and would
    # otherwise re-grow every cycle; rmtree only unlinks the symlink, never follows it
    shutil.rmtree(cycle_work_dir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--demux-dir', required=True, help="demux_watcher.py's persistent output dir")
    ap.add_argument('--batches-dir', required=True, help='where batchNN/ folders are created')
    ap.add_argument('--gtdb-db', required=True)
    ap.add_argument('--metabuli-env', default=None)
    ap.add_argument('--metabuli-bin', default=None, help='metabuli executable (if not using --metabuli-env)')
    ap.add_argument('--min-reads', type=int, default=1,
                     help='per-barcode read floor before a barcode is included in a cycle '
                          '(default: 1 -- classify whatever exists yet, live-monitoring '
                          'style, rather than rt_rep_seq\'s normal 10000-read quality gate)')
    ap.add_argument('--qc-min-reads', type=int, default=10000,
                     help='QC read threshold: start value of the report slider and cut-off for '
                          'the static PNG heatmap/histograms (default: 10000)')
    ap.add_argument('--run-name', default=None, help='run name shown in the report header')
    ap.add_argument('--threads', type=int, default=None,
                     help='Metabuli threads (default: all cores of the CPU budget)')
    ap.add_argument('--poll-interval', type=int, default=10,
                     help='seconds between checks for growth in --demux-dir (default: 10)')
    ap.add_argument('--latest-dir', default=None,
                     help='if given, also mirror each completed cycle into this single fixed '
                          'path (wiped and replaced every time) -- a stable "current results" '
                          'location alongside the never-overwritten batchNN/ history')
    ap.add_argument('--cpu-frac', type=float, default=0.5,
                     help='share of CPU cores this process (and Metabuli, which it spawns '
                          'each cycle) is pinned to (default: 0.5)')
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
    print(f'[analysis-watcher] resource budget: CPU {len(rt.TOOL_RESOURCES.cpus)}/{os.cpu_count()} '
          f'cores (cap {args.cpu_frac * 100:.0f}%), RAM cap '
          f'{rt.human_bytes(rt.TOOL_RESOURCES.mem_limit)}{cap_note}', flush=True)

    os.makedirs(args.batches_dir, exist_ok=True)
    results_root = os.path.dirname(os.path.abspath(args.batches_dir))
    rt.status_init(os.path.join(results_root, '.status', 'analysis.json'), state='waiting for data', pid=os.getpid(),
                   latest=os.path.abspath(args.latest_dir) if args.latest_dir else None)
    print(f'[analysis-watcher] observing {args.demux_dir} for growth -- checking every '
          f'{args.poll_interval}s -- batches -> {args.batches_dir}', flush=True)

    last_state = None
    while True:
        if os.path.exists(os.path.join(args.demux_dir, STOP_MARKER)):
            print('[analysis-watcher] stop marker found -- exiting', flush=True)
            break

        current = demux_state(args.demux_dir)
        counts_ready = os.path.exists(os.path.join(args.demux_dir, 'demux_read_counts.tsv'))
        if not current or current == last_state or not counts_ready:
            # demux_watcher writes/grows repBC*.fastq incrementally WHILE still processing
            # a raw file, but only (re)writes demux_read_counts.tsv once that file is fully
            # done -- without this check we'd race ahead of a file still being demuxed.
            time.sleep(args.poll_interval)
            continue

        batch_dir = next_batch_dir(args.batches_dir)
        batch_name = os.path.basename(batch_dir)
        print(f'\n[analysis-watcher] demux data grew -- starting {batch_name} '
              f'({len(current)} barcode file(s) present)', flush=True)
        state_at_start = current  # snapshot taken BEFORE the slow cycle runs
        rt.status_update(force=True, state='running', batch=batch_name, batch_start=time.time(),
                         batch_barcodes=len(current), step=None, progress=None, running=None, error=None)

        try:
            run_one_cycle(args, batch_dir)
            last_state = state_at_start
            print(f'[analysis-watcher] {batch_name} complete -> {batch_dir}/results', flush=True)
            rt.status_update(force=True, state='up to date', last_batch=batch_name, last_done=time.time(),
                             last_barcodes=len(state_at_start), step=None, progress=None, running=None)
            if args.latest_dir:
                # build the new copy next to latest/, then swap it in with two renames, so a
                # browser auto-reloading latest/report.html never finds the folder missing/half-copied
                tmp_dir, old_dir = args.latest_dir + '.new', args.latest_dir + '.old'
                shutil.rmtree(tmp_dir, ignore_errors=True)
                shutil.rmtree(old_dir, ignore_errors=True)
                shutil.copytree(os.path.join(batch_dir, 'results'), tmp_dir)
                if os.path.exists(args.latest_dir):
                    os.rename(args.latest_dir, old_dir)
                os.rename(tmp_dir, args.latest_dir)
                shutil.rmtree(old_dir, ignore_errors=True)
                print(f'[analysis-watcher] mirrored -> {args.latest_dir}', flush=True)
        except Exception as e:
            print(f'[analysis-watcher] {batch_name} failed ({e}); will retry next cycle',
                  file=sys.stderr, flush=True)
            rt.status_update(force=True, state='error -- retrying', error=str(e)[:300], step=None, progress=None)
            shutil.rmtree(batch_dir, ignore_errors=True)  # don't leave a half-written batch
            time.sleep(args.poll_interval)


if __name__ == '__main__':
    main()
