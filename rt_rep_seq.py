#!/usr/bin/env python3
"""RT-Rep-Seq — a compact, one-script, real-time analysis pipeline for the
On-rep-seq method (Krych et al. 2019, Communications Biology), originally
rewritten from the lab's `pyons3` reference implementation.

Steps (run in order, each writing into --work-dir / --results-dir):

  1. demux    Spacer/barcode demultiplexing. Asks how many barcodes were used in the
              run (or --num-barcodes), then stops early as soon as ceil(stop-fraction
              * num-barcodes) barcodes have reached --min-reads reads (defaults: 90%,
              10000) -- instead of always processing every input file.
  2. classify Pools every passed barcode's reads, classifies them with Metabuli against
              a GTDB-style reference DB, and resolves each classified read down to its
              lowest assigned taxonomic rank ("species" in GTDB). Writes
              results/classification/species_by_barcode.tsv: species (rows) x
              barcode (columns), read counts.
  3. lcps     Per-barcode length-count profile (LCP): how many reads fall at each read
              length. Writes one *_LCP.txt per barcode plus a single PNG with one
              histogram panel per barcode, all sharing the same x (length) / y (count)
              axes -> results/LCPs/lcp_histograms.png.
  4. heatmap  Barcode x barcode LCP-similarity heatmap (Pearson correlation of each
              barcode's length distribution, hierarchically clustered) ->
              results/LCPs/lcp_similarity_heatmap.png.
  5. report   Single self-contained results/report.html (images embedded, no external
              files needed) with the reads-per-barcode table, the species-classification
              table, the heatmap, and the histogram page(s), under the ON-rep-seq /
              University of Copenhagen banner.

On launch, prints an ASCII banner (--no-banner to skip) and the resource budget it's
running under. By default the whole tool -- including Metabuli, which it spawns -- is
pinned to half the machine's CPU cores and capped to half its RAM (a self-re-exec under
a systemd-run --user --scope cgroup; skipped for RAM, with a one-line notice, if
systemd-run isn't usable on this machine); override with --cpu-frac/--mem-frac. Each
step prints a "[i/5] <step>" banner with a live CPU/RAM usage line when it starts;
steps with a per-item loop (files, barcodes) show a progress bar, and the one step with
unpredictable duration (the Metabuli subprocess in `classify`) prints an elapsed-time/
resource heartbeat every --heartbeat-interval seconds instead of going silent.

Reuses the lab's torchlex spacer/barcode-matching engine as a plain library import
(Demultiplexer/FastqFile) -- NOT its `torchlex` command-line tool, which has a known
merge/triage bug for custom (non-default) barcode sets. This script owns its own
demux loop instead, so that bug never applies here.

See README.txt for the one-shot conda environment, full method background, and how
to read report.html.

Usage:
  python rt_rep_seq.py --input-dir /path/to/fastq_pass --gtdb-db /path/to/gtdb_db \\
      --metabuli-env metabuli [--num-barcodes 96] [--threads 8]
"""
import argparse
import base64
import json
import glob
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict

import numpy as np

try:
    import psutil
except ImportError:
    psutil = None

HERE = os.path.dirname(os.path.realpath(__file__))

VERSION = '1.1.0'
REPO_URL = 'https://github.com/krychu84/RT-Rep-Seq'
RANK_PREFIX = re.compile(r'^[a-z]__')

N_STEPS = 5  # demux, classify, lcps, heatmap, report -- for the "[i/N]" step banners

TOOL_RESOURCES = None  # set by Resources.install() at the top of main()


# ---------------------------------------------------------------------------
# startup banner
# ---------------------------------------------------------------------------

def _tty():
    return sys.stdout.isatty() and os.environ.get('TERM', '') not in ('', 'dumb')


def colour(text, code):
    """ANSI colour on a terminal, plain text when logged/redirected."""
    return f'\033[{code}m{text}\033[0m' if _tty() else text


def banner():
    """RT-Rep-Seq logo + one-line credit, printed once at launch."""
    path = os.path.join(HERE, 'assets', 'banner.txt')
    try:
        with open(path) as fh:
            logo = fh.read().rstrip('\n')
    except OSError:
        logo = 'RT-Rep-Seq'
    print()
    for line in logo.split('\n'):
        print('  ' + colour(line, '1;31'))
    spark = '\u2581\u2582\u2586\u2581\u2583\u2588\u2582\u2581\u2585\u2581\u2582\u2587\u2583\u2581\u2581\u2584\u2581\u2582\u2581\u2585\u2582\u2581\u2581\u2583\u2581'
    print('  ' + colour(spark * 2 + '\u2581' * 13, '31'))
    print(f'  {colour("Real-time ON-rep-seq analysis", "1")} \u00b7 v{VERSION} \u00b7 {REPO_URL}')
    print('  Method: Krych et al. 2019, Communications Biology 2:369, doi:10.1038/s42003-019-0617-x\n')


# ---------------------------------------------------------------------------
# resource management: CPU/RAM budget + live usage reporting
# ---------------------------------------------------------------------------

def human_bytes(n):
    """12345678 -> '11.8MB'. None -> '?' (e.g. psutil unavailable)."""
    if n is None:
        return '?'
    n = float(n)
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if abs(n) < 1024 or unit == 'TB':
            return f'{n:.1f}{unit}'
        n /= 1024


class Resources:
    """Keeps the tool within a CPU/RAM budget (default: half the machine) and
    reports live usage.

    CPU is capped by pinning this process to the upper `cpu_frac` share of the
    machine's cores (`os.sched_setaffinity`) plus a nice(10) -- since affinity is
    inherited, every subprocess this tool spawns (Metabuli included) is capped the
    same way automatically, with no per-command wrapping needed.

    RAM is capped by re-executing the ENTIRE process tree once, at startup, under a
    `systemd-run --user --scope` cgroup with a MemoryMax -- this covers the current
    process's own in-process work (the demux loop, LCP/plot building, ...), not just
    external subprocesses. If systemd-run isn't available or usable (e.g. no user
    cgroup delegation on this machine), RAM capping is skipped with a one-line notice
    and the tool keeps running uncapped for memory -- CPU capping still applies.
    """
    REEXEC_SENTINEL = 'RT_REP_SEQ_MEM_CAPPED'

    def __init__(self, cpu_frac=0.5, mem_frac=0.5):
        self.cpu_frac = cpu_frac
        self.mem_frac = mem_frac
        n_cpu = os.cpu_count() or 1
        self.n_threads = max(1, round(n_cpu * cpu_frac))
        self.cpus = list(range(n_cpu - self.n_threads, n_cpu))
        self.mem_total = psutil.virtual_memory().total if psutil else None
        self.mem_limit = int(self.mem_total * mem_frac) if self.mem_total else None
        self.mem_capped = os.environ.get(self.REEXEC_SENTINEL) == '1'
        self._proc = psutil.Process() if psutil else None
        if psutil:
            psutil.cpu_percent(None, percpu=True)  # prime the counter (first call is always 0)

    def reexec_with_mem_cap(self, argv):
        """Re-exec this same invocation under a memory-capped systemd scope, once.
        Returns normally (a no-op) if already capped, there's no RAM budget to
        enforce, or systemd-run isn't usable -- the caller just continues uncapped.
        Never returns if the re-exec succeeds (the process image is replaced)."""
        if self.mem_capped or not self.mem_limit or not shutil.which('systemd-run'):
            if not self.mem_capped and self.mem_limit:
                print('[RT-Rep-Seq] note: systemd-run not found -- RAM cap skipped '
                      '(CPU cap still applies)', flush=True)
            return
        probe = subprocess.run(['systemd-run', '--user', '--scope', '--quiet',
                                 '-p', 'MemoryMax=1G', 'true'],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if probe.returncode != 0:
            print('[RT-Rep-Seq] note: systemd-run user scope not available on this '
                  'machine -- RAM cap skipped (CPU cap still applies)', flush=True)
            return
        env = dict(os.environ)
        env[self.REEXEC_SENTINEL] = '1'
        cmd = ['systemd-run', '--user', '--scope', '--quiet',
               '-p', f'MemoryMax={self.mem_limit}', '-p', 'MemorySwapMax=0',
               sys.executable] + list(argv)
        os.execvpe(cmd[0], cmd, env)

    def apply_cpu_cap(self):
        try:
            os.sched_setaffinity(0, self.cpus)
        except (AttributeError, OSError):
            pass
        try:
            os.nice(10)
        except OSError:
            pass

    def _children(self):
        try:
            return self._proc.children(recursive=True)
        except psutil.Error:
            return []

    def usage(self):
        if not psutil:
            return None
        n_cpu = os.cpu_count() or 1
        per = psutil.cpu_percent(None, percpu=True)
        pinned = [per[i] for i in self.cpus if i < len(per)]
        rss = 0
        for p in [self._proc] + self._children():
            try:
                rss += p.memory_info().rss
            except psutil.Error:
                pass
        vm = psutil.virtual_memory()
        return {
            'tool_cpu': (sum(pinned) / n_cpu) if per else 0.0,
            'tool_rss': rss,
            'sys_cpu': (sum(per) / len(per)) if per else 0.0,
            'sys_mem': vm.percent,
        }

    def line(self):
        """One-line usage summary for step banners/heartbeats."""
        u = self.usage()
        if u is None:
            return _basic_resource_snapshot()
        cap = human_bytes(self.mem_limit) + ('' if self.mem_capped else ', not enforced')
        return (f"CPU tool {u['tool_cpu']:4.1f}% (cap {self.cpu_frac * 100:.0f}%) · "
                f"sys {u['sys_cpu']:4.1f}%  |  RAM tool {human_bytes(u['tool_rss'])} "
                f"/ cap {cap} · sys {u['sys_mem']:.0f}%")


# ---------------------------------------------------------------------------
# progress / resource reporting
# ---------------------------------------------------------------------------

def _meminfo_gb():
    """(used, total) system RAM in GB, from /proc/meminfo (Linux only). None elsewhere."""
    try:
        info = {}
        with open('/proc/meminfo') as f:
            for line in f:
                key, val = line.split(':', 1)
                info[key] = int(val.strip().split()[0])  # kB
        total = info['MemTotal'] / 1e6
        avail = info.get('MemAvailable', info.get('MemFree', 0)) / 1e6
        return total - avail, total
    except Exception:
        return None


def _basic_resource_snapshot():
    """One-line 'RAM x/y GB | CPU load a/b cores' fallback when psutil (and
    therefore a real Resources budget) isn't available."""
    parts = []
    mem = _meminfo_gb()
    if mem:
        used, total = mem
        parts.append(f'RAM {used:.1f}/{total:.1f} GB')
    try:
        load1 = os.getloadavg()[0]
        ncpu = os.cpu_count() or 1
        parts.append(f'CPU load {load1:.1f}/{ncpu} cores')
    except (OSError, AttributeError):
        pass
    return ' | '.join(parts) if parts else 'resource info unavailable'


def resource_snapshot():
    """One-line resource summary for step banners/heartbeats -- the real CPU/RAM-
    budget line once Resources.install() has run (main()), else a basic fallback
    (e.g. a step function imported and called directly, outside this tool's CLI)."""
    return TOOL_RESOURCES.line() if TOOL_RESOURCES else _basic_resource_snapshot()


# Live status for the launcher's terminal view: when STATUS['path'] is set (by the
# watchers), every step banner / progress bar / heartbeat also updates a small JSON file.
STATUS = {'path': None, 'data': {}, 'last': 0.0}


def status_init(path, **base):
    STATUS.update(path=path, data=dict(base), last=0.0)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    status_update(force=True)


def status_update(force=False, **kw):
    if not STATUS['path']:
        return
    STATUS['data'].update(kw)
    now = time.time()
    STATUS['data']['updated'] = now
    if not force and now - STATUS['last'] < 0.5:
        return
    STATUS['last'] = now
    tmp = STATUS['path'] + '.tmp'
    try:
        with open(tmp, 'w') as f:
            json.dump(STATUS['data'], f)
        os.replace(tmp, STATUS['path'])
    except OSError:
        pass


STEP_NAMES = {1: 'Demultiplexing', 2: 'Species classification (Metabuli)', 3: 'Read-length profiles (LCPs)',
              4: 'LCP similarity + strain groups', 5: 'HTML report'}


def step_banner(step_no, name):
    print(f'\n{colour(f"[{step_no}/{N_STEPS}] {name}", "1;36")}  ({resource_snapshot()})')
    status_update(force=True, step=step_no, n_steps=N_STEPS, step_name=name, step_start=time.time(),
                  progress=None, running=None)


def progress_bar(current, total, prefix='', width=30):
    """Single-line progress bar; overwrites in place on a TTY, prints one line per
    update otherwise (so it stays readable in a redirected/logged run)."""
    if total <= 0:
        return
    frac = min(current / total, 1.0)
    filled = int(width * frac)
    tty = sys.stdout.isatty()
    bar = ('\u2588' * filled + '\u2591' * (width - filled)) if tty else ('#' * filled + '-' * (width - filled))
    end = '\n' if (current >= total or not tty) else ''
    print(f'\r  {prefix} [{bar}] {current}/{total} ({frac * 100:5.1f}%)', end=end, flush=True)
    status_update(force=current >= total, progress=[current, total, prefix])


def run_with_heartbeat(cmd, interval=15):
    """Like subprocess.run(cmd, check=True), but prints an elapsed-time + resource
    heartbeat every `interval` seconds while the (unpredictable-duration) subprocess
    runs, instead of going silent until it finishes."""
    start = time.time()
    proc = subprocess.Popen(cmd)
    status_update(force=True, progress=None, running=os.path.basename(cmd[0]) if 'conda' not in cmd[0] else 'metabuli',
                  running_start=start)
    while proc.poll() is None:
        time.sleep(interval)
        if proc.poll() is None:
            elapsed = time.time() - start
            print(f'  ... still running ({elapsed:.0f}s elapsed, {resource_snapshot()})')
            status_update(force=True)
    status_update(force=True, running=None)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)


# ---------------------------------------------------------------------------
# small shared helpers
# ---------------------------------------------------------------------------

def bc_sort_key(bc):
    m = re.sub(r'\D', '', bc)
    return int(m) if m else 0


def load_nodes(path):
    """GTDB/NCBI-style nodes.dmp -> ({taxid: parent_taxid}, {taxid: rank})."""
    par, rank = {}, {}
    for line in open(path):
        c = [x.strip() for x in line.split('|')]
        if len(c) >= 3:
            par[c[0]] = c[1]
            rank[c[0]] = c[2]
    return par, rank


def load_names(path):
    """GTDB/NCBI-style names.dmp -> {taxid: scientific_name}."""
    names = {}
    for line in open(path):
        c = [x.strip() for x in line.split('|')]
        if len(c) >= 2:
            cls = c[3] if len(c) >= 4 else ''
            if c[0] not in names or cls == 'scientific name':
                names[c[0]] = c[1]
    return names


def to_species(taxid, par, rank):
    """Walk up the taxonomy from `taxid` to the nearest ancestor ranked 'species'."""
    seen = set()
    t = taxid
    while t and t not in seen and t != '1':
        seen.add(t)
        if rank.get(t) == 'species':
            return t
        t = par.get(t)
    return ''


def length_profile(fastq_path):
    """Counter: read length -> number of reads, for one FASTQ file."""
    lengths = Counter()
    with open(fastq_path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lengths[len(line.rstrip('\n'))] += 1
    return lengths


# ---------------------------------------------------------------------------
# step 1 -- demux (early-stopping)
# ---------------------------------------------------------------------------

def step_demux(args):
    step_banner(1, 'Demultiplexing raw reads (torchlex)')

    from torchlex.config import load_config
    from torchlex.demultiplexer import Demultiplexer
    from torchlex.fastq_reader import FastqFile

    n_barcodes = args.num_barcodes
    if n_barcodes is None:
        n_barcodes = int(input("How many barcodes were used in this run? ").strip())

    need_passed = math.ceil(args.stop_fraction * n_barcodes)
    print(f"[demux] target: {need_passed}/{n_barcodes} barcodes "
          f"({args.stop_fraction * 100:.0f}%) reaching >= {args.min_reads} reads")

    config = load_config(args.onrep_config)
    config['amplicones']['barcodes_path'] = args.onrep_barcodes
    demux = Demultiplexer(config)

    out_dir = os.path.join(args.work_dir, 'demultiplexed_fastq')
    os.makedirs(out_dir, exist_ok=True)

    fastq_paths = sorted(glob.glob(os.path.join(args.input_dir, '*.fastq')))
    fastq_paths += sorted(glob.glob(os.path.join(args.input_dir, '*.fastq.gz')))
    if not fastq_paths:
        sys.exit(f"[demux] no *.fastq(.gz) files found in {args.input_dir}")

    counts = Counter()
    processed = 0
    for i, path in enumerate(fastq_paths, 1):
        fq = FastqFile(path)
        classified = demux.demultiplex_file(fq, out_dir)  # appends {barcode}.fastq in out_dir
        for row in classified.rows:
            counts[row.barcode] += 1
        processed = i
        n_passed = sum(1 for c in counts.values() if c >= args.min_reads)
        lead = '\n' if sys.stdout.isatty() else ''  # keep the overwriting bar below on a TTY
        print(f"{lead}[demux] {os.path.basename(path)} "
              f"({len(fq)} reads, {len(classified.rows)} classified) | "
              f"barcodes >= {args.min_reads} reads: {n_passed}/{n_barcodes} "
              f"({resource_snapshot()})")
        progress_bar(i, len(fastq_paths), prefix='input files')
        if n_passed >= need_passed:
            print(f"[demux] stopping early after {i}/{len(fastq_paths)} files: "
                  f"{n_passed}/{n_barcodes} barcodes reached {args.min_reads} reads "
                  f"(target was {need_passed})")
            break
    else:
        n_passed = sum(1 for c in counts.values() if c >= args.min_reads)
        print(f"[demux] processed all {len(fastq_paths)} input files "
              f"({n_passed}/{n_barcodes} barcodes reached the threshold)")

    passed = {bc: n for bc, n in counts.items() if n >= args.min_reads}

    counts_path = os.path.join(out_dir, 'demux_read_counts.tsv')
    with open(counts_path, 'w') as f:
        f.write('barcode\treads\tpassed\n')
        for bc, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            f.write(f"{bc}\t{n}\t{'yes' if n >= args.min_reads else 'no'}\n")

    print(f"[demux] {len(passed)} barcodes passed (>= {args.min_reads} reads) "
          f"of {len(counts)} seen, {processed} files read -> {counts_path}")
    args._files_read = processed
    return passed, dict(counts)


def load_passed_from_counts(work_dir, min_reads=None):
    """Barcodes from a prior/external demux pass's demux_read_counts.tsv.

    min_reads=None (default) trusts the file's own 'passed' column as-is -- whatever
    threshold whoever wrote it used. Pass an explicit min_reads to instead recompute
    'passed' from the raw per-barcode read counts against *this* threshold: needed
    when the reader's quality bar differs from the writer's (e.g. analysis_watcher.py
    wants a low live-monitoring floor like 1, but demux_watcher.py bakes in its own,
    usually much higher, production threshold like 10000 -- trusting the file's
    column in that case would silently keep every barcode filtered out no matter what
    --min-reads this run was given).
    """
    counts_path = os.path.join(work_dir, 'demultiplexed_fastq', 'demux_read_counts.tsv')
    passed = {}
    with open(counts_path) as f:
        next(f)
        for line in f:
            bc, n, p = line.rstrip('\n').split('\t')
            n = int(n)
            ok = (n >= min_reads) if min_reads is not None else (p == 'yes')
            if ok:
                passed[bc] = n
    return passed


# ---------------------------------------------------------------------------
# step 2 -- classify (Metabuli)
# ---------------------------------------------------------------------------

def step_classify(args, passed):
    step_banner(2, 'Classifying reads (Metabuli)')

    demux_dir = os.path.join(args.work_dir, 'demultiplexed_fastq')
    class_dir = os.path.join(args.work_dir, 'classify')
    results_dir = os.path.join(args.results_dir, 'classification')
    os.makedirs(class_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    barcodes = sorted(passed, key=bc_sort_key)
    tagged_fasta = os.path.join(class_dir, 'all_tagged.fasta')
    n_reads = n_short = 0
    with open(tagged_fasta, 'w') as out:
        for i, bc in enumerate(barcodes, 1):
            fq_path = os.path.join(demux_dir, f'{bc}.fastq')
            if not os.path.exists(fq_path):
                continue
            with open(fq_path) as fh:
                while True:
                    h = fh.readline()
                    if not h:
                        break
                    s = fh.readline()
                    fh.readline()
                    fh.readline()
                    if len(s.rstrip('\n')) < args.min_read_length:
                        n_short += 1
                        continue
                    rid = h[1:].split()[0].strip()
                    out.write(f'>{bc}__{rid}\n{s}')
                    n_reads += 1
            progress_bar(i, len(barcodes), prefix='tagging barcodes')
    print(f'[classify] tagged {n_reads} reads >= {args.min_read_length} bp from {len(barcodes)} barcodes '
          f'({n_short} shorter reads not classified) -> {tagged_fasta}')

    jobid = 'rt_rep_seq'
    if args.metabuli_env:
        cmd = ['conda', 'run', '-n', args.metabuli_env, 'metabuli', 'classify']
    else:
        cmd = [args.metabuli_bin, 'classify']
    threads = args.threads or (len(TOOL_RESOURCES.cpus) if TOOL_RESOURCES else 8)
    cmd += ['--seq-mode', '3', '--threads', str(threads)]
    max_ram = args.metabuli_max_ram
    if max_ram is None and TOOL_RESOURCES and TOOL_RESOURCES.mem_limit:
        # leave ~25% of the RT-Rep-Seq RAM budget for Python + the demux watcher
        max_ram = max(4, int(TOOL_RESOURCES.mem_limit / 2 ** 30 * 0.75))
    if max_ram:
        cmd += ['--max-ram', str(max_ram)]
    cmd += [tagged_fasta, args.gtdb_db, class_dir, jobid]
    print(f"[classify] running: {' '.join(cmd)}  ({resource_snapshot()})")
    run_with_heartbeat(cmd, interval=args.heartbeat_interval)

    classifications = os.path.join(class_dir, f'{jobid}_classifications.tsv')
    par, rank = load_nodes(args.taxonomy_nodes)
    names = load_names(args.taxonomy_names)

    per_bc = defaultdict(Counter)
    for line in open(classifications):
        if line.startswith('#'):
            continue
        c = line.rstrip('\n').split('\t')
        if len(c) < 3 or '__' not in c[1]:
            continue
        bc = c[1].split('__', 1)[0]
        if c[0] != '1' or c[2] == '0':
            continue
        sp = to_species(c[2], par, rank)
        if not sp:
            continue
        per_bc[bc][RANK_PREFIX.sub('', names.get(sp, 'taxid:' + sp))] += 1

    barcodes = sorted(per_bc, key=bc_sort_key)
    all_species = sorted({sp for c in per_bc.values() for sp in c})

    # Per-barcode noise floor: this is pure-isolate data (one organism per barcode), so any
    # species reaching less than min-species-pct of that barcode's classified reads is treated
    # as index-hopping/misclassification noise, not a real call, and zeroed out.
    thresholded = defaultdict(dict)
    barcode_totals = {bc: sum(per_bc[bc].values()) for bc in barcodes}
    for bc in barcodes:
        tot = barcode_totals[bc] or 1
        for sp, n in per_bc[bc].items():
            thresholded[bc][sp] = n if (100 * n / tot) >= args.min_species_pct else 0

    species = [sp for sp in all_species if any(thresholded[bc].get(sp, 0) > 0 for bc in barcodes)]

    table_path = os.path.join(results_dir, 'species_by_barcode.tsv')
    with open(table_path, 'w') as f:
        f.write('species\t' + '\t'.join(barcodes) + '\n')
        for sp in species:
            f.write(sp + '\t' + '\t'.join(str(thresholded[bc].get(sp, 0)) for bc in barcodes) + '\n')

    print(f'[classify] {len(species)} species x {len(barcodes)} barcodes '
          f'(species below {args.min_species_pct:.0f}% of a barcode\'s reads zeroed as noise) -> {table_path}')

    # same table WITHOUT the noise filter -- every species-level call -- for the report's
    # per-barcode composition bars
    all_path = os.path.join(results_dir, 'species_by_barcode_all.tsv')
    with open(all_path, 'w') as f:
        f.write('species\t' + '\t'.join(barcodes) + '\n')
        for sp in all_species:
            f.write(sp + '\t' + '\t'.join(str(per_bc[bc].get(sp, 0)) for bc in barcodes) + '\n')
    print(f'[classify] unfiltered table ({len(all_species)} species) -> {all_path}')

    # dominant species per barcode, from the RAW (pre-threshold) counts -- used for the
    # console preview and as the heatmap's left-side species legend
    dominant = {bc: per_bc[bc].most_common(1)[0][0] for bc in barcodes if per_bc[bc]}
    print('[classify] dominant (lowest-level) taxon per barcode:')
    for bc in barcodes:
        tot = barcode_totals[bc]
        top_sp, top_n = per_bc[bc].most_common(1)[0]
        print(f'  {bc:10} {top_sp:35} {top_n:>8} / {tot:>8} classified reads')
    return table_path, dominant


# ---------------------------------------------------------------------------
# step 3 -- LCPs (length-count profiles) + unified-axis histogram grid
# ---------------------------------------------------------------------------

def bin_profile(counter, binsize):
    """Sum a {length: count} Counter into binsize-wide bins, keyed by each bin's start."""
    binned = Counter()
    for length, n in counter.items():
        binned[(length // binsize) * binsize] += n
    return binned


def auto_xmax(binned_profiles, binsize, min_count=20):
    """Longest bin (bp) that reaches >= min_count reads in at least one sample.

    A handful of long outlier reads (each seen only once or twice) shouldn't be able
    to stretch every panel's x-axis out and squash the real signal into a sliver on
    the left, so the axis only has to reach as far as the longest *well-represented*
    peak, in any one sample. Falls back to the plain data max if no bin anywhere
    reaches min_count (e.g. a very shallow sample), so the axis is never degenerate.
    """
    qualifying = [length for p in binned_profiles.values()
                  for length, n in p.items() if n >= min_count]
    all_bins = [b for p in binned_profiles.values() for b in p]
    longest = max(qualifying) if qualifying else max(all_bins)
    return longest + binsize


def auto_ymax(binned_profiles):
    """Highest single-bin read count across all samples."""
    return max(n for p in binned_profiles.values() for n in p.values())


def similarity_matrix_and_order(profiles, min_len=500):
    """Barcode x barcode LCP similarity from band (peak) presence -- weighted Dice with a
    position tolerance (lcp_peaks.py) -- on reads >= min_len, ordered by UPGMA. Returns
    (barcodes_ord, sim_ord); sim_ord is None if fewer than 2 barcodes have data."""
    import lcp_peaks
    barcodes = sorted(profiles, key=bc_sort_key)
    if len(barcodes) < 2:
        return barcodes, None
    peaks = {bc: lcp_peaks.call_peaks({L: c for L, c in profiles[bc].items() if L >= min_len},
                                      min_len=min_len) for bc in barcodes}
    sim = np.clip(lcp_peaks.similarity_matrix(peaks, barcodes), 0, 1)
    order = lcp_peaks.upgma_order(sim, list(range(len(barcodes))))
    return [barcodes[i] for i in order], sim[np.ix_(order, order)]


def group_by_species_columns(barcodes, dominant):
    """Group barcodes into species -> [barcodes] columns, in order of each species'
    first appearance in `barcodes`; row order within a column also follows `barcodes`.
    Barcodes with no dominant-species entry get their own 'unclassified' column."""
    columns = {}
    col_order = []
    for bc in barcodes:
        species = dominant.get(bc) or 'unclassified'
        if species not in columns:
            columns[species] = []
            col_order.append(species)
        columns[species].append(bc)
    return [(sp, columns[sp]) for sp in col_order]


def _render_grid_page(binned_profiles, grid, nrows, ncols, binsize, xmin, xmax, ymax,
                       dominant, panel_width, panel_height, species_font_size, png_path):
    """Draw and save one page (a fixed grid of barcode/species columns) of the LCP
    histogram figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_width * ncols, panel_height * nrows),
                              sharex=True, sharey=True, squeeze=False)
    for row_idx in range(nrows):
        for col_idx in range(ncols):
            ax = axes[row_idx, col_idx]
            bc = grid[row_idx][col_idx]
            if bc is None:
                ax.axis('off')
                continue
            p = binned_profiles[bc]
            bins = sorted(p)
            ax.bar(bins, [p[b] for b in bins], width=binsize, align='edge', color='#3b6fa0')
            ax.set_title(bc, fontsize=8)
            ax.tick_params(labelsize=6)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(0, ymax)
            species = dominant.get(bc)
            if species:
                ax.set_xlabel(species, fontsize=species_font_size, fontstyle='italic')
    fig.supxlabel(f'read length (bp, {binsize}bp bins)')
    fig.supylabel('read count')
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return png_path


def plot_lcp_histogram_grid(profiles, binsize, png_path, xmax=None, ymax=None, xmin=0,
                             panel_width=4.5, panel_height=1.8, order=None, dominant=None,
                             min_count_for_xmax=20, cols_per_page=None, species_font_size=11):
    """Render the LCP histogram grid for {barcode: {length: count}} profiles, one column
    per bacterial species: every barcode dominated by the same species is stacked in the
    same column (top to bottom), and a barcode with a new dominant species starts a new
    column. Falls back to a plain 8-wide grid if `dominant` isn't given (no species info
    to group by).

    xmax/ymax=None auto-derive robust, non-outlier-driven bounds (see auto_xmax/auto_ymax)
    instead of either the raw data range or a hardcoded number, and are shared across
    every page so pages stay visually comparable.
    panel_width/panel_height (inches, per panel) control the figure's aspect ratio --
    raise panel_width to stretch the x-axis so individual bins are easier to tell apart.
    order fixes the barcode sequence used to build columns (species order = first
    appearance in `order`; row order within a column also follows `order`) -- pass the
    heatmap's hierarchical-clustering order so barcodes of the same/related species land
    next to each other. Barcodes missing from profiles are skipped, and any barcode in
    profiles but not in `order` is appended at the end. Defaults to numeric/alphabetic
    barcode order.
    dominant maps barcode -> its top-hit species name; also printed under each panel.
    cols_per_page splits the species columns into separate image files (page breaks) of
    at most that many columns each -- e.g. 4 -> lcp_histograms_page1.png, _page2.png, ...
    -- instead of one wide sparse image. None (default) keeps everything on one page.

    Returns the list of PNG paths written (length 1 unless paginated).
    """
    binned_profiles = {bc: bin_profile(p, binsize) for bc, p in profiles.items()}

    if order is not None:
        barcodes = [bc for bc in order if bc in binned_profiles]
        barcodes += [bc for bc in binned_profiles if bc not in barcodes]
    else:
        barcodes = sorted(binned_profiles, key=bc_sort_key)

    if xmax is None:
        xmax = auto_xmax(binned_profiles, binsize, min_count=min_count_for_xmax)
    if ymax is None:
        ymax = auto_ymax(binned_profiles)
    if xmin is None:
        xmin = min(b for p in binned_profiles.values() for b in p)

    dominant = dominant or {}

    if dominant:
        species_cols = group_by_species_columns(barcodes, dominant)
    else:
        # no species info: keep the old plain 8-wide, row-filled grid, as a single
        # "column group" of up to 8 pseudo-columns so the paging logic below still works
        species_cols = [(None, [bc]) for bc in barcodes]
        cols_per_page = cols_per_page or 8

    if cols_per_page is None:
        pages = [species_cols]
    else:
        pages = [species_cols[i:i + cols_per_page] for i in range(0, len(species_cols), cols_per_page)]

    root, ext = os.path.splitext(png_path)
    paths = []
    for page_num, page_cols in enumerate(pages, start=1):
        ncols = len(page_cols)
        nrows = max(len(bcs) for _, bcs in page_cols)
        grid = [[None] * ncols for _ in range(nrows)]
        for col_idx, (_, bcs) in enumerate(page_cols):
            for row_idx, bc in enumerate(bcs):
                grid[row_idx][col_idx] = bc
        page_path = png_path if len(pages) == 1 else f'{root}_page{page_num}{ext}'
        _render_grid_page(binned_profiles, grid, nrows, ncols, binsize, xmin, xmax, ymax,
                           dominant, panel_width, panel_height, species_font_size, page_path)
        paths.append(page_path)
    return paths


def step_lcps(args, passed, dominant=None, plot_set=None):
    step_banner(3, 'Building per-barcode length-count profiles (LCPs)')

    demux_dir = os.path.join(args.work_dir, 'demultiplexed_fastq')
    out_dir = os.path.join(args.results_dir, 'LCPs')
    raw_dir = os.path.join(out_dir, 'raw_LCPs')
    os.makedirs(raw_dir, exist_ok=True)

    barcodes = sorted(passed, key=bc_sort_key)
    profiles = {}
    for i, bc in enumerate(barcodes, 1):
        fq_path = os.path.join(demux_dir, f'{bc}.fastq')
        if not os.path.exists(fq_path):
            progress_bar(i, len(barcodes), prefix='barcodes')
            continue
        lengths = length_profile(fq_path)
        profiles[bc] = lengths
        with open(os.path.join(raw_dir, f'{bc}_LCP.txt'), 'w') as f:
            for length in sorted(lengths):
                f.write(f'{lengths[length]} {length}\n')
        progress_bar(i, len(barcodes), prefix='barcodes')

    print(f'[lcps] {len(profiles)} barcode LCPs -> {raw_dir}/*_LCP.txt')
    all_profiles = profiles
    if plot_set is not None:  # static PNG grid = QC-passing barcodes only; the HTML report
        profiles = {bc: p for bc, p in profiles.items() if bc in plot_set}  # covers the rest
    barcodes = [bc for bc in barcodes if bc in profiles]
    if not barcodes:
        print(f'[lcps] no barcode >= {args.min_reads} reads -- static histogram PNG skipped '
              f'(the interactive report.html still shows every barcode)')
        return all_profiles

    binsize = args.lcp_bin_bp
    profiles = {bc: {L: c for L, c in p.items() if L >= args.min_read_length} for bc, p in profiles.items()}
    profiles = {bc: p for bc, p in profiles.items() if p}
    order, _ = similarity_matrix_and_order(profiles, args.min_read_length)
    cols_per_page = args.lcp_cols_per_page or None  # 0 -> None (single unpaginated image)
    png_path = os.path.join(out_dir, 'lcp_histograms.png')
    png_paths = plot_lcp_histogram_grid(
        profiles, binsize, png_path, xmax=args.lcp_xmax, ymax=args.lcp_ymax,
        order=order, dominant=dominant, min_count_for_xmax=args.lcp_min_count_for_xmax,
        cols_per_page=cols_per_page)

    print(f'[lcps] combined histogram grid ({len(barcodes)} barcodes >= {args.min_reads} reads,  ({binsize}bp bins, x max auto >= '
          f'{args.lcp_min_count_for_xmax} reads/bin, y max auto, ordered by LCP similarity, '
          f'{cols_per_page or "all"} species columns/page) -> {", ".join(png_paths)}')
    return all_profiles


# ---------------------------------------------------------------------------
# step 4 -- LCP similarity heatmap
# ---------------------------------------------------------------------------

def load_dominant_from_table(results_dir):
    """Fallback: recover each barcode's top species from species_by_barcode.tsv, for a run
    where --skip classify reused a prior classification but still wants the heatmap legend."""
    table_path = os.path.join(results_dir, 'classification', 'species_by_barcode.tsv')
    if not os.path.exists(table_path):
        return {}
    with open(table_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        barcodes = header[1:]
        best = {bc: (0, 'unclassified') for bc in barcodes}
        for line in f:
            c = line.rstrip('\n').split('\t')
            sp, counts = c[0], c[1:]
            for bc, v in zip(barcodes, counts):
                n = int(v)
                if n > best[bc][0]:
                    best[bc] = (n, sp)
    return {bc: sp for bc, (n, sp) in best.items()}


def step_heatmap(args, profiles, dominant=None):
    step_banner(4, 'Clustering barcodes by LCP similarity (heatmap)')

    out_dir = os.path.join(args.results_dir, 'LCPs')
    os.makedirs(out_dir, exist_ok=True)

    barcodes_ord, sim_ord = similarity_matrix_and_order(profiles, args.min_read_length)
    if sim_ord is None:
        print('[heatmap] fewer than 2 barcodes with data -- skipping')
        return barcodes_ord

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    side = max(0.25 * len(barcodes_ord) + 3, 9)  # floor so title/colorbar never get clipped
    fig, ax = plt.subplots(figsize=(side, side))
    im = ax.imshow(sim_ord, cmap='viridis', vmin=0.5, vmax=1)
    dominant = dominant or {}
    y_labels = [f"{bc} — {dominant[bc]}" if dominant.get(bc) else bc for bc in barcodes_ord]

    ax.set_xticks(range(len(barcodes_ord)))
    ax.set_xticklabels(barcodes_ord, rotation=90, fontsize=6)
    ax.set_yticks(range(len(barcodes_ord)))
    ax.set_yticklabels(y_labels, fontsize=6)
    ax.set_title(f'LCP similarity between barcodes (band presence, weighted Dice; UPGMA; reads >= {args.min_read_length} bp)')
    fig.colorbar(im, ax=ax, shrink=0.8, label='similarity')
    fig.tight_layout()
    png_path = os.path.join(out_dir, 'lcp_similarity_heatmap.png')
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    tsv_path = os.path.join(out_dir, 'lcp_similarity_matrix.tsv')
    with open(tsv_path, 'w') as f:
        f.write('barcode\t' + '\t'.join(barcodes_ord) + '\n')
        for i, bc in enumerate(barcodes_ord):
            f.write(bc + '\t' + '\t'.join(f'{v:.3f}' for v in sim_ord[i]) + '\n')

    print(f'[heatmap] {len(barcodes_ord)} x {len(barcodes_ord)} similarity matrix -> {tsv_path}')
    print(f'[heatmap] heatmap -> {png_path}')
    return barcodes_ord


def load_profiles(args, passed):
    demux_dir = os.path.join(args.work_dir, 'demultiplexed_fastq')
    profiles = {}
    for bc in passed:
        fq_path = os.path.join(demux_dir, f'{bc}.fastq')
        if os.path.exists(fq_path):
            profiles[bc] = length_profile(fq_path)
    return profiles


# ---------------------------------------------------------------------------
# step 5 -- HTML report
# ---------------------------------------------------------------------------

def _b64_image(path):
    """data: URI for an image file (embedded so the report is one portable file), or
    None if the file doesn't exist."""
    if not path or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lstrip('.').lower()
    mime = {'png': 'image/png', 'jpg': 'image/jpeg', 'jpeg': 'image/jpeg'}.get(ext, 'application/octet-stream')
    with open(path, 'rb') as f:
        data = base64.b64encode(f.read()).decode('ascii')
    return f'data:{mime};base64,{data}'


def _html_escape(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def _reads_per_barcode_table(work_dir, min_reads):
    counts_path = os.path.join(work_dir, 'demultiplexed_fastq', 'demux_read_counts.tsv')
    if not os.path.exists(counts_path):
        return None
    rows = []
    with open(counts_path) as f:
        next(f)  # header
        for line in f:
            bc, n, p = line.rstrip('\n').split('\t')
            rows.append((bc, int(n), p == 'yes'))
    rows.sort(key=lambda r: -r[1])
    body = ''.join(
        f'<tr class="{"pass" if ok else "fail"}"><td>{_html_escape(bc)}</td>'
        f'<td>{n:,}</td><td>{"passed" if ok else f"below {min_reads:,}"}</td></tr>'
        for bc, n, ok in rows)
    return ('<table><thead><tr><th>Barcode</th><th>Reads</th><th>Status</th></tr></thead>'
            f'<tbody>{body}</tbody></table>')


def _species_summary_table(results_dir, dominant):
    """One row per barcode: its dominant species, that species' read count, and the
    barcode's total classified reads -- read straight from species_by_barcode.tsv."""
    table_path = os.path.join(results_dir, 'classification', 'species_by_barcode.tsv')
    if not os.path.exists(table_path) or not dominant:
        return None
    with open(table_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        barcodes = header[1:]
        totals = dict.fromkeys(barcodes, 0)
        top = {bc: (0, dominant.get(bc, '')) for bc in barcodes}
        for line in f:
            c = line.rstrip('\n').split('\t')
            sp, counts = c[0], c[1:]
            for bc, v in zip(barcodes, counts):
                n = int(v)
                totals[bc] += n
                if n > top[bc][0]:
                    top[bc] = (n, sp)
    body = ''.join(
        f'<tr><td>{_html_escape(bc)}</td><td><i>{_html_escape(top[bc][1] or "unclassified")}</i></td>'
        f'<td>{top[bc][0]:,}</td><td>{totals[bc]:,}</td></tr>'
        for bc in sorted(barcodes, key=bc_sort_key))
    return ('<table><thead><tr><th>Barcode</th><th>Dominant species</th>'
            '<th>Reads (dominant species)</th><th>Reads (total classified)</th></tr></thead>'
            f'<tbody>{body}</tbody></table>')


def step_report(args, dominant=None):
    """Single self-contained interactive report.html (see report_html.py): read-count
    threshold slider, 2 x 96-well plate maps, species table, live re-clustered heatmap,
    LCP histograms, TSV/PNG downloads, A4 print layout."""
    step_banner(5, 'Writing HTML report')
    import report_html

    results_dir = args.results_dir
    counts_src = os.path.join(args.work_dir, 'demultiplexed_fastq', 'demux_read_counts.tsv')
    counts_dst = os.path.join(results_dir, 'demux_read_counts.tsv')
    if os.path.exists(counts_src):
        shutil.copyfile(counts_src, counts_dst)  # results/ stays self-contained

    meta = {
        'version': VERSION,
        'run': args.run_name or os.path.basename(os.path.normpath(
            os.path.dirname(os.path.abspath(args.input_dir))
            if os.path.basename(os.path.normpath(args.input_dir)).startswith(('barcode', 'fastq'))
            else args.input_dir)),
        'input_dir': os.path.abspath(args.input_dir),
        'num_barcodes': args.num_barcodes,
        'files_read': getattr(args, '_files_read', None),
        'repo_url': REPO_URL,
    }
    n_ref = 0
    try:
        n_ref = len({l[1:].split('_')[0] for l in open(args.onrep_barcodes) if l.startswith('>')})
    except Exception:
        pass
    meta['methods'] = report_html.methods_meta(
        onrep_config=args.onrep_config, gtdb_db=args.gtdb_db, metabuli_env=args.metabuli_env,
        metabuli_bin=args.metabuli_bin, version=VERSION,
        tool_config=dict({'repo_url': REPO_URL}, **report_html.load_tool_config(HERE)),
        params={'min_reads': args.min_reads, 'stop_fraction': args.stop_fraction,
                'min_read_length': args.min_read_length,
                'min_species_pct': args.min_species_pct, 'n_reference_barcodes': n_ref or None})
    payload = report_html.build_payload(
        results_dir, counts_path=counts_dst if os.path.exists(counts_dst) else None,
        default_threshold=args.report_min_reads if args.report_min_reads is not None else args.min_reads,
        binsize=args.lcp_bin_bp, meta=meta, min_len=args.min_read_length,
        logo_ucph=args.logo_ucph, logo_onrepseq=args.logo_onrepseq)
    report_path = report_html.write_report(os.path.join(results_dir, 'report.html'), payload)
    print(f'[report] -> {report_path}  ({len(payload["reads"])} barcodes, slider starts at '
          f'{payload["default_threshold"]:,} reads)')
    return report_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    ap = argparse.ArgumentParser(
        description='RT-Rep-Seq: On-rep-seq demux + Metabuli classify + LCP QC, one script.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input-dir', required=True, help='directory of ONT *.fastq(.gz) pass reads')
    ap.add_argument('--gtdb-db', required=True, help='Metabuli reference DB directory')
    ap.add_argument('--taxonomy-names', help='taxonomyDB_names.dmp '
                     '(default: <gtdb-db>/taxonomy/taxonomyDB_names.dmp)')
    ap.add_argument('--taxonomy-nodes', help='taxonomyDB_nodes.dmp '
                     '(default: <gtdb-db>/taxonomy/taxonomyDB_nodes.dmp)')
    ap.add_argument('--work-dir', default='work', help='intermediate files (default: ./work)')
    ap.add_argument('--results-dir', default='results', help='final outputs (default: ./results)')
    ap.add_argument('--num-barcodes', type=int, default=None,
                     help='how many barcodes were used; skips the interactive prompt')
    ap.add_argument('--min-reads', type=int, default=10000,
                     help='per-barcode read target for the demux stop rule (default: 10000)')
    ap.add_argument('--analysis-min-reads', type=int, default=1,
                     help='barcodes with at least this many reads are classified and get an LCP, '
                          'so the report slider can go below --min-reads (default: 1 = every '
                          'barcode seen). --min-reads stays the QC threshold: it drives the '
                          'demux stop rule, the static PNG figures and the slider start value')
    ap.add_argument('--min-read-length', type=int, default=500,
                     help='only reads >= this length (bp) are used for species classification and '
                          'LCP analysis (default: 500); read-length statistics use all reads')
    ap.add_argument('--report-min-reads', type=int, default=None,
                     help='start value of the report\'s read-threshold slider (default: --min-reads)')
    ap.add_argument('--run-name', default=None,
                     help='run name shown in the report header (default: derived from --input-dir)')
    ap.add_argument('--stop-fraction', type=float, default=0.9,
                     help='fraction of barcodes that must reach --min-reads (default: 0.9)')
    ap.add_argument('--min-species-pct', type=float, default=10.0,
                     help='zero out a species in species_by_barcode.tsv if it is below this %% '
                          'of that barcode\'s classified reads -- treated as noise (default: 10)')
    ap.add_argument('--lcp-bin-bp', type=int, default=10,
                     help='bin width, in bp, for the LCP histogram grid (default: 10)')
    ap.add_argument('--lcp-ymax', type=int, default=None,
                     help='y-axis max (read count) for the LCP histogram grid; default auto-sets '
                          'it to the highest single-bin count seen in any one sample')
    ap.add_argument('--lcp-xmax', type=int, default=None,
                     help='x-axis max (bp) for the LCP histogram grid; default auto-sets it to the '
                          'longest bin that reaches --lcp-min-count-for-xmax reads in any one '
                          'sample, so a few long outlier reads can\'t stretch the axis out and '
                          'squash the real signal into a sliver on the left')
    ap.add_argument('--lcp-min-count-for-xmax', type=int, default=20,
                     help='read-count threshold a bin must reach (in at least one sample) to count '
                          'toward the auto x-axis max (default: 20)')
    ap.add_argument('--lcp-cols-per-page', type=int, default=4,
                     help='max species columns per LCP histogram image before breaking onto a new '
                          'page (lcp_histograms_page1.png, _page2.png, ...); default: 4. Pass 0 for '
                          'a single unpaginated image.')
    ap.add_argument('--onrep-barcodes', default=os.path.join(HERE, 'assets', 'onrep_barcodes.fasta'))
    ap.add_argument('--onrep-config', default=os.path.join(HERE, 'assets', 'onrep_config.json'))
    ap.add_argument('--threads', type=int, default=None,
                     help='Metabuli threads (default: all cores of the CPU budget)')
    ap.add_argument('--metabuli-max-ram', type=int, default=None,
                     help='Metabuli --max-ram in GiB (default: 75%% of the RAM budget; Metabuli\'s '
                          'own default of 128 GiB would exceed the budget on most sequencing PCs)')
    ap.add_argument('--metabuli-bin', default='metabuli', help='metabuli executable if already on PATH')
    ap.add_argument('--metabuli-env', default=None,
                     help='conda env name to run metabuli in (via `conda run -n`)')
    ap.add_argument('--heartbeat-interval', type=int, default=15,
                     help='seconds between "still running" elapsed-time/resource-usage lines '
                          'while a long external command (e.g. Metabuli) runs (default: 15)')
    ap.add_argument('--logo-ucph', default=os.path.join(HERE, 'assets', 'ucph_logo.png'),
                     help='University of Copenhagen logo for the HTML report header')
    ap.add_argument('--logo-onrepseq', default=os.path.join(HERE, 'assets', 'onrepseq_logo.png'),
                     help='ON-rep-seq logo for the HTML report header')
    ap.add_argument('--skip', nargs='*', default=[],
                     choices=['demux', 'classify', 'lcps', 'heatmap', 'report'],
                     help='skip these stages and reuse existing --work-dir/--results-dir output')
    ap.add_argument('--watch', action='store_true',
                     help='live-run mode for a --input-dir that starts empty and fills in as '
                          'sequencing proceeds: wait for fastq(.gz) files, run the full pipeline, '
                          'and redo it (overwriting results) whenever new files have appeared '
                          'since the last completed pass. Runs until interrupted (Ctrl-C).')
    ap.add_argument('--poll-interval', type=int, default=20,
                     help='seconds between checks for new input files in --watch mode (default: 20)')
    ap.add_argument('--cpu-frac', type=float, default=0.5,
                     help='share of CPU cores this tool (and everything it spawns, e.g. '
                          'Metabuli) is pinned to (default: 0.5 -- half the machine)')
    ap.add_argument('--mem-frac', type=float, default=0.5,
                     help='share of system RAM this tool is capped to via a systemd-run '
                          'user scope (default: 0.5); skipped with a notice if systemd-run '
                          'is unavailable on this machine')
    ap.add_argument('--no-mem-cap', action='store_true',
                     help='apply the CPU cap only; never try the systemd-run RAM cap re-exec')
    ap.add_argument('--no-banner', action='store_true', help='skip the startup banner')
    args = ap.parse_args()

    if not args.taxonomy_names:
        args.taxonomy_names = os.path.join(args.gtdb_db, 'taxonomy', 'taxonomyDB_names.dmp')
    if not args.taxonomy_nodes:
        args.taxonomy_nodes = os.path.join(args.gtdb_db, 'taxonomy', 'taxonomyDB_nodes.dmp')
    return args


def run_pipeline(args):
    """One full demux -> classify -> lcps -> heatmap -> report pass."""
    run_start = time.time()
    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)
    print(f'[RT-Rep-Seq] starting  ({resource_snapshot()})')

    if 'demux' not in args.skip:
        passed, counts = step_demux(args)
    else:
        step_banner(1, 'Demultiplexing raw reads (torchlex) -- SKIPPED')
        counts = load_passed_from_counts(args.work_dir, min_reads=0)
        passed = {bc: n for bc, n in counts.items() if n >= args.min_reads}
        print(f'[demux] reusing a prior demux pass: {len(passed)} of {len(counts)} barcodes '
              f'>= {args.min_reads} reads')
    floor = max(1, min(args.analysis_min_reads, args.min_reads))
    analysed = {bc: n for bc, n in counts.items() if n >= floor}
    print(f'[RT-Rep-Seq] classifying + profiling {len(analysed)} barcodes with >= {floor} reads '
          f'({len(passed)} of them pass the {args.min_reads}-read QC threshold)')

    dominant = None
    if 'classify' not in args.skip:
        _, dominant = step_classify(args, analysed)
    else:
        step_banner(2, 'Classifying reads (Metabuli) -- SKIPPED')
    if dominant is None:
        dominant = load_dominant_from_table(args.results_dir)  # reuse a prior run's table, if any

    profiles = None
    if 'lcps' not in args.skip:
        profiles = step_lcps(args, analysed, dominant, plot_set=passed)
    else:
        step_banner(3, 'Building per-barcode length-count profiles (LCPs) -- SKIPPED')

    if 'heatmap' not in args.skip:
        if profiles is None:
            profiles = load_profiles(args, passed)
        step_heatmap(args, {bc: p for bc, p in profiles.items() if bc in passed}, dominant)
    else:
        step_banner(4, 'Clustering barcodes by LCP similarity (heatmap) -- SKIPPED')

    if 'report' not in args.skip:
        step_report(args, dominant)
    else:
        step_banner(5, 'Writing HTML report -- SKIPPED')

    elapsed = time.time() - run_start
    print(f'\n[RT-Rep-Seq] done in {elapsed / 60:.1f} min  ({resource_snapshot()}).')
    print(f'  {args.results_dir}/classification/species_by_barcode.tsv')
    print(f'  {args.results_dir}/LCPs/lcp_histograms.png (or _page1.png, _page2.png, ... if paginated)')
    print(f'  {args.results_dir}/LCPs/lcp_similarity_heatmap.png')
    print(f'  {args.results_dir}/report.html  (interactive: threshold slider, plates, downloads, A4 print)')


def _input_fastqs(input_dir):
    return set(glob.glob(os.path.join(input_dir, '*.fastq')) +
               glob.glob(os.path.join(input_dir, '*.fastq.gz')))


def run_watch(args):
    """Live-run mode for a MinKNOW output directory that starts empty and fills in as
    sequencing proceeds. Waits for fastq(.gz) files to appear, then repeatedly redoes
    the FULL pipeline (from scratch, over every file present at that moment) each time
    new files have shown up since the last completed pass -- never concurrently with
    a pass still in progress, since this loop is single-threaded and simply doesn't
    check again until run_pipeline() returns. Every pass wipes --work-dir first so
    reads from files already seen in a previous pass are never double-counted, and
    --results-dir so stale files from a since-shrunk output (e.g. fewer LCP histogram
    pages) never linger -- both are fully disposable/reproducible from the input
    fastqs, and neither lives inside the (read-only) input directory.
    """
    print(f'[watch] monitoring {args.input_dir} for *.fastq(.gz) -- checking every '
          f'{args.poll_interval}s until sequencing produces data')
    last_seen = frozenset()
    while True:
        current = _input_fastqs(args.input_dir)
        if not current or current == last_seen:
            time.sleep(args.poll_interval)
            continue

        n_new = len(current - last_seen)
        print(f'\n[watch] {len(current)} fastq(.gz) file(s) present ({n_new} new since the last '
              f'pass) -- starting a full re-analysis')
        last_seen = frozenset(current)

        for d in (args.work_dir, args.results_dir):
            shutil.rmtree(d, ignore_errors=True)

        try:
            run_pipeline(args)
        except Exception as e:
            print(f'[watch] analysis pass failed ({e}); will retry once more files arrive',
                  file=sys.stderr)
        # loop straight back to the top: if more files arrived while that pass ran,
        # `current` will already differ and the next pass starts immediately, no sleep


def main():
    # Line-buffer stdout even when it isn't a TTY (e.g. redirected to a log file or
    # piped), so status/progress output shows up as it happens instead of sitting in a
    # block buffer for however long a --watch run keeps the process alive.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass

    args = parse_args()

    global TOOL_RESOURCES
    TOOL_RESOURCES = Resources(args.cpu_frac, args.mem_frac)
    if not args.no_mem_cap:
        TOOL_RESOURCES.reexec_with_mem_cap(sys.argv)  # re-execs & restarts main(); never
        # returns if it actually re-executed -- everything below only runs once, in
        # whichever process (capped or not) ends up actually doing the work
    TOOL_RESOURCES.apply_cpu_cap()

    if not args.no_banner:
        banner()
    cap_note = '' if (args.no_mem_cap or TOOL_RESOURCES.mem_capped) else ' (RAM cap not enforced)'
    print(f'[RT-Rep-Seq] resource budget: CPU {len(TOOL_RESOURCES.cpus)}/{os.cpu_count()} '
          f'cores (cap {args.cpu_frac * 100:.0f}%), RAM cap '
          f'{human_bytes(TOOL_RESOURCES.mem_limit)}{cap_note}\n')

    if args.watch:
        if args.num_barcodes is None:  # ask once, up front, rather than once per pass
            args.num_barcodes = int(input('How many barcodes will be used in this run? ').strip())
        run_watch(args)
    else:
        run_pipeline(args)


if __name__ == '__main__':
    main()
