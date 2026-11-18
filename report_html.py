#!/usr/bin/env python3
"""report_html.py -- builds RT-Rep-Seq's single, self-contained, interactive report.html.

Reads only what is already in a results/ folder (plus, optionally, the demux read-count
table), embeds it as JSON, and lets the browser do the rest -- no server, works offline,
one file:

  * read-count threshold slider (log scale) + number box + quick-pick buttons
  * two 96-well plate maps (repBC01 -> plate 1 A1, repBC02 -> B1, ... column-wise;
    repBC97 -> plate 2 A1 ... repBC192 -> plate 2 H12), wells green when the barcode has
    >= threshold reads, red when below
  * species table, LCP similarity heatmap (correlation subset + average-linkage
    re-clustering) and LCP histogram grid, all recomputed live for the barcodes that
    pass the current threshold
  * downloads: tables as TSV, figures as PNG; print / save-as-PDF formatted for A4

Inputs (all optional except at least one of them):
  <results>/demux_read_counts.tsv           barcode, reads[, passed]   (reads per barcode)
  <results>/classification/species_by_barcode.tsv
  <results>/LCPs/raw_LCPs/<barcode>_LCP.txt  "count length" per line
"""
import base64
import datetime
import glob
import json
import os
import re
from collections import Counter

import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import lcp_peaks  # noqa: E402

PLATE_WELLS = 96


def _bc_num(bc):
    m = re.sub(r'\D', '', bc)
    return int(m) if m else 0


def _b64(path):
    if not path or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lstrip('.').lower()
    mime = {'png': 'image/png', 'jpg': 'image/jpeg', 'jpeg': 'image/jpeg',
            'svg': 'image/svg+xml'}.get(ext, 'application/octet-stream')
    with open(path, 'rb') as f:
        return f'data:{mime};base64,' + base64.b64encode(f.read()).decode('ascii')


def load_counts(path):
    counts = {}
    if path and os.path.exists(path):
        with open(path) as f:
            next(f, None)
            for line in f:
                c = line.rstrip('\n').split('\t')
                if len(c) >= 2 and c[1].strip().isdigit():
                    counts[c[0]] = int(c[1])
    return counts


def load_species(results_dir, name='species_by_barcode.tsv'):
    """barcode -> {species: reads} (non-zero only), from classification/<name>."""
    path = os.path.join(results_dir, 'classification', name)
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        barcodes = f.readline().rstrip('\n').split('\t')[1:]
        for bc in barcodes:
            out[bc] = {}
        for line in f:
            c = line.rstrip('\n').split('\t')
            for bc, v in zip(barcodes, c[1:]):
                n = int(v)
                if n:
                    out[bc][c[0]] = n
    return out


def load_profiles(results_dir):
    """barcode -> Counter{length: reads}, from LCPs/raw_LCPs/*_LCP.txt."""
    profiles = {}
    for path in glob.glob(os.path.join(results_dir, 'LCPs', 'raw_LCPs', '*_LCP.txt')):
        bc = os.path.basename(path)[:-len('_LCP.txt')]
        prof = Counter()
        with open(path) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    prof[int(parts[1])] = int(parts[0])
        profiles[bc] = prof
    return profiles


def correlation_matrix(profiles, barcodes):
    """Pearson r of relative exact-length distributions -- same math as
    rt_rep_seq.similarity_matrix_and_order, full (unclustered) matrix."""
    if not barcodes:
        return np.zeros((0, 0))
    lengths = sorted({ln for bc in barcodes for ln in profiles[bc]})
    idx = {ln: i for i, ln in enumerate(lengths)}
    mat = np.zeros((len(barcodes), len(lengths)))
    for r, bc in enumerate(barcodes):
        tot = sum(profiles[bc].values()) or 1
        for ln, n in profiles[bc].items():
            mat[r, idx[ln]] = n / tot
    if len(barcodes) < 2:
        return np.ones((1, 1))
    with np.errstate(invalid='ignore', divide='ignore'):
        sim = np.corrcoef(mat)
    return np.nan_to_num(sim, nan=0.0)


def length_stats(profile, min_len):
    n = sum(profile.values())
    if not n:
        return None
    mean = sum(L * c for L, c in profile.items()) / n
    sd = (sum(c * (L - mean) ** 2 for L, c in profile.items()) / n) ** 0.5
    return {'n': n, 'mean': round(mean, 1), 'sd': round(sd, 1),
            'n_min': sum(c for L, c in profile.items() if L >= min_len)}


def build_payload(results_dir, counts_path=None, default_threshold=10000, binsize=10,
                  meta=None, logo_ucph=None, logo_onrepseq=None, min_len=500):
    counts_path = counts_path or os.path.join(results_dir, 'demux_read_counts.tsv')
    counts = load_counts(counts_path)
    species = load_species(results_dir)
    species_all = load_species(results_dir, 'species_by_barcode_all.tsv')  # unfiltered, if present
    profiles = load_profiles(results_dir)

    reads_source = 'demultiplexed reads'
    if not counts:  # older results folder without the demux table -> best available proxy
        reads_source = 'classified reads (no demux read-count table found)'
        for bc in set(species) | set(profiles):
            counts[bc] = sum(species.get(bc, {}).values()) or sum(profiles.get(bc, {}).values())
    if not counts:
        raise SystemExit(f'[report] nothing to report under {results_dir}')

    max_num = max([_bc_num(bc) for bc in counts] + [PLATE_WELLS])
    n_plates = max(2, -(-max_num // PLATE_WELLS))

    # read-length statistics from ALL demultiplexed reads; everything below uses reads >= min_len
    lenstats = {bc: length_stats(p, min_len) for bc, p in profiles.items()}
    trimmed = {bc: {L: c for L, c in p.items() if L >= min_len} for bc, p in profiles.items()}
    trimmed = {bc: p for bc, p in trimmed.items() if p}
    lcp_bcs = sorted(trimmed, key=_bc_num)
    sim = correlation_matrix(trimmed, lcp_bcs)
    peaks = {bc: lcp_peaks.call_peaks(trimmed[bc], min_len=min_len) for bc in lcp_bcs}
    psim = lcp_peaks.similarity_matrix(peaks, lcp_bcs) if lcp_bcs else np.zeros((0, 0))
    bins = {}
    for bc in lcp_bcs:
        b = Counter()
        for ln, n in trimmed[bc].items():
            b[(ln // binsize) * binsize] += n
        bins[bc] = sorted(b.items())

    # species colour ranking: by total reads over all barcodes where it is dominant
    dom_tot = Counter()
    for bc, spc in species.items():
        if spc:
            sp, n = max(spc.items(), key=lambda kv: kv[1])
            dom_tot[sp] += counts.get(bc, n)

    meta = dict(meta or {})
    meta.setdefault('generated', datetime.datetime.now().strftime('%Y-%m-%d %H:%M'))
    meta['stamp'] = datetime.datetime.now().strftime('%Y%m%dT%H%M%S') + '-' + os.urandom(3).hex()
    meta['reads_source'] = reads_source

    return {
        'meta': meta,
        'default_threshold': int(default_threshold),
        'n_plates': n_plates,
        'plate_wells': PLATE_WELLS,
        'reads': counts,
        'species': species,
        'species_all': species_all,
        'species_rank': [sp for sp, _ in dom_tot.most_common()],
        'lcp_barcodes': lcp_bcs,
        'matrix': [[round(float(v), 3) for v in row] for row in sim],
        'peak_matrix': [[round(float(min(1.0, max(0.0, v))), 3) for v in row] for row in psim],
        'peaks': {bc: [[p[0], p[3]] for p in pk] for bc, pk in peaks.items()},   # [position bp, weight]
        'peak_params': lcp_peaks.PARAMS,
        'lenstats': lenstats,
        'min_len': min_len,
        'bins': bins,
        'binsize': binsize,
        'logo_ucph': _b64(logo_ucph),
        'logo_onrepseq': _b64(logo_onrepseq),
    }


def methods_meta(onrep_config=None, gtdb_db=None, metabuli_env=None, metabuli_bin=None,
                 params=None, tool_config=None, version=None):
    """Parameters, software versions and reference-DB details for the report's Materials
    and methods panel. Every lookup is best-effort: anything not found is simply left out."""
    import importlib.metadata as md
    import platform
    import shutil
    import subprocess

    out = {'params': dict(params or {}), 'versions': {}, 'db': {}}
    v = out['versions']
    v['python'] = platform.python_version()
    if version:
        v['rt_rep_seq'] = version
    for pkg in ('numpy', 'scipy', 'matplotlib', 'biopython', 'polyleven', 'torchlex'):
        try:
            v[pkg] = md.version(pkg)
        except Exception:
            pass
    try:
        if metabuli_env:
            r = subprocess.run(['conda', 'list', '-n', metabuli_env, '^metabuli$'],
                               capture_output=True, text=True, timeout=120)
            for line in r.stdout.splitlines():
                c = line.split()
                if c and c[0] == 'metabuli' and len(c) > 1:
                    v['metabuli'] = c[1]
        elif metabuli_bin:
            import glob as _g
            exe = shutil.which(metabuli_bin) or metabuli_bin
            prefix = os.path.dirname(os.path.dirname(os.path.realpath(exe)))
            for meta_json in _g.glob(os.path.join(prefix, 'conda-meta', 'metabuli-*.json')):
                v['metabuli'] = os.path.basename(meta_json).split('-')[1]
    except Exception:
        pass

    try:
        cfg = json.load(open(onrep_config)) if onrep_config else {}
        seq, amp = cfg.get('sequence', {}), cfg.get('amplicones', {})
        out['params'].update({
            'min_qual': seq.get('minimum_average_phred_quality'),
            'min_len': seq.get('min_sequence_length'), 'max_len': seq.get('max_sequence_length'),
            'window': seq.get('barcode_window_size'), 'max_distance': cfg.get('max_barcode_distance', 5),
            'barcode_length': amp.get('barcode_length'), 'spacer_length': amp.get('spacer_length'),
        })
    except Exception:
        pass

    if gtdb_db:
        db = out['db']
        db['name'] = os.path.basename(os.path.normpath(gtdb_db))
        try:
            for line in open(os.path.join(gtdb_db, 'db.parameters')):
                k, _, val = line.rstrip('\n').partition('\t')
                if k == 'Creation_date':
                    y, m, d = (val.strip().split('-') + ['', '', ''])[:3]
                    db['created'] = f'{y}-{m.zfill(2)}-{d.zfill(2)}' if d else val.strip()
        except Exception:
            pass
    tc = tool_config or {}
    if tc.get('db_description'):
        out['db']['description'] = tc['db_description']
    if tc.get('repo_url'):
        out['repo'] = tc['repo_url']
    return out


def load_tool_config(here):
    try:
        return json.load(open(os.path.join(here, 'assets', 'rt-rep-seq.config.json')))
    except Exception:
        return {}


def write_report(out_path, payload):
    html = HTML_TEMPLATE.replace('__PAYLOAD__', json.dumps(payload, separators=(',', ':'))
                                 .replace('</', '<\\/'))
    tmp = out_path + '.tmp'
    with open(tmp, 'w') as f:
        f.write(html)
    os.replace(tmp, out_path)              # never expose a half-written report
    # written LAST: an open report reloads itself when this stamp changes (live runs)
    stamp_path = os.path.join(os.path.dirname(os.path.abspath(out_path)), 'report_stamp.js')
    with open(stamp_path + '.tmp', 'w') as f:
        f.write('window.RTREPSEQ_STAMP = ' + json.dumps(payload.get('meta', {}).get('stamp', '')) + ';\n')
    os.replace(stamp_path + '.tmp', stamp_path)
    return out_path


HTML_TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RT-Rep-Seq report</title>
<style>
  :root { --red:#901a1e; --ok:#2e9d4f; --bad:#d0453b; --ink:#222; --muted:#666; --line:#ddd; }
  * { box-sizing: border-box; }
  body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: var(--ink);
         background: #fff; max-width: 1180px; margin: 0 auto; padding: 0 16px 40px; }
  header { display: flex; align-items: center; justify-content: space-between; gap: 1.5em;
           border-bottom: 3px solid var(--red); padding: 18px 0 14px; }
  header img { max-height: 64px; }
  header .title { flex: 1; }
  h1 { font-size: 1.45em; margin: 0 0 4px; }
  .credit, .meta { color: var(--muted); font-size: 0.88em; }
  .meta span { margin-right: 1.2em; white-space: nowrap; }
  .controls { position: sticky; top: 0; z-index: 20; background: #fff; border-bottom: 1px solid var(--line);
              padding: 10px 0; display: flex; flex-wrap: wrap; align-items: center; gap: 10px 16px; }
  .controls label { font-weight: 600; }
  #thrRange { width: 300px; max-width: 60vw; accent-color: var(--red); }
  #thrNum { width: 100px; padding: 3px 6px; font-size: 1em; }
  .chips button, .btn { border: 1px solid var(--line); background: #f7f7f7; border-radius: 4px;
                        padding: 4px 10px; cursor: pointer; font-size: 0.85em; }
  .chips button:hover, .btn:hover { background: #ececec; }
  .btn.primary { background: var(--red); color: #fff; border-color: var(--red); }
  .stats { display: flex; flex-wrap: wrap; gap: 10px; margin: 16px 0 4px; }
  .stat { border: 1px solid var(--line); border-radius: 6px; padding: 8px 14px; min-width: 150px; }
  .stat b { display: block; font-size: 1.35em; }
  .stat small { color: var(--muted); }
  h2 { font-size: 1.08em; margin: 0; }
  .hint { margin: 6px 0 0; }
  .grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; margin-top: 14px; }
  .card.wide { grid-column: 1 / -1; }
  @media (max-width: 900px) { .grid { grid-template-columns: minmax(0, 1fr); } }
  .card { border: 1px solid var(--line); border-radius: 10px; background: #fff; min-width: 0;
          box-shadow: 0 1px 2px rgba(0,0,0,.05); cursor: zoom-in; transition: box-shadow .15s, border-color .15s; }
  .card:hover { border-color: #cfa9aa; box-shadow: 0 4px 14px rgba(144,26,30,.13); }
  .card-head { display: flex; align-items: center; gap: 10px; padding: 10px 14px; border-bottom: 1px solid var(--line); }
  .card-head .sub, .ov-bar .sub { color: var(--muted); font-size: .85em; flex: 1; }
  .ico { width: 32px; height: 32px; border-radius: 8px; background: #f6e8e8; color: var(--red); flex: none;
         display: inline-flex; align-items: center; justify-content: center; }
  .ico svg { width: 21px; height: 21px; }
  .expand { border: 1px solid var(--line); background: #fff; border-radius: 6px; width: 30px; height: 30px;
            cursor: pointer; color: #555; display: inline-flex; align-items: center; justify-content: center; flex: none; }
  .expand svg { width: 16px; height: 16px; }
  .expand:hover { color: var(--red); border-color: var(--red); }
  .card-body { padding: 6px 14px 12px; max-height: 430px; overflow: hidden; position: relative; }
  .card-body::after { content: ''; position: absolute; left: 0; right: 0; bottom: 0; height: 64px;
                      background: linear-gradient(rgba(255,255,255,0), #fff 85%); pointer-events: none; }
  .card .table-scroll { max-height: none; overflow: hidden; }
  #hist.by-species { display: flex; gap: 18px; align-items: flex-start; overflow-x: auto; padding-bottom: 10px; }
  #hist.by-species canvas { max-width: none; flex: none; }
  .live { font-size: .85em; display: inline-flex; align-items: center; gap: 4px; white-space: nowrap; }
  .live[hidden] { display: none; }
  .live .dot { width: 8px; height: 8px; border-radius: 50%; background: var(--ok); display: inline-block; }
  .simctl { display: inline-flex; align-items: center; gap: 6px; padding-left: 10px; border-left: 1px solid var(--line); }
  .simctl label { font-weight: 600; }
  #simRange { width: 150px; accent-color: #2a6f97; }
  #simVal { min-width: 46px; color: #2a6f97; }
  .seg { display: flex; flex-wrap: wrap; align-items: center; gap: 6px 10px; margin: 4px 0 8px; font-size: .88em; }
  .seg > span { font-weight: 600; }
  .seg button { border: 1px solid var(--line); background: #f7f7f7; padding: 3px 10px; border-radius: 4px; cursor: pointer; font-size: 1em; }
  .seg button.on { background: #2a6f97; color: #fff; border-color: #2a6f97; }
  .seg select { font-size: 1em; padding: 2px 4px; }
  .note.warn { background: #fff8e6; border-left: 3px solid #e0a100; padding: 6px 10px; }
  .gchip { display: inline-block; min-width: 30px; text-align: center; padding: 1px 6px; border-radius: 10px; color: #fff;
           font-weight: 700; font-size: .85em; }
  .gchip.uniq { background: #9e9e9e; }
  .mem { display: inline-block; margin: 1px 6px 1px 0; white-space: nowrap; font-size: .92em; }
  #grpTable td { vertical-align: top; }
  #heatmaps canvas { margin-bottom: 10px; }
  .compcol { min-width: 220px; width: 30%; }
  .comp { display: flex; height: 14px; border-radius: 3px; overflow: hidden; background: #eee; min-width: 160px; }
  .comp span { display: block; height: 100%; border-right: 1px solid rgba(255,255,255,.85); cursor: default; }
  .comp span:last-child { border-right: none; }
  .comp span:hover { filter: brightness(0.85); }
  .btn.back { font-weight: 600; background: #fff; border-color: var(--red); color: var(--red); }
  .btn.back:hover { background: #f6e8e8; }
  .ov-foot { max-width: 1500px; margin: 28px auto 0; }
  .goback { width: 100%; border: 1px dashed #c9a3a4; background: #fbf5f5; color: var(--red); border-radius: 10px;
            padding: 16px; font-size: 1em; font-weight: 600; cursor: pointer; }
  .goback:hover { background: #f6e8e8; border-style: solid; }
  .goback .arr { font-size: 1.2em; margin-right: 6px; }
  #overlay { position: fixed; inset: 0; z-index: 40; background: #fff; display: flex; flex-direction: column; }
  #overlay[hidden] { display: none; }
  .ov-bar { display: flex; align-items: center; flex-wrap: wrap; gap: 8px 14px; padding: 8px 18px;
            border-bottom: 3px solid var(--red); background: #fff; }
  .ov-bar h2 { font-size: 1.2em; white-space: nowrap; }
  #ovCtl { flex: 1; min-width: 0; }
  #ovCtl .controls { position: static; border: none; padding: 0; }
  #ovCtl #printBtn { display: none; }
  .ov-body { flex: 1; overflow: auto; padding: 10px 24px 48px; }
  .ov-body .card-body { max-height: none; overflow: visible; padding: 0; max-width: 1500px; margin: 0 auto; }
  .ov-body .card-body::after { display: none; }
  .ov-body .table-scroll { max-height: none; overflow: visible; }
  body.ov-open { overflow: hidden; }
  .methods { font-family: Georgia, "Times New Roman", serif; font-size: 0.98em; line-height: 1.6; max-width: 860px; color: #1d1d1d; }
  .methods h3 { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; font-size: 0.98em; margin: 1.3em 0 0.3em; }
  .methods p { margin: 0 0 0.7em; text-align: justify; hyphens: auto; }
  .methods code { font-size: 0.88em; background: #f4f4f4; padding: 0 3px; border-radius: 3px; }
  .methods ol.refs { padding-left: 1.6em; font-size: 0.9em; }
  .methods ol.refs li { margin-bottom: 0.35em; }
  .methods a { color: #1f5f99; word-break: break-word; }
  .methods table { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; margin: 0.4em 0 1em; }
  .note { color: var(--muted); font-size: 0.9em; }
  canvas { display: block; max-width: 100%; height: auto; }
  .legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 0.85em; margin: 8px 0; }
  .legend i { display: inline-block; width: 12px; height: 12px; border-radius: 50%; vertical-align: -1px; margin-right: 5px; }
  .legend i.sq { border-radius: 2px; }
  .dl { display: flex; flex-wrap: wrap; gap: 8px; margin: 8px 0; }
  table { border-collapse: collapse; width: 100%; font-size: 0.88em; }
  th, td { padding: 4px 8px; border-bottom: 1px solid var(--line); text-align: left; }
  th { background: #f5f5f5; position: sticky; top: 0; }
  td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
  .table-scroll { max-height: 440px; overflow: auto; border: 1px solid var(--line); }
  #tip { position: fixed; pointer-events: none; background: rgba(30,30,30,0.92); color: #fff;
         font-size: 12px; padding: 6px 8px; border-radius: 4px; display: none; z-index: 50; line-height: 1.4; }
  .print-only { display: none; }
  footer { margin-top: 36px; padding-top: 10px; border-top: 1px solid var(--line); color: #888;
           font-size: 0.82em; text-align: center; }
  @media (max-width: 700px) { header { flex-wrap: wrap; } header img { max-height: 44px; } }
  @page { size: A4 portrait; margin: 12mm; }
  @media print {
    body { max-width: none; padding: 0; font-size: 10pt; }
    .no-print, #tip { display: none !important; }
    .print-only { display: block; }
    .controls { position: static; }
    .table-scroll { max-height: none; overflow: visible; border: none; }
    th { position: static; }
    #overlay { display: none !important; }
    .grid { display: block; }
    .card { border: none; box-shadow: none; margin: 0 0 10px; cursor: auto; }
    .card + .card { break-before: page; }
    .card-head { padding: 0 0 6px; border-bottom: 2px solid var(--red); margin-bottom: 6px; }
    .card-body { max-height: none; overflow: visible; padding: 0; }
    .card-body::after { display: none; }
    .card .table-scroll { overflow: visible; }
    .methods { font-size: 9.5pt; max-width: none; }
    .methods h3, .methods p { break-inside: avoid; }
    canvas { max-width: 100% !important; break-inside: avoid; }
    #plates { height: 158mm !important; width: auto !important; margin: 0 auto; }
    .stats { margin: 6px 0 0; } .stat b { font-size: 1.1em; }
    header img { max-height: 48px; }
    tr { break-inside: avoid; }
    .stat { padding: 4px 8px; }
  }
  * { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
</style></head>
<body>
<header>
  <img id="logoOn" alt="ON-rep-seq logo" hidden>
  <div class="title">
    <h1>RT-Rep-Seq run report</h1>
    <div class="credit">On-rep-seq method by Lukasz Krych &middot; Department of Food Science, University of Copenhagen</div>
    <div class="meta" id="metaLine"></div>
  </div>
  <img id="logoKU" alt="University of Copenhagen logo" hidden>
</header>

<div class="controls no-print" id="controls">
  <label for="thrRange">Min reads per barcode</label>
  <input type="range" id="thrRange" min="0" max="1000" step="1">
  <input type="number" id="thrNum" min="0" step="1">
  <span class="chips" id="chips"></span>
  <span class="simctl" title="Barcodes whose LCP band patterns are at least this similar (UPGMA, average linkage) form one strain group">
    <label for="simRange">Strain similarity &ge;</label>
    <input type="range" id="simRange" min="90" max="100" step="0.5" value="90"><b id="simVal">90%</b></span>
  <span style="flex:1"></span>
  <label class="live" id="liveBox" hidden title="Reload automatically when the analysis writes a newer report (keeps threshold, open panel and scroll position)">
    <input type="checkbox" id="liveChk" checked> Auto-update <span class="note" id="liveNote"></span></label>
  <button class="btn primary" id="printBtn">Print / PDF (A4)</button>
</div>
<div class="print-only" id="printThr"></div>

<div class="stats" id="stats"></div>
<p class="note no-print hint">Click any panel to open it as a full page (Esc or &larr; Go back to return). The read threshold and the strain-similarity slider update every panel.</p>

<main class="grid">
<section class="card wide" id="sec-plates">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><rect x="2.5" y="5" width="19" height="14" rx="2.5"/><circle cx="7" cy="9.5" r="1.4" fill="currentColor" stroke="none"/><circle cx="12" cy="9.5" r="1.4" fill="currentColor" stroke="none"/><circle cx="17" cy="9.5" r="1.4"/><circle cx="7" cy="14.5" r="1.4"/><circle cx="12" cy="14.5" r="1.4" fill="currentColor" stroke="none"/><circle cx="17" cy="14.5" r="1.4" fill="currentColor" stroke="none"/></svg></span><h2>Barcode plates</h2><span class="sub" id="sub-plates"></span><button class="expand no-print" title="Open full page" aria-label="Open Barcode plates full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
  <p class="note">repBC01&ndash;96 on plate 1 and repBC97&ndash;192 on plate 2, filled column-wise
    (A1, B1 &hellip; H1, A2 &hellip; H12). Each well shows the barcode number and its read count;
    hover a well for details.</p>
  <div class="seg no-print" id="plateMode"><span>Colour wells by</span>
    <button data-m="reads" class="on">read threshold</button><button data-m="strain">strain group</button></div>
  <div class="legend" id="plateLegend"></div>
  <canvas id="plates"></canvas>
  <div class="dl no-print"><button class="btn" data-png="plates">Plates (PNG)</button>
    <button class="btn" data-tsv="barcodes">Barcode summary (TSV)</button></div>
  </div>
</section>

<section class="card wide" id="sec-species">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6"><rect x="3" y="7.5" width="13" height="7" rx="3.5" transform="rotate(-25 9.5 11)"/><path d="M15.5 8.2c1.6-.9 2.4.6 3.6-.2s1.2-2.2 2.4-2.4"/><circle cx="8" cy="12" r=".9" fill="currentColor" stroke="none"/><circle cx="11" cy="10.4" r=".9" fill="currentColor" stroke="none"/></svg></span><h2>Species classification</h2><span class="sub" id="sub-species"></span><button class="expand no-print" title="Open full page" aria-label="Open Species classification full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
  <p class="note">Barcodes passing the threshold: dominant species from Metabuli and its share of the
    barcode's classified reads. <b>Composition</b> shows the relative distribution of all species classified
    in the barcode <span id="compSrc"></span>; point at a segment to see the species, its reads and share.
    Only reads &ge; <span class="minlen"></span> bp are classified; read length is given as mean &plusmn; SD of all demultiplexed reads.</p>
  <p class="note warn"><b>Note on per-read classification.</b> Metabuli assigns every ON-rep-seq read independently. A Rep-PCR
    amplicon is a random genomic fragment (often only 0.5&ndash;1.5 kb): fragments from conserved genes, mobile elements or plasmids
    shared between related taxa can be assigned to a neighbouring species, and very closely related species (e.g.
    <i>Lactococcus lactis</i> / <i>L. cremoris</i>) split each other's reads. The dominant species of a barcode is robust; minor species
    in the composition bar should be read as indications (possible contamination, mixed culture or shared DNA) rather than proof.
    Strain identity is judged from the LCP fingerprints, not from the species calls.</p>
  <div class="table-scroll"><table id="spTable"><thead><tr>
    <th>Barcode</th><th>Plate</th><th>Well</th><th class="num">Reads</th><th class="num">Read length (bp)</th><th>Dominant species</th>
    <th class="num">% of classified</th><th class="compcol">Composition (classified reads)</th></tr></thead><tbody></tbody></table></div>
  <div class="dl no-print"><button class="btn" data-tsv="species">Species &times; barcode (TSV)</button>
    <button class="btn" data-tsv="speciesAll" id="dlSpAll" hidden>All classified species, unfiltered (TSV)</button></div>
  </div>
</section>

<section class="card" id="sec-heatmap">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" stroke="none"><rect x="3" y="3" width="5.5" height="5.5" rx="1" fill="currentColor"/><rect x="9.25" y="3" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".45"/><rect x="15.5" y="3" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".2"/><rect x="3" y="9.25" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".45"/><rect x="9.25" y="9.25" width="5.5" height="5.5" rx="1" fill="currentColor"/><rect x="15.5" y="9.25" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".6"/><rect x="3" y="15.5" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".2"/><rect x="9.25" y="15.5" width="5.5" height="5.5" rx="1" fill="currentColor" opacity=".6"/><rect x="15.5" y="15.5" width="5.5" height="5.5" rx="1" fill="currentColor"/></svg></span><h2>LCP similarity heatmap</h2><span class="sub" id="sub-heatmap"></span><button class="expand no-print" title="Open full page" aria-label="Open LCP similarity heatmap full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
  <p class="note">LCP similarity between passing barcodes, one heatmap per dominant species, ordered by UPGMA
    clustering. Boxes on the diagonal are <b>strain groups</b>: barcodes whose band patterns are at least as similar as the
    <b>strain-similarity</b> slider (top bar). Similarity is based on the <b>presence/absence of LCP bands</b> (peaks), not
    their intensities.</p>
  <div class="seg no-print"><span>Similarity</span>
    <select id="simMetric"><option value="peaks" selected>band presence (weighted Dice, default)</option>
      <option value="pearson">whole-profile correlation (Pearson)</option></select>
    <label><input type="checkbox" id="bySpecies" checked> separate heatmap per species</label></div>
  <div id="heatmaps"></div>
  <div class="dl no-print"><button class="btn" data-png="heatmap">Heatmap (PNG)</button>
    <button class="btn" data-tsv="matrix">Similarity matrix (TSV)</button></div>
  </div>
</section>

<section class="card" id="sec-groups">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7" stroke-linecap="round" stroke-linejoin="round"><path d="M4 4v5h5M4 9v5h3M9 9v-5M7 14v6M7 14h4v6M11 9h6v11M17 9V4"/><circle cx="4" cy="4" r="1.2" fill="currentColor" stroke="none"/><circle cx="9" cy="4" r="1.2" fill="currentColor" stroke="none"/><circle cx="17" cy="4" r="1.2" fill="currentColor" stroke="none"/></svg></span><h2>Strain groups</h2><span class="sub" id="sub-groups"></span><button class="expand no-print" title="Open full page" aria-label="Open Strain groups full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
  <p class="note" id="groupSummary"></p>
  <div class="table-scroll"><table id="grpTable"><thead><tr><th>Group</th><th>Species</th><th class="num">n</th>
    <th>Barcodes (plate well)</th><th class="num">Similarity mean (min)</th><th>Shared bands (bp)</th></tr></thead><tbody></tbody></table></div>
  <div class="dl no-print"><button class="btn" data-tsv="groups">Strain groups (TSV)</button></div>
  </div>
</section>

<section class="card wide" id="sec-hist">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round"><path d="M3 20.5h18"/><path d="M5.5 20V15M8.5 20V6M11.5 20v-8M14.5 20V9M17.5 20v-3.5M20 20v-6" stroke-width="2.2"/></svg></span><h2>LCP histograms</h2><span class="sub" id="sub-hist"></span><button class="expand no-print" title="Open full page" aria-label="Open LCP histograms full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
  <p class="note">Read-length profiles (reads &ge; <span class="minlen"></span> bp) of the passing barcodes in heatmap order, i.e.
    grouped by species and strain group. Dashed vertical lines mark the <b>bands shared by the strain group</b> (present in &ge; half of
    its members), small triangles the bands detected in each barcode, and faint grey lines every 100 bp help to follow bands
    between samples. Full page: one column per species, strain groups stacked one under another.
    Y axis of every panel = 110% of that barcode's own highest peak (reads per <span id="binW"></span>-bp bin).</p>
  <div class="legend" id="spLegend2"></div>
  <div id="hist"></div>
  <div class="dl no-print"><button class="btn" data-png="hist">LCP histograms (PNG)</button>
    <button class="btn" data-tsv="lcp">Read-length profiles (TSV)</button></div>
  </div>
</section>

<section class="card wide" id="sec-methods">
  <div class="card-head"><span class="ico"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round"><path d="M9 3h6M10 3v6.2L4.6 18.3A1.8 1.8 0 0 0 6.2 21h11.6a1.8 1.8 0 0 0 1.6-2.7L14 9.2V3"/><path d="M7 15h10" /><circle cx="10.5" cy="17.6" r=".9" fill="currentColor" stroke="none"/><circle cx="13.8" cy="16.6" r=".7" fill="currentColor" stroke="none"/></svg></span><h2>Materials and methods</h2><span class="sub" id="sub-methods"></span><button class="expand no-print" title="Open full page" aria-label="Open Materials and methods full page"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7"/></svg></button></div>
  <div class="card-body">
    <div class="dl no-print"><button class="btn" id="copyMethods">Copy text</button>
      <button class="btn" id="dlMethods">Download (TXT)</button></div>
    <article class="methods" id="methodsText"></article>
  </div>
</section>
</main>

<div id="overlay" hidden>
  <div class="ov-bar"><button class="btn back" id="ovClose" title="Back to all panels (Esc)">&#8592; Go back</button>
    <span class="ico" id="ovIco"></span><h2 id="ovTitle"></h2><div id="ovCtl"></div></div>
  <div class="ov-body" id="ovBody"><div class="ov-foot" id="ovFoot">
    <button class="goback" id="ovBack2"><span class="arr">&#8592;</span> Go back to all panels</button></div></div>
</div>

<footer>Generated by RT-Rep-Seq <span id="ver"></span> &middot; <a id="repoLink" href="https://github.com/krychu84/RT-Rep-Seq">github.com/krychu84/RT-Rep-Seq</a>
  &middot; ON-rep-seq method developed by Lukasz Krych, Department of Food Science, University of Copenhagen</footer>
<div id="tip"></div>

<script>
const D = __PAYLOAD__;
const SCALE = 2;                       // canvas pixel density (crisp print + PNG export)
const OK = '#2e9d4f', BAD = '#d0453b';
const PAL = ['#3b6fa0','#e07b39','#4a9c5d','#b5485c','#7b5ea7','#c9a227','#3fa7b5','#8c6d4f','#d36fa8','#4f5d75'];
const ROWS = 'ABCDEFGH';

// ---------------------------------------------------------------- data helpers
const num = bc => parseInt(String(bc).replace(/\D/g, ''), 10) || 0;
const reads = bc => D.reads[bc] || 0;
const allBcs = Object.keys(D.reads).sort((a, b) => num(a) - num(b));
const prefix = (allBcs.find(b => /^\D+/.test(b)) || 'repBC').replace(/\d+$/, '');
const width2 = Math.max(2, String(Math.max(...allBcs.map(num), 99)).length);
const nameFor = n => allBcs.find(b => num(b) === n) || prefix + String(n).padStart(width2, '0');
function well(n) {               // 1-based barcode number -> {plate, well}
  const i = n - 1, p = Math.floor(i / D.plate_wells) + 1, k = i % D.plate_wells;
  return { plate: p, row: k % 8, col: Math.floor(k / 8), well: ROWS[k % 8] + (Math.floor(k / 8) + 1) };
}
function dominant(bc) {
  const s = D.species[bc]; if (!s) return null;
  const e = Object.entries(s).sort((a, b) => b[1] - a[1]);
  if (!e.length) return null;
  // dominant species = top call after the noise filter; its share is reported against ALL
  // classified (species-level) reads of the barcode, i.e. the unfiltered table when present
  const all = D.species_all && D.species_all[bc];
  if (all && Object.keys(all).length) {
    const ea = Object.entries(all).sort((a, b) => b[1] - a[1]);
    const tot = ea.reduce((t, x) => t + x[1], 0), n = all[e[0][0]] || e[0][1];
    const sec = ea.find(x => x[0] !== e[0][0]);
    return { sp: e[0][0], n, pct: 100 * n / tot, second: sec ? sec[0] : '', tot };
  }
  const tot = e.reduce((t, x) => t + x[1], 0);
  return { sp: e[0][0], n: e[0][1], pct: 100 * e[0][1] / tot, second: e[1] ? e[1][0] : '', tot };
}
const spColor = {};
D.species_rank.forEach((sp, i) => { spColor[sp] = i < PAL.length ? PAL[i] : '#9aa0a6'; });
const colorOf = bc => { const d = dominant(bc); return d ? (spColor[d.sp] || '#9aa0a6') : '#9aa0a6'; };
const fmtK = n => n >= 1e6 ? (n / 1e6).toFixed(1) + 'M' : n >= 1e4 ? Math.round(n / 1e3) + 'k'
                : n >= 1e3 ? (n / 1e3).toFixed(1) + 'k' : String(n);
const esc = s => String(s).replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));

function ctxFor(cv, w, h, scale) {
  const sc = scale || Math.min(SCALE, 30000 / Math.max(w, h));   // stay under browser canvas limits
  cv.width = Math.round(w * sc); cv.height = Math.round(h * sc);
  cv.style.width = w + 'px'; cv._w = w; cv._h = h;
  const ctx = cv.getContext('2d');
  ctx.setTransform(sc, 0, 0, sc, 0, 0);
  ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, w, h);
  return ctx;
}

// ---------------------------------------------------------------- threshold control
const maxReads = Math.max(1, ...allBcs.map(reads));
const rng = document.getElementById('thrRange'), box = document.getElementById('thrNum');
const posToVal = p => p <= 0 ? 0 : Math.round(Math.pow(10, (p / 1000) * Math.log10(maxReads + 1)) - 1);
const valToPos = v => v <= 0 ? 0 : Math.round(1000 * Math.log10(v + 1) / Math.log10(maxReads + 1));
box.max = maxReads;
let threshold = Math.min(D.default_threshold, maxReads);
const chips = document.getElementById('chips');
[100, 500, 1000, 2000, 5000, 10000].filter(v => v <= maxReads || v === D.default_threshold).forEach(v => {
  const b = document.createElement('button'); b.textContent = v.toLocaleString();
  b.onclick = () => setThreshold(v); chips.appendChild(b);
});
let pending = false;
function setThreshold(v, fromSlider) {
  threshold = Math.max(0, Math.round(+v || 0));
  box.value = threshold;
  if (!fromSlider) rng.value = valToPos(threshold);
  if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; renderAll(); }); }
}
rng.addEventListener('input', () => setThreshold(posToVal(+rng.value), true));
box.addEventListener('input', () => setThreshold(box.value));
const passing = () => allBcs.filter(bc => reads(bc) >= threshold);

// ---------------------------------------------------------------- header + stats
(function header() {
  if (D.logo_onrepseq) { const i = document.getElementById('logoOn'); i.src = D.logo_onrepseq; i.hidden = false; }
  if (D.logo_ucph) { const i = document.getElementById('logoKU'); i.src = D.logo_ucph; i.hidden = false; }
  const m = D.meta || {}, bits = [];
  if (m.run) bits.push('Run: <b>' + esc(m.run) + '</b>');
  if (m.input_dir) bits.push('Input: ' + esc(m.input_dir));
  if (m.files_read) bits.push(m.files_read + ' fastq files');
  bits.push('Generated ' + esc(m.generated || ''));
  document.getElementById('metaLine').innerHTML = bits.map(b => '<span>' + b + '</span>').join('');
  document.getElementById('ver').textContent = m.version ? 'v' + m.version : '';
  if (m.repo_url) { const a = document.getElementById('repoLink'); a.href = m.repo_url; a.textContent = m.repo_url.replace(/^https?:\/\//, ''); }
  if (m.repo_url) bits.push('<a href="' + esc(m.repo_url) + '">RT-Rep-Seq on GitHub</a>');
  document.getElementById('metaLine').innerHTML = bits.map(b => '<span>' + b + '</span>').join('');
})();

function renderStats(pass) {
  const tot = allBcs.reduce((t, b) => t + reads(b), 0);
  const inPass = pass.reduce((t, b) => t + reads(b), 0);
  const sp = new Set(pass.map(b => (dominant(b) || {}).sp).filter(Boolean));
  const seen = allBcs.filter(b => reads(b) > 0).length;
  const used = (D.meta || {}).num_barcodes;
  const cards = [
    [pass.length + ' / ' + seen, 'barcodes passing / with reads' + (used ? ' (' + used + ' used)' : '')],
    [inPass.toLocaleString(), 'reads in passing barcodes (' + (tot ? (100 * inPass / tot).toFixed(1) : 0) + '%)'],
    [tot.toLocaleString(), 'total ' + (D.meta.reads_source || 'reads')],
    [sp.size, 'dominant species among passing'],
  ];
  // pooled mean / SD of read length over the passing barcodes
  let N = 0, S1 = 0, S2 = 0;
  for (const b of pass) { const l = (D.lenstats || {})[b]; if (l) { N += l.n; S1 += l.n * l.mean; S2 += l.n * (l.sd * l.sd + l.mean * l.mean); } }
  if (N) { const m = S1 / N; cards.push([Math.round(m) + ' &plusmn; ' + Math.round(Math.sqrt(Math.max(0, S2 / N - m * m))) + ' bp', 'read length, mean &plusmn; SD (passing)']); }
  const multi = STR.groups.filter(g => g.size > 1);
  cards.push([multi.length + ' + ' + (STR.groups.length - multi.length), 'strain groups (&ge;2 barcodes) + unique, at &ge; ' + (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%']);
  document.getElementById('stats').innerHTML = cards.map(c => '<div class="stat"><b>' + c[0] + '</b><small>' + c[1] + '</small></div>').join('');
  document.getElementById('printThr').innerHTML = '<p><b>Threshold: &ge; ' + threshold.toLocaleString() +
    ' reads per barcode</b> &mdash; ' + pass.length + ' barcodes pass.</p>';
}

// ---------------------------------------------------------------- plates
const PL = { left: 26, top: 50, cell: 46, gap: 40, bottom: 10 };
let plateMode = 'reads';
PL.w = PL.left + 12 * PL.cell; PL.h = PL.top + 8 * PL.cell + PL.bottom;
let platesPerRow = 2;
function renderPlates() {
  const cv = document.getElementById('plates');
  const perRow = platesPerRow, rowsOfPlates = Math.ceil(D.n_plates / perRow);
  const W = perRow * PL.w + (perRow - 1) * PL.gap + 8, H = rowsOfPlates * PL.h + (rowsOfPlates - 1) * 20;
  const ctx = ctxFor(cv, W, H);
  cv._wells = [];
  for (let p = 0; p < D.n_plates; p++) {
    const ox = (p % perRow) * (PL.w + PL.gap), oy = Math.floor(p / perRow) * (PL.h + 20);
    const first = p * D.plate_wells + 1, last = first + D.plate_wells - 1;
    const pass = Array.from({ length: D.plate_wells }, (_, k) => reads(nameFor(first + k)) >= threshold).filter(Boolean).length;
    const strain = plateMode === 'strain';
    ctx.fillStyle = '#222'; ctx.font = 'bold 14px sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
    ctx.fillText('Plate ' + (p + 1) + '  (' + nameFor(first) + ' – ' + nameFor(last) + ')', ox, oy + 16);
    ctx.font = '12px sans-serif'; ctx.fillStyle = '#666'; ctx.textAlign = 'right';
    ctx.fillText(pass + ' / ' + D.plate_wells + ' pass', ox + PL.w, oy + 16);
    ctx.lineWidth = 1;
    // plate outline
    ctx.strokeStyle = '#bbb'; ctx.lineWidth = 1.2;
    roundRect(ctx, ox + PL.left - 6, oy + PL.top - 6, 12 * PL.cell + 12, 8 * PL.cell + 12, 10); ctx.stroke();
    ctx.fillStyle = '#555'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    for (let c = 0; c < 12; c++) ctx.fillText(c + 1, ox + PL.left + c * PL.cell + PL.cell / 2, oy + PL.top - 16);
    for (let r = 0; r < 8; r++) ctx.fillText(ROWS[r], ox + 10, oy + PL.top + r * PL.cell + PL.cell / 2);
    for (let k = 0; k < D.plate_wells; k++) {
      const n = first + k, bc = nameFor(n), w = well(n), rd = reads(bc);
      const cx = ox + PL.left + w.col * PL.cell + PL.cell / 2, cy = oy + PL.top + w.row * PL.cell + PL.cell / 2;
      const ok = rd >= threshold, g = STR.groupOf[bc];
      ctx.beginPath(); ctx.arc(cx, cy, PL.cell / 2 - 2.5, 0, 2 * Math.PI);
      if (strain) {
        if (ok && g) { ctx.fillStyle = g.color; ctx.fill(); }
        else { ctx.fillStyle = '#f3f3f3'; ctx.fill(); ctx.strokeStyle = '#cfcfcf'; ctx.lineWidth = 1; ctx.stroke(); }
        ctx.fillStyle = ok && g ? (g.size > 1 ? '#fff' : '#333') : '#aaa'; ctx.font = 'bold 12px sans-serif';
        ctx.fillText(String(n).padStart(width2, '0'), cx, cy - 7);
        ctx.font = (ok && g && g.size > 1 ? 'bold ' : '') + '10.5px sans-serif';
        ctx.fillText(ok && g ? g.label : (ok ? 'no LCP' : fmtK(rd)), cx, cy + 8);
      } else {
        ctx.fillStyle = ok ? OK : BAD; ctx.globalAlpha = rd > 0 ? 1 : 0.45; ctx.fill(); ctx.globalAlpha = 1;
        ctx.fillStyle = '#fff'; ctx.font = 'bold 12px sans-serif';
        ctx.fillText(String(n).padStart(width2, '0'), cx, cy - 7);
        ctx.font = '10.5px sans-serif'; ctx.fillText(fmtK(rd), cx, cy + 8);
      }
      cv._wells.push({ cx, cy, r: PL.cell / 2, bc, n, w, rd });
    }
  }
}
function renderPlateLegend() {
  const el = document.getElementById('plateLegend');
  if (plateMode === 'reads') {
    el.innerHTML = '<span><i style="background:var(--ok)"></i>&ge; read threshold (included)</span><span><i style="background:var(--bad)"></i>below threshold (excluded)</span>';
    return;
  }
  const multi = STR.groups.filter(g => g.size > 1);
  el.innerHTML = multi.slice(0, 40).map(g => '<span><i style="background:' + g.color + '"></i><b>' + g.label + '</b> <em>' + esc(shortSp(g.species)) + '</em> (' + g.size + ')</span>').join('') +
    (multi.length > 40 ? '<span>&hellip; ' + (multi.length - 40) + ' more</span>' : '') +
    '<span><i style="background:#c9c9c9"></i>unique strain</span><span><i style="background:#f3f3f3;border:1px solid #ccc"></i>below read threshold</span>';
}
document.querySelectorAll('#plateMode button').forEach(b => b.addEventListener('click', e => {
  e.stopPropagation(); plateMode = b.dataset.m;
  document.querySelectorAll('#plateMode button').forEach(x => x.classList.toggle('on', x === b));
  renderPlates(); renderPlateLegend(); saveState();
}));
function roundRect(ctx, x, y, w, h, r) {
  ctx.beginPath(); ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
  ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
}
const tip = document.getElementById('tip');
function canvasPoint(cv, e) { const b = cv.getBoundingClientRect(); return [(e.clientX - b.left) * cv._w / b.width, (e.clientY - b.top) * cv._h / b.height]; }
function showTip(e, html) { tip.innerHTML = html; tip.style.display = 'block';
  tip.style.left = Math.min(e.clientX + 14, innerWidth - tip.offsetWidth - 8) + 'px'; tip.style.top = (e.clientY + 14) + 'px'; }
document.getElementById('plates').addEventListener('mousemove', e => {
  const cv = e.currentTarget, [x, y] = canvasPoint(cv, e);
  const h = (cv._wells || []).find(o => (x - o.cx) ** 2 + (y - o.cy) ** 2 <= o.r * o.r);
  if (!h) { tip.style.display = 'none'; return; }
  const d = dominant(h.bc), g = STR.groupOf[h.bc], ls = (D.lenstats || {})[h.bc];
  showTip(e, '<b>' + esc(h.bc) + '</b> &middot; plate ' + h.w.plate + ' ' + h.w.well + '<br>' + h.rd.toLocaleString() + ' reads &middot; ' +
    (h.rd >= threshold ? '<span style="color:#7fe39a">passes</span>' : '<span style="color:#ff9b93">below threshold</span>') +
    (ls ? '<br>read length ' + Math.round(ls.mean) + ' &plusmn; ' + Math.round(ls.sd) + ' bp' : '') +
    (d ? '<br><i>' + esc(d.sp) + '</i> (' + d.pct.toFixed(0) + '%)' : '') +
    (g && h.rd >= threshold ? '<br>strain group <b>' + g.label + '</b>' + (g.size > 1 ? ' (' + g.size + ' barcodes)' : ' (unique)') : ''));
});
document.getElementById('plates').addEventListener('mouseleave', () => { tip.style.display = 'none'; });

// ---------------------------------------------------------------- species table
const HAS_ALL = D.species_all && Object.keys(D.species_all).length > 0;
const EXTRA = ['#9ecae1','#fdae6b','#a1d99b','#fc9272','#bcbddc','#e6c86e','#8fd3d9','#c9a98b','#f1a7cf','#a7b0c0',
               '#6baed6','#fd8d3c','#74c476','#ef6548','#9e9ac8','#d4b106','#41b6c4','#a07850','#dd6fae','#7d8798'];
const segColorMap = {};
(function assignColors() {
  const tot = {};
  for (const src of [D.species_all || {}, D.species]) for (const bc in src) for (const sp in src[bc]) tot[sp] = (tot[sp] || 0) + src[bc][sp];
  let k = 0;
  Object.keys(tot).sort((a, b) => tot[b] - tot[a]).forEach(sp => {
    segColorMap[sp] = spColor[sp] && spColor[sp] !== '#9aa0a6' ? spColor[sp] : EXTRA[k++ % EXTRA.length];
  });
})();
const OTHER_PCT = 1;   // species below 1% of a barcode are merged into one "other species" segment
function composition(bc) {
  const src = (HAS_ALL ? D.species_all[bc] : null) || D.species[bc] || {};
  const e = Object.entries(src).filter(x => x[1] > 0).sort((a, b) => b[1] - a[1]);
  const tot = e.reduce((t, x) => t + x[1], 0);
  if (!tot) return { tot: 0, segs: [] };
  const segs = [], rest = [];
  for (const [sp, n] of e) (100 * n / tot >= OTHER_PCT ? segs : rest).push({ sp, n, pct: 100 * n / tot });
  if (rest.length) segs.push({ sp: rest.length + ' other species (<' + OTHER_PCT + '% each)', n: rest.reduce((t, x) => t + x.n, 0),
                               pct: rest.reduce((t, x) => t + x.pct, 0), other: rest });
  return { tot, segs };
}
function compBar(bc) {
  const c = composition(bc);
  if (!c.tot) return '<span class="note">no classified reads</span>';
  return '<div class="comp" data-bc="' + esc(bc) + '">' + c.segs.map((g, i) =>
    '<span data-i="' + i + '" style="width:' + g.pct.toFixed(3) + '%;background:' + (g.other ? '#d0d0d0' : (segColorMap[g.sp] || '#9aa0a6')) + '"></span>').join('') + '</div>';
}
document.getElementById('compSrc').textContent = HAS_ALL
  ? '(all Metabuli species-level calls, before the ' + (((D.meta.methods || {}).params || {}).min_species_pct || 10) + '% noise filter used for the dominant-species call)'
  : '(species above the noise filter only; rerun with RT-Rep-Seq ≥ 1.1 for the unfiltered composition)';
if (HAS_ALL) document.getElementById('dlSpAll').hidden = false;
document.getElementById('spTable').addEventListener('mousemove', e => {
  const seg = e.target.closest('.comp span'); if (!seg) { tip.style.display = 'none'; return; }
  const bc = seg.parentElement.dataset.bc, c = composition(bc), g = c.segs[+seg.dataset.i]; if (!g) return;
  let html = '<b>' + esc(bc) + '</b> &middot; ' + c.tot.toLocaleString() + ' classified reads<br>' +
    '<span style="display:inline-block;width:9px;height:9px;background:' + seg.style.background + ';margin-right:5px"></span>' +
    (g.other ? esc(g.sp) : '<i>' + esc(g.sp) + '</i>') + ': ' + g.n.toLocaleString() + ' reads (' + g.pct.toFixed(1) + '%)';
  if (g.other) html += '<br>' + g.other.slice(0, 12).map(o => '&nbsp;&nbsp;<i>' + esc(o.sp) + '</i> ' + o.n.toLocaleString() + ' (' + o.pct.toFixed(2) + '%)').join('<br>') +
    (g.other.length > 12 ? '<br>&nbsp;&nbsp;&hellip; and ' + (g.other.length - 12) + ' more' : '');
  showTip(e, html);
});
document.getElementById('spTable').addEventListener('mouseleave', () => { tip.style.display = 'none'; });
function renderTable(pass) {
  const tb = document.querySelector('#spTable tbody');
  tb.innerHTML = pass.slice().sort((a, b) => num(a) - num(b)).map(bc => {
    const d = dominant(bc), w = well(num(bc));
    const ls = (D.lenstats || {})[bc];
    return '<tr><td>' + esc(bc) + '</td><td>' + w.plate + '</td><td>' + w.well + '</td><td class="num">' + reads(bc).toLocaleString() +
      '</td><td class="num">' + (ls ? Math.round(ls.mean) + ' &plusmn; ' + Math.round(ls.sd) : '') +
      '</td><td>' + (d ? '<i style="color:' + (spColor[d.sp] || '#555') + '">&#9679;</i> <i>' + esc(d.sp) + '</i>' : '<span class="note">not classified</span>') +
      '</td><td class="num">' + (d ? d.pct.toFixed(1) : '') + '</td><td>' + compBar(bc) + '</td></tr>';
  }).join('') || '<tr><td colspan="8" class="note">no barcode passes this threshold</td></tr>';
}
function renderLegends(pass) {
  const sps = [...new Set(pass.map(b => (dominant(b) || {}).sp).filter(Boolean))]
    .sort((a, b) => D.species_rank.indexOf(a) - D.species_rank.indexOf(b));
  const html = sps.map(sp => '<span><i class="sq" style="background:' + (spColor[sp] || '#9aa0a6') + '"></i><em>' + esc(sp) + '</em></span>').join('');
  ['spLegend1', 'spLegend2'].forEach(id => { const el = document.getElementById(id); if (el) el.innerHTML = html; });
}

// ---------------------------------------------------------------- strain grouping
const lcpIndex = {}; D.lcp_barcodes.forEach((b, i) => { lcpIndex[b] = i; });
const PP = D.peak_params || { tol_bp: 20, tol_frac: 0.02 };
const tolAt = pos => Math.max(PP.tol_bp, PP.tol_frac * pos);
let simThr = 0.90, simMetric = 'peaks', bySpeciesHM = true;
const simOf = (a, b) => { const i = lcpIndex[a], j = lcpIndex[b]; return simMetric === 'pearson' ? D.matrix[i][j] : D.peak_matrix[i][j]; };
const GPAL = ['#1f77b4','#d62728','#2ca02c','#9467bd','#ff7f0e','#17becf','#8c564b','#e377c2','#bcbd22','#393b79',
              '#ad494a','#637939','#7b4173','#3182bd','#e6550d','#31a354','#756bb1','#843c39','#5254a3','#8c6d31'];
const shortSp = sp => { const w = String(sp || '').split(' '); return w.length > 1 ? w[0][0] + '. ' + w.slice(1).join(' ') : sp; };
function upgma(bcs) {
  // average linkage on 1 - similarity; returns full leaf order + clusters formed at distance <= 1 - simThr
  const n = bcs.length;
  if (n === 1) return { order: bcs.slice(), clusters: [bcs.slice()] };
  const d = bcs.map(a => Float64Array.from(bcs, b => a === b ? 0 : 1 - simOf(a, b)));
  const size = new Array(n).fill(1), order = bcs.map(b => [b]), alive = new Array(n).fill(true);
  const cut = 1 - simThr + 1e-9; let clusters = null;
  for (let step = 0; step < n - 1; step++) {
    let bi = -1, bj = -1, best = Infinity;
    for (let i = 0; i < n; i++) if (alive[i]) for (let j = i + 1; j < n; j++) if (alive[j] && d[i][j] < best) { best = d[i][j]; bi = i; bj = j; }
    if (!clusters && best > cut) clusters = order.filter((_, k) => alive[k]).map(o => o.slice());
    const si = size[bi], sj = size[bj];
    for (let k = 0; k < n; k++) if (alive[k] && k !== bi && k !== bj) { const v = (si * d[bi][k] + sj * d[bj][k]) / (si + sj); d[bi][k] = v; d[k][bi] = v; }
    size[bi] = si + sj; order[bi] = order[bi].concat(order[bj]); alive[bj] = false;
  }
  const full = order[alive.indexOf(true)];
  if (!clusters) clusters = [full.slice()];
  const pos = {}; full.forEach((b, i) => { pos[b] = i; });
  clusters.forEach(c => c.sort((x, y) => pos[x] - pos[y]));
  clusters.sort((x, y) => pos[x[0]] - pos[y[0]]);
  return { order: full, clusters };
}
function consensusBands(members) {
  // bands present in >= half of the members (weights ignored), position = median
  const pts = [];
  members.forEach((b, k) => (D.peaks[b] || []).forEach(([p, w]) => pts.push([p, k])));
  pts.sort((x, y) => x[0] - y[0]);
  const out = []; let cur = [];
  const flush = () => { if (!cur.length) return; const who = new Set(cur.map(x => x[1]));
    if (who.size >= Math.max(1, Math.ceil(members.length / 2)) && (members.length > 1 || true)) {
      const ps = cur.map(x => x[0]).sort((x, y) => x - y); out.push(ps[Math.floor(ps.length / 2)]); } cur = []; };
  for (const pt of pts) { if (cur.length && pt[0] - cur[0][0] > tolAt(cur[0][0])) flush(); cur.push(pt); }
  flush();
  return out;
}
let STR = { groups: [], groupOf: {}, order: [], blocks: [] };
function computeStrains(pass) {
  const bcs = pass.filter(b => b in lcpIndex);
  const parts = new Map();
  for (const b of bcs) { const k = bySpeciesHM ? ((dominant(b) || {}).sp || 'unclassified') : 'All passing barcodes'; if (!parts.has(k)) parts.set(k, []); parts.get(k).push(b); }
  const rank = sp => { const i = D.species_rank.indexOf(sp); return i < 0 ? 1e9 : i; };
  const keys = [...parts.keys()].sort((x, y) => (parts.get(y).length - parts.get(x).length) || (rank(x) - rank(y)));
  const groups = [], groupOf = {}, order = [], blocks = [];
  let gi = 0, ci = 0;
  for (const k of keys) {
    const u = upgma(parts.get(k)), blockGroups = [];
    for (const members of u.clusters) {
      let sum = 0, cnt = 0, mn = 1;
      for (let i = 0; i < members.length; i++) for (let j = i + 1; j < members.length; j++) { const v = simOf(members[i], members[j]); sum += v; cnt++; mn = Math.min(mn, v); }
      const sp = bySpeciesHM ? k : ((dominant(members[0]) || {}).sp || 'unclassified');
      const g = { id: ++gi, label: 'G' + gi, species: sp, members, size: members.length, mean: cnt ? sum / cnt : 1, min: cnt ? mn : 1,
                  color: members.length > 1 ? GPAL[ci++ % GPAL.length] : '#c9c9c9', bands: consensusBands(members) };
      groups.push(g); blockGroups.push(g); members.forEach(b => { groupOf[b] = g; order.push(b); });
    }
    blocks.push({ key: k, members: blockGroups.flatMap(g => g.members), groups: blockGroups });
  }
  STR = { groups, groupOf, order, blocks };
  currentOrder = order;
}
function viridis(t) {
  const s = [[68,1,84],[72,36,117],[65,68,135],[53,95,141],[42,120,142],[33,145,140],[34,168,132],[68,191,112],[122,209,81],[189,223,38],[253,231,37]];
  t = Math.min(1, Math.max(0, t)); const p = t * (s.length - 1), i = Math.floor(p), f = p - i, a = s[i], b = s[Math.min(i + 1, s.length - 1)];
  return 'rgb(' + a.map((v, k) => Math.round(v + (b[k] - v) * f)).join(',') + ')';
}
const SIM_LO = 0.5;                       // colour scale 0.5 .. 1
const simColor = v => viridis((v - SIM_LO) / (1 - SIM_LO));

// ---------------------------------------------------------------- heatmaps (one per species, strain groups boxed)
let currentOrder = [];
function renderHeatmap() {
  const host = document.getElementById('heatmaps'); host.innerHTML = '';
  if (!STR.order.length) { host.innerHTML = '<p class="note">No passing barcode with read-length data.</p>'; return; }
  for (const blk of STR.blocks) {
    const order = blk.members, n = order.length;
    const cell = Math.max(5, Math.min(24, Math.floor(760 / Math.max(n, 1))));
    const font = Math.max(6, Math.min(11, cell - 2));
    const left = Math.ceil(Math.max(...order.map(b => (b + '  G' + (STR.groupOf[b] || {}).id).length)) * font * 0.62) + 16;
    const top = 30, grid = n * cell, colLab = Math.ceil(Math.max(...order.map(b => b.length)) * font * 0.62) + 8;
    const W = Math.max(left + grid + 110, 420), H = top + grid + colLab + 8;
    const cv = document.createElement('canvas'); host.appendChild(cv);
    const ctx = ctxFor(cv, W, H);
    // title
    const ng = blk.groups.filter(g => g.size > 1).length, nu = blk.groups.length - ng;
    ctx.textBaseline = 'middle'; ctx.textAlign = 'left';
    if (bySpeciesHM) { ctx.fillStyle = spColor[blk.key] || '#9aa0a6'; ctx.fillRect(0, 9, 12, 12); }
    ctx.fillStyle = '#222'; ctx.font = (bySpeciesHM ? 'italic ' : '') + 'bold 13px sans-serif';
    ctx.fillText(blk.key, bySpeciesHM ? 18 : 0, 15);
    const tw = ctx.measureText(blk.key).width;
    ctx.font = '11px sans-serif'; ctx.fillStyle = '#666';
    ctx.fillText('  ' + n + ' barcode' + (n > 1 ? 's' : '') + ' · ' + ng + ' strain group' + (ng === 1 ? '' : 's') + ' + ' + nu + ' unique at ≥ ' + (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%', (bySpeciesHM ? 18 : 0) + tw, 15);
    // cells
    for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {
      ctx.fillStyle = r === c ? simColor(1) : simColor(simOf(order[r], order[c]));
      ctx.fillRect(left + c * cell, top + r * cell, cell + 0.5, cell + 0.5);
    }
    // strain-group boxes on the diagonal
    let at = 0;
    for (const g of blk.groups) {
      const x = left + at * cell, sz = g.size * cell;
      if (g.size > 1) {
        ctx.strokeStyle = '#fff'; ctx.lineWidth = 4; ctx.strokeRect(x, top + at * cell, sz, sz);
        ctx.strokeStyle = g.color; ctx.lineWidth = 2.2; ctx.strokeRect(x, top + at * cell, sz, sz);
      }
      at += g.size;
    }
    // row labels with group chip
    ctx.font = font + 'px sans-serif';
    order.forEach((b, r) => {
      const g = STR.groupOf[b], y = top + r * cell + cell / 2;
      ctx.fillStyle = g.color; ctx.fillRect(left - 7, top + r * cell + 0.5, 5, cell - 1);
      ctx.textAlign = 'right'; ctx.fillStyle = '#222'; ctx.fillText(b, left - 10 - ctx.measureText('  G' + g.id).width, y);
      ctx.fillStyle = g.size > 1 ? g.color : '#999'; ctx.font = 'bold ' + font + 'px sans-serif'; ctx.fillText('G' + g.id, left - 10, y); ctx.font = font + 'px sans-serif';
    });
    order.forEach((b, c) => {
      ctx.save(); ctx.translate(left + c * cell + cell / 2, top + grid + 5); ctx.rotate(Math.PI / 2);
      ctx.fillStyle = '#222'; ctx.textAlign = 'left'; ctx.fillText(b, 0, 0); ctx.restore();
    });
    // colour bar with the strain threshold marked
    const bx = left + grid + 24, bh = Math.max(80, Math.min(grid, 220)), barW = 14;
    for (let i = 0; i < bh; i++) { ctx.fillStyle = viridis(1 - i / bh); ctx.fillRect(bx, top + i, barW, 1.5); }
    ctx.strokeStyle = '#999'; ctx.lineWidth = 1; ctx.strokeRect(bx, top, barW, bh);
    ctx.fillStyle = '#333'; ctx.font = '10px sans-serif'; ctx.textAlign = 'left';
    [[1, 0], [0.75, 0.5], [0.5, 1]].forEach(([v, f]) => ctx.fillText(v.toFixed(2).replace(/0$/, ''), bx + barW + 4, top + f * bh));
    const ty = top + (1 - (simThr - SIM_LO) / (1 - SIM_LO)) * bh;
    ctx.strokeStyle = '#d62728'; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(bx - 4, ty); ctx.lineTo(bx + barW + 4, ty); ctx.stroke();
    ctx.fillStyle = '#d62728'; ctx.font = 'bold 10px sans-serif'; ctx.textAlign = 'right'; ctx.fillText('≥' + Math.round(simThr * 1000) / 10 + '%', bx - 6, ty);
    cv._cells = { left, top, cell, n, order };
    cv.addEventListener('mousemove', e => {
      const g = cv._cells, [x, y] = canvasPoint(cv, e), c = Math.floor((x - g.left) / g.cell), r = Math.floor((y - g.top) / g.cell);
      if (c < 0 || r < 0 || c >= g.n || r >= g.n) { tip.style.display = 'none'; return; }
      const a = g.order[r], b = g.order[c], ga = STR.groupOf[a], gb = STR.groupOf[b];
      showTip(e, esc(a) + ' (G' + ga.id + ') &times; ' + esc(b) + ' (G' + gb.id + ')<br>similarity ' + (100 * simOf(a, b)).toFixed(1) + '%' +
        (ga === gb ? ' &middot; same strain group' : ''));
    });
    cv.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
  }
}
document.getElementById('simMetric').addEventListener('change', e => { simMetric = e.target.value; renderAll(); saveState(); });
document.getElementById('bySpecies').addEventListener('change', e => { bySpeciesHM = e.target.checked; renderAll(); saveState(); });
['simMetric', 'bySpecies'].forEach(id => document.getElementById(id).addEventListener('click', e => e.stopPropagation()));

// ---------------------------------------------------------------- strain groups table
function renderGroups() {
  const multi = STR.groups.filter(g => g.size > 1), uniq = STR.groups.filter(g => g.size === 1);
  const thr = (simThr * 100).toFixed(1).replace(/\.0$/, '');
  document.getElementById('groupSummary').innerHTML = STR.order.length
    ? 'At &ge; <b>' + thr + '%</b> band-pattern similarity the ' + STR.order.length + ' passing barcodes form <b>' + multi.length +
      ' strain group' + (multi.length === 1 ? '' : 's') + '</b> with &ge; 2 barcodes (' + multi.reduce((t, g) => t + g.size, 0) +
      ' barcodes: replicates or isolates of the same strain) and <b>' + uniq.length + ' unique</b> profile' + (uniq.length === 1 ? '' : 's') +
      '. Move the strain-similarity slider to split (higher %) or merge (lower %) groups.'
    : 'No passing barcode with read-length data.';
  const row = g => '<tr><td><span class="gchip' + (g.size > 1 ? '' : ' uniq') + '" style="' + (g.size > 1 ? 'background:' + g.color : '') + '">' + g.label + '</span></td>' +
    '<td><i>' + esc(g.species) + '</i></td><td class="num">' + g.size + '</td><td>' +
    g.members.map(b => { const w = well(num(b)); return '<span class="mem">' + esc(b) + ' <span class="note">' + w.plate + '&thinsp;' + w.well + '</span></span>'; }).join(' ') +
    '</td><td class="num">' + (g.size > 1 ? (100 * g.mean).toFixed(1) + ' (' + (100 * g.min).toFixed(1) + ')' : '&ndash;') +
    '</td><td class="note">' + g.bands.join(', ') + '</td></tr>';
  document.querySelector('#grpTable tbody').innerHTML = multi.map(row).join('') + uniq.map(row).join('') ||
    '<tr><td colspan="6" class="note">no groups</td></tr>';
}

// ---------------------------------------------------------------- LCP histograms
const HC = { cols: 6, pw: 176, ph: 88, gapX: 8, gapY: 30, rowsPerCanvas: 8, left: 30,   // dashboard grid
             colW: 330, colPh: 92, colGap: 26, head: 30 };                            // full page, by species
function niceStep(max) { const r = max / 4, p = Math.pow(10, Math.floor(Math.log10(r))); return [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= r) || p * 10; }
function peakOf(bc, xmax) { let m = 0; for (const [s, n] of D.bins[bc] || []) if (s < xmax && n > m) m = n; return m; }
const XMIN = D.min_len || 0;
function histXmax(bcs) {
  let xmax = 0;
  for (const b of bcs) for (const [s, n] of D.bins[b]) if (n >= 20) xmax = Math.max(xmax, s);
  if (!xmax) for (const b of bcs) for (const [s] of D.bins[b]) xmax = Math.max(xmax, s);
  return xmax + D.binsize;
}
function speciesGroups(bcs) {                    // Map species -> barcodes (input order kept = heatmap order)
  const g = new Map();
  for (const b of bcs) { const sp = (dominant(b) || {}).sp || 'unclassified'; if (!g.has(sp)) g.set(sp, []); g.get(sp).push(b); }
  return g;
}
function yLimits(groups, xmax) {                 // barcode -> y max = 110% of its own highest peak
  const lim = new Map();
  for (const members of groups.values()) for (const b of members) lim.set(b, Math.max(2, Math.ceil(Math.max(1, peakOf(b, xmax)) * 1.1)));
  return lim;
}
function drawPanel(ctx, bc, x0, y0, w, h, xmax, ymax, xs, label2, opts) {
  const L = HC.left, pw = w - L, X = v => x0 + L + 1 + (v - XMIN) / (xmax - XMIN) * (pw - 2);   // plot area: x0+L .. x0+w
  const g = STR.groupOf[bc];
  ctx.fillStyle = '#222'; ctx.font = 'bold 10px sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'alphabetic';
  ctx.fillText(bc, x0 + L, y0 - 4);
  let lx = x0 + L + ctx.measureText(bc).width + 6;
  if (g) { ctx.fillStyle = g.size > 1 ? g.color : '#999'; ctx.fillText(g.label, lx, y0 - 4); }
  ctx.font = 'italic 9px sans-serif'; ctx.fillStyle = '#555'; ctx.textAlign = 'right';
  ctx.fillText(label2, x0 + w, y0 - 4, pw - 80);
  // faint 100-bp reference grid
  ctx.fillStyle = '#efefef';
  for (let t = Math.ceil(XMIN / 100) * 100; t < xmax; t += 100) ctx.fillRect(X(t), y0 + 1, 0.7, h - 2);
  ctx.strokeStyle = '#c8c8c8'; ctx.lineWidth = 1; ctx.strokeRect(x0 + L + 0.5, y0 + 0.5, pw - 1, h - 1);
  if (g && g.size > 1) { ctx.fillStyle = g.color; ctx.fillRect(x0 + L - 4, y0, 3, h); }
  // y axis: 0, half, max
  ctx.fillStyle = '#777'; ctx.font = '8px sans-serif'; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  [[0, 0], [0.5, Math.round(ymax / 2)], [1, ymax]].forEach(([f, v]) => {
    const ty = y0 + h - 1 - f * (h - 3); ctx.fillRect(x0 + L - 7, ty, 3, 0.6); ctx.fillText(fmtK(v), x0 + L - 9, ty);
  });
  ctx.fillStyle = colorOf(bc);
  const bw = Math.max(0.8, D.binsize / (xmax - XMIN) * (pw - 2));
  for (const [s, n] of D.bins[bc]) {
    if (s >= xmax || s < XMIN) continue;
    const bh = Math.min(1, n / ymax) * (h - 3);
    ctx.fillRect(X(s), y0 + h - 1 - bh, bw, bh);
  }
  // group consensus bands (dashed) -- drawn per panel unless the caller draws them across a stack
  if (g && g.size > 1 && !(opts && opts.noBands)) {
    ctx.save(); ctx.strokeStyle = g.color; ctx.globalAlpha = 0.75; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
    for (const p of g.bands) if (p < xmax) { ctx.beginPath(); ctx.moveTo(X(p) + 0.5, y0 + 1); ctx.lineTo(X(p) + 0.5, y0 + h - 1); ctx.stroke(); }
    ctx.restore();
  }
  // this barcode's detected bands: small triangles at the top
  ctx.fillStyle = '#444';
  for (const [p, wgt] of D.peaks[bc] || []) {
    if (p >= xmax) continue; const px = X(p);
    ctx.globalAlpha = 0.35 + 0.65 * wgt; ctx.beginPath(); ctx.moveTo(px - 3, y0 + 1); ctx.lineTo(px + 3, y0 + 1); ctx.lineTo(px, y0 + 6); ctx.closePath(); ctx.fill();
  }
  ctx.globalAlpha = 1;
  ctx.fillStyle = '#777'; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.font = '8px sans-serif';
  for (let t = Math.ceil(XMIN / xs) * xs; t <= xmax; t += xs) {
    const tx = X(t); ctx.fillRect(tx, y0 + h, 0.6, 3);
    if (tx < x0 + w - 8) ctx.fillText(t >= 1000 ? (t / 1000) + 'k' : t, tx, y0 + h + 3);
  }
}
function renderHist(order) {
  const host = document.getElementById('hist'); host.innerHTML = '';
  const bcs = order.filter(b => D.bins[b]);
  const bySpecies = !!(openCard && openCard.id === 'sec-hist');
  host.className = bySpecies ? 'by-species' : '';
  if (!bcs.length) { host.innerHTML = '<p class="note">No passing barcode has read-length data.</p>'; return; }
  const xmax = histXmax(bcs), xs = niceStep(xmax - XMIN);
  const spGroups = speciesGroups(bcs), lim = yLimits(spGroups, xmax);
  if (bySpecies) {
    // one column per species; inside it the strain groups stacked one under another, in heatmap order
    for (const [sp, members] of spGroups) {
      const gs = []; members.forEach(b => { const g = STR.groupOf[b]; if (!gs.includes(g)) gs.push(g); });
      const GH = 20;
      const H = HC.head + gs.length * GH + members.length * (HC.colPh + HC.colGap) + 6;
      const cv = document.createElement('canvas'); host.appendChild(cv);
      const ctx = ctxFor(cv, HC.colW, H);
      ctx.fillStyle = spColor[sp] || '#9aa0a6'; ctx.fillRect(HC.left, 6, 12, 12);
      ctx.fillStyle = '#222'; ctx.font = 'italic bold 12px sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText(sp, HC.left + 18, 12, HC.colW - HC.left - 80);
      ctx.font = '10px sans-serif'; ctx.fillStyle = '#666'; ctx.textAlign = 'right';
      const ngr = gs.filter(g => g.size > 1).length;
      ctx.fillText(members.length + ' barcode' + (members.length > 1 ? 's' : '') + ' · ' + ngr + ' group' + (ngr === 1 ? '' : 's'), HC.colW, 12);
      let y = HC.head;
      for (const g of gs) {
        const mem = members.filter(b => STR.groupOf[b] === g);
        // group header strip
        ctx.fillStyle = g.size > 1 ? g.color : '#bdbdbd'; ctx.globalAlpha = 0.15; ctx.fillRect(HC.left, y, HC.colW - HC.left, GH - 4); ctx.globalAlpha = 1;
        ctx.fillStyle = g.size > 1 ? g.color : '#888'; ctx.font = 'bold 11px sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        ctx.fillText(g.label + (g.size > 1 ? '  · ' + g.size + ' barcodes · mean ' + (100 * g.mean).toFixed(1) + '%' : '  · unique'), HC.left + 4, y + (GH - 4) / 2);
        y += GH;
        const yTop = y + 14;
        mem.forEach((bc, i) => {
          drawPanel(ctx, bc, 0, y + 14, HC.colW, HC.colPh, xmax, lim.get(bc), xs, fmtK(reads(bc)) + ' reads', { noBands: true });
          y += HC.colPh + HC.colGap;
        });
        if (g.size > 1) {   // shared bands as continuous dashed lines through the whole stack
          const yBot = y - HC.colGap + 14, X = v => HC.left + 1 + (v - XMIN) / (xmax - XMIN) * (HC.colW - HC.left - 2);
          ctx.save(); ctx.strokeStyle = g.color; ctx.globalAlpha = 0.7; ctx.lineWidth = 1; ctx.setLineDash([4, 3]);
          for (const p of g.bands) if (p < xmax) { ctx.beginPath(); ctx.moveTo(X(p) + 0.5, yTop); ctx.lineTo(X(p) + 0.5, yBot); ctx.stroke(); }
          ctx.restore();
        }
      }
    }
    return;
  }
  const W = HC.cols * HC.pw + (HC.cols - 1) * HC.gapX, perCv = HC.cols * HC.rowsPerCanvas;
  for (let start = 0; start < bcs.length; start += perCv) {
    const chunk = bcs.slice(start, start + perCv), rows = Math.ceil(chunk.length / HC.cols);
    const cv = document.createElement('canvas'); host.appendChild(cv);
    const ctx = ctxFor(cv, W, rows * (HC.ph + HC.gapY) + 4);
    chunk.forEach((bc, i) => {
      const x0 = (i % HC.cols) * (HC.pw + HC.gapX), y0 = Math.floor(i / HC.cols) * (HC.ph + HC.gapY) + 16, d = dominant(bc);
      drawPanel(ctx, bc, x0, y0, HC.pw, HC.ph, xmax, lim.get(bc), xs, (d ? shortSp(d.sp) : 'unclassified') + ' · ' + fmtK(reads(bc)));
    });
  }
}
document.getElementById('binW').textContent = D.binsize;
document.querySelectorAll('.minlen').forEach(e => { e.textContent = D.min_len || 0; });
function histOrder() {
  const p = passing().filter(b => D.bins[b]);
  return STR.order.length ? STR.order.concat(p.filter(b => !STR.groupOf[b])) : p;
}

// ---------------------------------------------------------------- downloads
function stamp() { return ((D.meta && D.meta.run) ? D.meta.run.replace(/[^\w.-]+/g, '_') + '_' : '') + 'min' + threshold; }
function save(blob, name) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }
const tsv = rows => new Blob([rows.map(r => r.join('\t')).join('\n') + '\n'], { type: 'text/tab-separated-values' });
const TSV = {
  barcodes() {
    const rows = [['barcode', 'plate', 'well', 'reads', 'status_at_' + threshold, 'mean_read_length_bp', 'sd_read_length_bp', 'dominant_species', 'dominant_pct_of_classified', 'classified_reads', 'second_species', 'strain_group_at_' + simThr]];
    const n = D.n_plates * D.plate_wells;
    for (let i = 1; i <= n; i++) { const bc = nameFor(i), w = well(i), d = dominant(bc), l = (D.lenstats || {})[bc], g = STR.groupOf[bc];
      rows.push([bc, w.plate, w.well, reads(bc), reads(bc) >= threshold ? 'pass' : 'below', l ? l.mean : '', l ? l.sd : '', d ? d.sp : '', d ? d.pct.toFixed(1) : '', d ? d.tot : '', d ? d.second : '', g ? g.label + (g.size > 1 ? '' : ' (unique)') : '']); }
    save(tsv(rows), stamp() + '_barcode_summary.tsv');
  },
  species() {
    const bcs = passing().filter(b => D.species[b]);
    const sps = [...new Set(bcs.flatMap(b => Object.keys(D.species[b])))].sort();
    save(tsv([['species', ...bcs], ...sps.map(sp => [sp, ...bcs.map(b => D.species[b][sp] || 0)])]), stamp() + '_species_by_barcode.tsv');
  },
  speciesAll() {
    const bcs = passing().filter(b => D.species_all[b]);
    const sps = [...new Set(bcs.flatMap(b => Object.keys(D.species_all[b])))].sort();
    save(tsv([['species', ...bcs], ...sps.map(sp => [sp, ...bcs.map(b => D.species_all[b][sp] || 0)])]), stamp() + '_species_by_barcode_unfiltered.tsv');
  },
  matrix() {
    const o = STR.order;
    save(tsv([['barcode', ...o], ...o.map(r => [r, ...o.map(c => simOf(r, c).toFixed(3))])]), stamp() + '_lcp_similarity_' + simMetric + '.tsv');
  },
  groups() {
    const rows = [['group', 'species', 'n_barcodes', 'barcodes', 'wells', 'mean_similarity', 'min_similarity', 'shared_bands_bp', 'similarity_threshold', 'metric']];
    for (const g of STR.groups) rows.push([g.label + (g.size > 1 ? '' : ' (unique)'), g.species, g.size, g.members.join(','),
      g.members.map(b => { const w = well(num(b)); return w.plate + ':' + w.well; }).join(','), g.mean.toFixed(3), g.min.toFixed(3), g.bands.join(','), simThr, simMetric]);
    save(tsv(rows), stamp() + '_strain_groups_' + Math.round(simThr * 1000) / 10 + 'pct.tsv');
  },
  lcp() {
    const rows = [['barcode', 'length_bin_start_bp', 'bin_width_bp', 'reads']];
    for (const b of histOrder()) for (const [s, n] of D.bins[b]) rows.push([b, s, D.binsize, n]);
    save(tsv(rows), stamp() + '_read_length_profiles.tsv');
  },
};
function canvasesToPng(list, name, sideBySide) {
  const gap = sideBySide ? 18 : 0;
  const W = sideBySide ? list.reduce((t, c) => t + c._w, 0) + gap * (list.length - 1) : Math.max(...list.map(c => c._w));
  const H = sideBySide ? Math.max(...list.map(c => c._h)) : list.reduce((t, c) => t + c._h, 0);
  const sc = Math.min(SCALE, 30000 / Math.max(W, H));
  const out = document.createElement('canvas'); out.width = Math.round(W * sc); out.height = Math.round(H * sc);
  const ctx = out.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, out.width, out.height);
  let pos = 0;
  for (const c of list) {
    if (sideBySide) { ctx.drawImage(c, pos * sc, 0, c._w * sc, c._h * sc); pos += c._w + gap; }
    else { ctx.drawImage(c, 0, pos * sc, c._w * sc, c._h * sc); pos += c._h; }
  }
  out.toBlob(b => save(b, name), 'image/png');
}
const PNG = {
  plates: () => canvasesToPng([document.getElementById('plates')], stamp() + '_plates.png'),
  heatmap: () => { const l = [...document.querySelectorAll('#heatmaps canvas')]; if (l.length) canvasesToPng(l, stamp() + '_lcp_heatmaps_' + Math.round(simThr * 1000) / 10 + 'pct.png'); },
  hist: () => { const l = [...document.querySelectorAll('#hist canvas')]; const by = document.getElementById('hist').classList.contains('by-species');
    if (l.length) canvasesToPng(l, stamp() + (by ? '_lcp_histograms_by_species.png' : '_lcp_histograms.png'), by); },
};
document.querySelectorAll('[data-tsv]').forEach(b => b.addEventListener('click', () => TSV[b.dataset.tsv]()));
document.querySelectorAll('[data-png]').forEach(b => b.addEventListener('click', () => PNG[b.dataset.png]()));
document.getElementById('printBtn').addEventListener('click', () => window.print());
// A4: stack the plates vertically so wells stay readable on paper
const printMQ = window.matchMedia('print');
function setPrintLayout(isPrint) { const n = isPrint ? 1 : 2; if (n !== platesPerRow) { platesPerRow = n; renderPlates(); } }
window.addEventListener('beforeprint', () => setPrintLayout(true));
window.addEventListener('afterprint', () => setPrintLayout(false));
if (printMQ.addEventListener) printMQ.addEventListener('change', e => setPrintLayout(e.matches));
if (printMQ.matches) platesPerRow = 1;

// ---------------------------------------------------------------- render
function renderAll() {
  const pass = passing();
  computeStrains(pass);
  renderStats(pass); renderPlates(); renderPlateLegend(); renderTable(pass); renderLegends(pass);
  renderHeatmap(); renderGroups(); renderHist(histOrder()); renderSubs(pass);
  buildMethods();
}
// strain-similarity slider
const simR = document.getElementById('simRange');
function setSim(v) {
  simThr = Math.min(1, Math.max(0.9, v / 100)); simR.value = simThr * 100;
  document.getElementById('simVal').textContent = (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%';
  if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; renderAll(); }); }
}
simR.addEventListener('input', () => setSim(+simR.value));
simR.addEventListener('change', () => saveState());
function renderSubs(pass) {
  const sp = new Set(pass.map(b => (dominant(b) || {}).sp).filter(Boolean));
  const nl = pass.filter(b => D.bins[b]).length;
  const set = (id, t) => { const e = document.getElementById(id); if (e) e.textContent = t; };
  set('sub-plates', pass.length + ' of ' + (D.n_plates * D.plate_wells) + ' wells pass (≥ ' + threshold.toLocaleString() + ' reads)');
  set('sub-species', pass.length + ' barcodes · ' + sp.size + ' dominant species');
  const ng = STR.groups.filter(g => g.size > 1).length, thrTxt = (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%';
  set('sub-heatmap', STR.blocks.length + ' species block' + (STR.blocks.length === 1 ? '' : 's') + ' · ' + ng + ' strain groups at ≥ ' + thrTxt);
  set('sub-groups', ng + ' groups + ' + (STR.groups.length - ng) + ' unique at ≥ ' + thrTxt);
  set('sub-hist', nl + ' read-length profiles, ordered by strain group');
  set('sub-methods', 'wet lab · sequencing · bioinformatics · references');
  if (openCard) document.getElementById('ovSub') && (document.getElementById('ovSub').textContent = (openCard.querySelector('.sub') || {}).textContent || '');
}

// ---------------------------------------------------------------- panels -> full page
const overlay = document.getElementById('overlay'), ovBody = document.getElementById('ovBody');
const ovCtl = document.getElementById('ovCtl'), controlsEl = document.getElementById('controls');
let openCard = null, cardPh = null, ctlPh = null;
function openPanel(card) {
  if (openCard) closePanel();
  openCard = card;
  const body = card.querySelector('.card-body');
  cardPh = document.createComment('panel'); body.before(cardPh); ovBody.insertBefore(body, document.getElementById('ovFoot'));
  ctlPh = document.createComment('controls'); controlsEl.before(ctlPh); ovCtl.appendChild(controlsEl);
  document.getElementById('ovIco').innerHTML = card.querySelector('.ico').innerHTML;
  document.getElementById('ovTitle').innerHTML = card.querySelector('h2').textContent +
    ' <span class="sub" id="ovSub"></span>';
  overlay.hidden = false; document.body.classList.add('ov-open'); ovBody.scrollTop = 0;
  renderSubs(passing());
  if (card.id === 'sec-hist') renderHist(histOrder());
  if (card.id === 'sec-plates') renderPlates();
  if (location.hash !== '#open-' + card.id) history.pushState({ panel: card.id }, '', '#open-' + card.id);
}
function closePanel(fromHistory) {
  if (!openCard) return;
  cardPh.replaceWith(ovBody.querySelector('.card-body')); ctlPh.replaceWith(controlsEl);
  overlay.hidden = true; document.body.classList.remove('ov-open');
  const c = openCard; openCard = null;
  if (c.id === 'sec-hist') renderHist(histOrder());
  if (!fromHistory && location.hash === '#open-' + c.id) history.back();
  c.scrollIntoView({ block: 'nearest' });
}
document.querySelectorAll('.card').forEach(card => card.addEventListener('click', e => {
  if (e.target.closest('input, label, a, select, textarea') || (e.target.closest('button') && !e.target.closest('.expand'))) return;
  if (window.getSelection && String(window.getSelection()).length) return;   // let people select text
  openPanel(card);
}));
document.getElementById('ovClose').addEventListener('click', () => closePanel());
document.getElementById('ovBack2').addEventListener('click', () => closePanel());
document.addEventListener('keydown', e => { if (e.key === 'Escape' && openCard) closePanel(); });
window.addEventListener('popstate', () => {
  const id = location.hash.replace(/^#(open-)?/, ''), card = id && document.getElementById(id);
  if (card && card.classList.contains('card')) { if (openCard !== card) openPanel(card); } else closePanel(true);
});
window.addEventListener('beforeprint', () => closePanel());

// ---------------------------------------------------------------- materials and methods
function buildMethods() {
  const M = (D.meta && D.meta.methods) || {}, P = M.params || {}, V = M.versions || {}, DB = M.db || {};
  const n = v => (v === undefined || v === null || v === '') ? null : (typeof v === 'number' ? v.toLocaleString('en-GB') : v);
  const ver = k => V[k] ? ' v' + V[k] : '';
  const nBc = P.n_reference_barcodes || 192;
  const minQ = n(P.min_qual) || 8, minL = n(P.min_len) || '300', maxL = n(P.max_len) || '6,000';
  const win = n(P.window) || 150, maxD = (P.max_distance || 5);
  const qc = n(P.min_reads) || n(D.default_threshold), stopF = P.stop_fraction ? Math.round(P.stop_fraction * 100) : 90;
  const spPct = n(P.min_species_pct) || 10, bin = n(D.binsize) || 10;
  const minLen = n(P.min_read_length) || D.min_len || 500;
  const dbTxt = DB.description ? esc(DB.description) : 'a GTDB-based Metabuli reference database' + (DB.name ? ' (<code>' + esc(DB.name) + '</code>)' : '');
  const dbDate = DB.created ? ', built ' + esc(DB.created) : '';
  const run = D.meta || {};
  const thisRun = [];
  if (run.run) thisRun.push('run <i>' + esc(run.run) + '</i>');
  if (run.files_read) thisRun.push(run.files_read + ' fastq files');
  const totalReads = allBcs.reduce((t, b) => t + reads(b), 0);
  const seen = allBcs.filter(b => reads(b) > 0).length;
  const R = (k) => '<sup>' + k + '</sup>';
  const html = `
<h3>ON-rep-seq library preparation</h3>
<p>Isolates were processed with the ON-rep-seq protocol${R(1)}, which enriches bacterial genomic DNA for repetitive
extragenic palindromic regions by Rep-PCR${R(2)} and tags each sample with a nanopore-compatible barcode, so that a
single sequencing library carries up to 96 (here up to ${nBc}) isolates. Genomic DNA (about 20&nbsp;ng&nbsp;&micro;l<sup>&minus;1</sup>)
was amplified in Rep-PCR-1 with the (GTG)<sub>5</sub> primer in 25&nbsp;&micro;l reactions containing 5&nbsp;&micro;l 5&times; PCRBIO HiFi
buffer, 0.25&nbsp;&micro;l PCRBIO HiFi Polymerase (PCR Biosystems), 4&nbsp;&micro;l (GTG)<sub>5</sub> primer (5&nbsp;&micro;M) and 1&nbsp;&micro;l DNA
(95&nbsp;&deg;C for 5&nbsp;min; 30 cycles of 95&nbsp;&deg;C for 30&nbsp;s, 45&nbsp;&deg;C for 1&nbsp;min and 62&nbsp;&deg;C for 4&nbsp;min; 72&nbsp;&deg;C for 5&nbsp;min).
Sample-specific barcodes were incorporated in a dual-stage Rep-PCR-2 (25&nbsp;&micro;l: 12&nbsp;&micro;l PCRBIO UltraMix, 2&nbsp;&micro;l of the
corresponding repBC primer (10&nbsp;&micro;M) and 1&nbsp;&micro;l Rep-PCR-1 product). Three cycles of 95&nbsp;&deg;C for 30&nbsp;s, 45&nbsp;&deg;C for 1&nbsp;min
and 62&nbsp;&deg;C for 4&nbsp;min allowed annealing of the (GTG)<sub>5</sub> region, followed by 10 cycles of 95&nbsp;&deg;C for 30&nbsp;s,
65&nbsp;&deg;C for 1&nbsp;min and 72&nbsp;&deg;C for 4&nbsp;min for hybridisation of the full-length primer, and a final extension at
72&nbsp;&deg;C for 5&nbsp;min. Each repBC primer consists (5&prime;&rarr;3&prime;) of a 15-bp spacer, a 24-bp barcode and the
(GTG)<sub>5</sub> pairing region; the spacer moves the barcode away from the low-quality start of the read and was designed
to avoid secondary structures. Barcoded products (10&nbsp;&micro;l per sample; not equimolar, because the amplicon-length
distribution is sample-specific) were pooled and purified with AMPure XP beads (Beckman Coulter; 100&nbsp;&micro;l pool to
50&nbsp;&micro;l beads, i.e. 0.5&times;), which removes primer dimers and shifts the library towards longer fragments. The pool was
quantified with the Qubit dsDNA HS assay (Thermo Fisher Scientific), and about 66&nbsp;ng was used for Oxford Nanopore Technologies
(ONT) ligation-based library preparation (originally the 1D amplicon-by-ligation protocol on R9.4.1 flow cells,
with 80% ethanol in all wash steps)${R(1)}. In the original validation, 48 isolates representing 38 species were all identified
to species level, and strains of <i>Salmonella enterica</i>, <i>Listeria monocytogenes</i> and <i>Bacillus cereus</i> were
discriminated with a sensitivity at least similar to whole-genome sequencing${R(1)}.</p>

<h3>Sequencing and real-time analysis</h3>
<p>Libraries were sequenced on ONT flow cells controlled by MinKNOW with live basecalling, and only reads passing the
MinKNOW quality filter (<code>fastq_pass</code>) were analysed. Data were analysed while sequencing was running with
RT-Rep-Seq${ver('rt_rep_seq')}, a real-time implementation of the ON-rep-seq analysis. RT-Rep-Seq watches the MinKNOW
output folder and runs two decoupled processes: each new fastq file is demultiplexed exactly once, and whenever the
demultiplexed data have grown, the downstream analysis below is repeated on all reads demultiplexed so far. Each cycle is archived
as a separate, time-stamped result set, and the newest is mirrored to a fixed location, so that the identification of
every isolate can be followed as reads accumulate. In single-pass mode (finished runs), demultiplexing stops once
${stopF}% of the barcodes used have reached ${qc} reads.${thisRun.length ? ' For the run reported here (' + thisRun.join(', ') + '), ' + totalReads.toLocaleString('en-GB') + ' reads were assigned to ' + seen + ' barcodes.' : ''}</p>

<h3>Demultiplexing</h3>
<p>Reads were demultiplexed with torchlex${ver('torchlex')} (G. Madjarov; part of the RubyRed pipeline, github.com/reillyeo/RubyRed), a
FASTQ demultiplexer for spacer&ndash;barcode&ndash;linker tagged nanopore amplicons, against the ${nBc} repBC barcodes. Reads with a
mean Phred quality below ${minQ} or a length outside ${minL}&ndash;${maxL}&nbsp;bp were discarded. In a ${win}-bp window at each end of
the read, the 15-bp spacer and the (GTG)<sub>5</sub> linker were located as alignment anchors, the 24-bp sequence between them was
extracted, and it was compared with all reference barcodes by Levenshtein distance (polyleven${ver('polyleven')}). A read was
assigned to the closest barcode if the edit distance was &le;&nbsp;${maxD - 1}; otherwise the reverse complement was tested in the same way.</p>

<h3>Taxonomic classification</h3>
<p>Only reads &ge;&nbsp;${minLen}&nbsp;bp were used for classification. They were tagged with the barcode name and
classified together in a single run of Metabuli${ver('metabuli')}${R(3)} in long-read mode (<code>--seq-mode 3</code>). Metabuli
combines amino-acid and DNA <i>k</i>-mer matching, which improves specificity at species level. The reference was ${dbTxt}${dbDate},
using the Genome Taxonomy Database (GTDB)${R(4)} taxonomy. Each classified read was assigned to the species-rank ancestor of its
Metabuli taxon, and reads resolved only above species rank were not counted. The dominant species of a barcode is the most
abundant species after species with fewer than ${spPct}% of the barcode's species-level reads had been removed as noise (barcode
cross-talk or misclassification), because every barcode is expected to hold a single isolate. Its share is reported against
all species-level reads of the barcode. Because each read (a single Rep-PCR amplicon) is classified independently, reads from
conserved genes, mobile elements or plasmids shared between related taxa may be assigned to neighbouring species, and very
closely related species can split each other's reads. Minor species should therefore be read as indications, not proof, of
contamination or mixed cultures.</p>

<h3>Read-length statistics</h3>
<p>Read length is reported per barcode and for all passing barcodes as mean &plusmn; standard deviation of all demultiplexed
reads (pooled over barcodes for the run summary).</p>

<h3>Read-length count profiles and strain grouping</h3>
<p>For each barcode, reads &ge;&nbsp;${minLen}&nbsp;bp were counted per read length to give the read-length count profile (LCP),
the sequencing equivalent of a Rep-PCR gel fingerprint${R(1)}. As in gel-based Rep-PCR typing, strain identity was judged from
the presence or absence of bands (LCP peaks) rather than from their intensities, which vary with PCR competition and
sequencing depth. Each LCP was re-binned to a ${PP.grid_bp}-bp grid and smoothed with a Gaussian kernel (&sigma;&nbsp;=&nbsp;${PP.sigma_bp}&nbsp;bp),
and local maxima were called as bands when their prominence was &ge;&nbsp;${Math.round(PP.min_rel_prom * 100)}% of the strongest band
and &ge;&nbsp;${PP.min_abs_prom} reads per bin (SciPy${ver('scipy')}${R(5)}). Each band was weighted by its prominence relative to the
strongest band, saturating at 1 for bands &ge;&nbsp;${Math.round(PP.w_full * 100)}% of it, so that faint bands near the detection limit
contribute less, and intensity has no influence above that level. Bands of two barcodes were paired (closest first) when
their positions differed by &le;&nbsp;max(${PP.tol_bp}&nbsp;bp, ${(PP.tol_frac * 100).toFixed(1)}% of the band length). Similarity was the
weighted Dice coefficient${R(9)}: the summed weight of paired bands divided by the summed weight of all bands of both barcodes.
Barcodes were clustered by UPGMA (average linkage) on 1&nbsp;&minus;&nbsp;similarity${R(10)}, separately within each dominant species, and
<b>strain groups</b> were defined by cutting the tree at the similarity threshold
(${(simThr * 100).toFixed(1).replace(/\.0$/, '')}% in this report; adjustable from 90 to 100%). With these settings, eight technical
replicates (separate Rep-PCR wells of the same DNA) formed a single group at &ge;&nbsp;90%, whereas barcodes of different species did
not exceed 0.77. Bands shared by a strain group are those present in at least half of its members. As an alternative, the report
can use the Pearson correlation of the relative exact-length distributions. Barcodes that share a dominant species and fall into
one strain group are interpreted as the same strain (or technical replicates), and barcodes of the same species in different
groups as different strains${R(1)}.</p>

<h3>Quality threshold and reporting</h3>
<p>Barcodes with fewer than ${qc} demultiplexed reads were treated as insufficiently sequenced (quality-control threshold).
All barcodes with reads are classified and profiled, so the threshold can be changed in the interactive report. The
classification table, strain groups, similarity heatmaps and LCP histograms are then recomputed for the barcodes that pass,
and the two 96-well plate maps show which positions pass or which strain group each well belongs to. Tables are exported as
tab-separated files and figures as PNG. Analyses used Python${ver('python')} with NumPy${ver('numpy')}${R(6)}, SciPy${ver('scipy')}${R(5)},
Matplotlib${ver('matplotlib')}${R(7)} and Biopython${ver('biopython')}${R(8)}.</p>

<h3>Code availability</h3>
<p>RT-Rep-Seq: ${M.repo ? '<a href="' + esc(M.repo) + '">' + esc(M.repo) + '</a>' : '<i>repository URL to be added on release</i>'}.
Original ON-rep-seq pipeline: <a href="https://github.com/lauramilena3/On-rep-seq">https://github.com/lauramilena3/On-rep-seq</a>
(documentation <a href="https://on-rep-seq.readthedocs.io">https://on-rep-seq.readthedocs.io</a>; version described in ref. 1:
<a href="https://doi.org/10.5281/zenodo.3384842">doi:10.5281/zenodo.3384842</a>).
Metabuli: <a href="https://github.com/steineggerlab/Metabuli">https://github.com/steineggerlab/Metabuli</a>.
GTDB: <a href="https://gtdb.ecogenomic.org">https://gtdb.ecogenomic.org</a>.
polyleven: <a href="https://github.com/fujimotos/polyleven">https://github.com/fujimotos/polyleven</a>.</p>

<h3>References</h3>
<ol class="refs">
<li>Krych, &#321;., Castro-Mej&iacute;a, J. L., Forero-Junco, L. M., Moesby, D. N., Mikkelsen, M. B., Rasmussen, M. A., Sykulski, M. &amp; Nielsen, D. S.
DNA enrichment and tagmentation method for species-level identification and strain-level differentiation using ON-rep-seq.
<i>Commun. Biol.</i> <b>2</b>, 369 (2019). <a href="https://doi.org/10.1038/s42003-019-0617-x">doi:10.1038/s42003-019-0617-x</a></li>
<li>Versalovic, J., Schneider, M., de Bruijn, F. J. &amp; Lupski, J. R. Genomic fingerprinting of bacteria using repetitive
sequence-based polymerase chain reaction. <i>Methods Mol. Cell. Biol.</i> <b>5</b>, 25&ndash;40 (1994).</li>
<li>Kim, J. &amp; Steinegger, M. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid
and DNA. <i>Nat. Methods</i> <b>21</b>, 971&ndash;973 (2024). <a href="https://doi.org/10.1038/s41592-024-02273-y">doi:10.1038/s41592-024-02273-y</a></li>
<li>Parks, D. H., Chuvochina, M., Rinke, C., Mussig, A. J., Chaumeil, P.-A. &amp; Hugenholtz, P. GTDB: an ongoing census of bacterial
and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy.
<i>Nucleic Acids Res.</i> <b>50</b>, D785&ndash;D794 (2022). <a href="https://doi.org/10.1093/nar/gkab776">doi:10.1093/nar/gkab776</a></li>
<li>Virtanen, P. <i>et al.</i> SciPy 1.0: fundamental algorithms for scientific computing in Python. <i>Nat. Methods</i>
<b>17</b>, 261&ndash;272 (2020). <a href="https://doi.org/10.1038/s41592-019-0686-2">doi:10.1038/s41592-019-0686-2</a></li>
<li>Harris, C. R. <i>et al.</i> Array programming with NumPy. <i>Nature</i> <b>585</b>, 357&ndash;362 (2020).
<a href="https://doi.org/10.1038/s41586-020-2649-2">doi:10.1038/s41586-020-2649-2</a></li>
<li>Hunter, J. D. Matplotlib: a 2D graphics environment. <i>Comput. Sci. Eng.</i> <b>9</b>, 90&ndash;95 (2007).
<a href="https://doi.org/10.1109/MCSE.2007.55">doi:10.1109/MCSE.2007.55</a></li>
<li>Cock, P. J. A. <i>et al.</i> Biopython: freely available Python tools for computational molecular biology and
bioinformatics. <i>Bioinformatics</i> <b>25</b>, 1422&ndash;1423 (2009). <a href="https://doi.org/10.1093/bioinformatics/btp163">doi:10.1093/bioinformatics/btp163</a></li>
<li>Dice, L. R. Measures of the amount of ecologic association between species. <i>Ecology</i> <b>26</b>, 297&ndash;302 (1945).
<a href="https://doi.org/10.2307/1932409">doi:10.2307/1932409</a></li>
<li>Sokal, R. R. &amp; Michener, C. D. A statistical method for evaluating systematic relationships. <i>Univ. Kansas Sci. Bull.</i>
<b>38</b>, 1409&ndash;1438 (1958).</li>
</ol>`;
  document.getElementById('methodsText').innerHTML = html;
}
function methodsPlain() {
  const el = document.getElementById('methodsText');
  const out = [];
  el.childNodes.forEach(nd => {
    if (nd.nodeType !== 1) return;
    if (nd.tagName === 'H3') out.push('\n' + nd.textContent.trim().toUpperCase());
    else if (nd.tagName === 'OL') [...nd.children].forEach((li, i) => out.push((i + 1) + '. ' + li.textContent.replace(/\s+/g, ' ').trim()));
    else out.push(nd.textContent.replace(/\s+/g, ' ').trim());
  });
  return 'Materials and methods\n' + out.join('\n') + '\n';
}
document.getElementById('copyMethods').addEventListener('click', e => {
  const t = methodsPlain(), b = e.currentTarget;
  const done = () => { b.textContent = 'Copied ✓'; setTimeout(() => { b.textContent = 'Copy text'; }, 1500); };
  if (navigator.clipboard) navigator.clipboard.writeText(t).then(done, () => { save(new Blob([t], { type: 'text/plain' }), 'materials_and_methods.txt'); });
  else save(new Blob([t], { type: 'text/plain' }), 'materials_and_methods.txt');
});
document.getElementById('dlMethods').addEventListener('click', () => save(new Blob([methodsPlain()], { type: 'text/plain' }), stamp() + '_materials_and_methods.txt'));
buildMethods();
// ---------------------------------------------------------------- live auto-update
// The pipeline writes report_stamp.js next to report.html (window.RTREPSEQ_STAMP = "..."). While a run is
// live, the page re-loads that tiny script every 30 s (works from file://, no server) and reloads itself
// when the stamp changes -- keeping threshold, open panel (URL hash) and scroll position.
const STAMP = (D.meta && D.meta.stamp) || '';
const store = (() => { try { const k = '__t'; sessionStorage.setItem(k, '1'); sessionStorage.removeItem(k); return sessionStorage; } catch (e) { return null; } })();
const SKEY = 'rtrepseq:' + location.pathname;
function saveState(extra) {
  if (!store) return;
  try {
    const st = Object.assign(JSON.parse(store.getItem(SKEY) || '{}'), {
      thr: threshold, sim: simThr, metric: simMetric, bySp: bySpeciesHM, plate: plateMode,
      live: document.getElementById('liveChk').checked, sy: openCard ? ovBody.scrollTop : window.scrollY }, extra || {});
    store.setItem(SKEY, JSON.stringify(st));
  } catch (e) {}
}
let restored = {};
try { restored = store ? JSON.parse(store.getItem(SKEY) || '{}') : {}; } catch (e) { restored = {}; }
if (restored.thr !== undefined) threshold = restored.thr;
if (restored.sim) { simThr = restored.sim; simR.value = simThr * 100; }
document.getElementById('simVal').textContent = (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%';
if (restored.metric) { simMetric = restored.metric; document.getElementById('simMetric').value = simMetric; }
if (restored.bySp === false) { bySpeciesHM = false; document.getElementById('bySpecies').checked = false; }
if (restored.plate) { plateMode = restored.plate; document.querySelectorAll('#plateMode button').forEach(x => x.classList.toggle('on', x.dataset.m === plateMode)); }
if (restored.live === false) document.getElementById('liveChk').checked = false;
rng.addEventListener('change', () => saveState()); box.addEventListener('change', () => saveState());
document.getElementById('liveChk').addEventListener('change', () => saveState());
let lastCheck = null;
function pollStamp() {
  const sc = document.createElement('script');
  sc.src = 'report_stamp.js?_=' + Date.now();
  sc.onload = () => {
    sc.remove();
    const fresh = window.RTREPSEQ_STAMP;
    document.getElementById('liveBox').hidden = false;
    lastCheck = new Date();
    document.getElementById('liveNote').textContent = '\u00b7 checked ' + lastCheck.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit' });
    if (fresh && STAMP && fresh !== STAMP && document.getElementById('liveChk').checked && restored.reloadedFor !== fresh) {
      saveState({ reloadedFor: fresh });   // never reload twice for the same stamp (no loops if a copy is half-written)
      location.reload();
    }
  };
  sc.onerror = () => { sc.remove(); };    // no stamp file (static copy) or being replaced: try again later
  document.head.appendChild(sc);
}
if (location.protocol === 'file:' || location.protocol.startsWith('http')) { setTimeout(pollStamp, 3000); setInterval(pollStamp, 30000); }
window.addEventListener('beforeunload', () => saveState());

// settings from the address, e.g. report.html?thr=1000&sim=95&plate=strain#open-sec-plates (shareable views)
try {
  const q = new URLSearchParams(location.search);
  if (q.has('thr')) threshold = Math.max(0, +q.get('thr') || 0);
  if (q.has('sim')) { simThr = Math.min(1, Math.max(0.9, (+q.get('sim') || 90) / 100)); simR.value = simThr * 100;
    document.getElementById('simVal').textContent = (simThr * 100).toFixed(1).replace(/\.0$/, '') + '%'; }
  if (q.get('plate') === 'strain' || q.get('plate') === 'reads') { plateMode = q.get('plate');
    document.querySelectorAll('#plateMode button').forEach(x => x.classList.toggle('on', x.dataset.m === plateMode)); }
} catch (e) {}
setThreshold(threshold);
if (location.hash.length > 1) { const c = document.getElementById(location.hash.replace(/^#(open-)?/, '')); if (c && c.classList.contains('card')) openPanel(c); }
if (restored.sy) requestAnimationFrame(() => requestAnimationFrame(() => { if (openCard) ovBody.scrollTop = restored.sy; else window.scrollTo(0, restored.sy); }));
</script>
</body></html>
"""
