#!/usr/bin/env python3
"""lcp_peaks.py -- peak calling on ON-rep-seq read-length count profiles (LCPs) and
peak-presence similarity between barcodes, for strain-level grouping.

ON-rep-seq fingerprints are read like Rep-PCR gel bands: what matters for "same strain or
not" is WHICH amplicon lengths are present, not exactly how strong each one is (PCR
competition and sequencing depth shift intensities between replicates). So:

  1. peaks   -- each LCP (reads >= min_len) is put on a 5-bp grid, Gaussian-smoothed, and
                local maxima are kept if their prominence is >= `min_rel_prom` of the
                profile's highest peak and >= `min_abs_prom` reads/bin (noise floor).
  2. weight  -- every peak gets a weight w = min(1, relative prominence / `w_full`): clear
                bands (>= 20% of the top band by default) count fully, faint bands near the
                detection limit count proportionally less -- so a faint band that drops
                below detection in a shallower replicate costs little similarity, while a
                clear band present in one profile and absent in the other costs a lot.
                Above `w_full`, intensity has no influence at all.
  3. similarity -- weighted Dice coefficient: peaks of two profiles are paired (greedy,
                closest first) when their positions differ by <= max(tol_bp, tol_frac x
                length); S = sum(w of paired peaks, both sides) / sum(w of all peaks).
                S = 1: identical band patterns; every unpaired clear band lowers S.

The same functions (and parameters) are mirrored in the interactive report's JavaScript.
"""
import numpy as np

# Tuned on technical replicates (8 wells of one DNA, an in-house ON-rep-seq run, repBC89-96):
# replicates stay together at >= 90% (UPGMA), different species never exceed ~0.77.
PARAMS = {
    'grid_bp': 5,          # LCP is re-binned to this grid before smoothing
    'sigma_bp': 10,        # Gaussian smoothing (bp)
    'min_rel_prom': 0.04,  # peak prominence >= 4% of the profile's highest peak ...
    'min_abs_prom': 2.0,   # ... and >= 2 reads per 5-bp bin after smoothing
    'w_full': 0.50,        # peaks with relative prominence >= 50% get full weight 1
    'tol_bp': 20,          # two peaks match if |a-b| <= max(tol_bp, tol_frac * length)
    'tol_frac': 0.02,
}


def _smooth(y, sigma_bins):
    r = int(4 * sigma_bins)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma_bins) ** 2)
    k /= k.sum()
    return np.convolve(y, k, mode='same')


def call_peaks(profile, min_len=500, max_len=6000, params=None):
    """profile: {read_length: count}. Returns list of (position_bp, height, rel_prominence,
    weight), sorted by position."""
    from scipy.signal import find_peaks
    p = dict(PARAMS, **(params or {}))
    g = p['grid_bp']
    n = (max_len - min_len) // g + 1
    y = np.zeros(n)
    for L, c in profile.items():
        if min_len <= L <= max_len:
            y[(L - min_len) // g] += c
    if y.sum() == 0:
        return []
    ys = _smooth(y, p['sigma_bp'] / g)
    top = ys.max()
    idx, props = find_peaks(ys, prominence=max(p['min_abs_prom'], p['min_rel_prom'] * top))
    out = []
    prom_top = props['prominences'].max() if len(idx) else 1
    for i, pr in zip(idx, props['prominences']):
        rel = pr / prom_top
        out.append((int(min_len + i * g + g // 2), round(float(ys[i]), 2), round(float(rel), 4),
                    round(float(min(1.0, rel / p['w_full'])), 4)))
    return out


def tolerance(pos, params=None):
    p = dict(PARAMS, **(params or {}))
    return max(p['tol_bp'], p['tol_frac'] * pos)


def peak_similarity(pa, pb, params=None):
    """Weighted Dice similarity of two peak lists (see module docstring)."""
    if not pa and not pb:
        return 1.0
    if not pa or not pb:
        return 0.0
    pairs = []
    for i, a in enumerate(pa):
        for j, b in enumerate(pb):
            d = abs(a[0] - b[0])
            if d <= tolerance((a[0] + b[0]) / 2, params):
                pairs.append((d, i, j))
    pairs.sort()
    ua, ub, matched = set(), set(), 0.0
    for d, i, j in pairs:
        if i in ua or j in ub:
            continue
        ua.add(i)
        ub.add(j)
        matched += pa[i][3] + pb[j][3]
    total = sum(x[3] for x in pa) + sum(x[3] for x in pb)
    return matched / total if total else 0.0


def similarity_matrix(peaks, barcodes, params=None):
    n = len(barcodes)
    m = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            m[i, j] = m[j, i] = peak_similarity(peaks[barcodes[i]], peaks[barcodes[j]], params)
    return m


def upgma_order(sim, labels):
    """Leaf order of an average-linkage (UPGMA) tree on distance 1 - sim."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform
    if len(labels) < 3:
        return list(labels)
    d = np.clip(1 - sim, 0, None)
    np.fill_diagonal(d, 0)
    return [labels[i] for i in leaves_list(linkage(squareform((d + d.T) / 2, checks=False), 'average'))]
