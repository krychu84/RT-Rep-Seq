#!/usr/bin/env python3
"""Simulate a small ON-rep-seq run for testing RT-Rep-Seq (this is how test/data/ was made).

Rep-PCR amplicons are cut from three public complete genomes, wrapped in the ON-rep-seq read
structure (spacer - repBC barcode - (GTG)5 ... (GTG)5 - barcode - spacer) and given ONT-like
errors and quality strings. Each simulated strain has its own set of amplicons ("bands"), so the
expected species and strain groups are known:

  repBC01-03  Lactiplantibacillus plantarum  strain A  (3 replicates)  -> one strain group
  repBC04-05  Lactiplantibacillus plantarum  strain B  (2 replicates)  -> second group
  repBC06-08  Escherichia coli               strain A  (3 replicates)  -> one strain group
  repBC09     Escherichia coli               strain B                   -> unique profile
  repBC10-11  Lactococcus lactis             strain A  (2 replicates)  -> one strain group
  repBC12     Lactococcus lactis             strain A, only ~150 reads -> fails a 500-read threshold

Genomes (NCBI RefSeq, downloaded with --download):
  NC_004567.2  L. plantarum WCFS1 · NC_000913.3  E. coli K-12 MG1655 · NC_002662.1  L. lactis IL1403

usage: python test/simulate_test_data.py --download --out test/data
"""
import argparse
import gzip
import os
import random
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
GENOMES = {'plantarum': 'NC_004567.2', 'coli': 'NC_000913.3', 'lactis': 'NC_002662.1'}
COMP = str.maketrans('ACGTN', 'TGCAN')

# barcode -> (genome, strain, number of reads)
DESIGN = {1: ('plantarum', 'A', 1300), 2: ('plantarum', 'A', 1100), 3: ('plantarum', 'A', 1200),
          4: ('plantarum', 'B', 1250), 5: ('plantarum', 'B', 1050),
          6: ('coli', 'A', 1200), 7: ('coli', 'A', 1350), 8: ('coli', 'A', 1000),
          9: ('coli', 'B', 1150),
          10: ('lactis', 'A', 1250), 11: ('lactis', 'A', 1100),
          12: ('lactis', 'A', 150)}


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    with open(path) as f:
        return ''.join(line.strip() for line in f if not line.startswith('>')).upper()


def load_barcodes(path):
    """repBCNN_1st / _2nd entries of 3 lines: spacer, barcode, linker."""
    entries, name, lines = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            if name:
                entries[name] = lines
            name, lines = line[1:], []
        elif line:
            lines.append(line)
    entries[name] = lines
    return entries


def bands_for(genome, rng, n):
    """n amplicons (start, length, weight) spread over the genome, 400-3200 bp."""
    out = []
    for _ in range(n):
        length = int(rng.choice([rng.uniform(450, 1200), rng.uniform(1200, 2200), rng.uniform(2200, 3200)]))
        start = rng.randrange(0, len(genome) - length)
        out.append((start, length, rng.uniform(0.2, 1.0)))
    return out


def mutate(seq, rng, err):
    """ONT-like errors: substitutions, insertions, deletions in roughly 2:1:1."""
    out = []
    for b in seq:
        r = rng.random()
        if r < err / 2:
            out.append(rng.choice('ACGT'))
        elif r < err * 3 / 4:
            out.append(b + rng.choice('ACGT'))
        elif r < err:
            continue
        else:
            out.append(b)
    return ''.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=os.path.join(HERE, 'data'))
    ap.add_argument('--genomes', default=os.path.join(HERE, 'genomes'))
    ap.add_argument('--download', action='store_true', help='fetch the three genomes from NCBI')
    ap.add_argument('--files', type=int, default=6, help='number of fastq.gz chunks (like MinKNOW)')
    ap.add_argument('--seed', type=int, default=2019)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    os.makedirs(args.genomes, exist_ok=True)
    genomes = {}
    for key, acc in GENOMES.items():
        path = os.path.join(args.genomes, acc + '.fa')
        if not os.path.exists(path):
            if not args.download:
                raise SystemExit(f'{path} missing -- run with --download')
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text'
            urllib.request.urlretrieve(url, path)
        genomes[key] = read_fasta(path)

    # strain A: 12 bands; strain B: strain A with 4 bands replaced (clearly a different strain)
    strains = {}
    for key, g in genomes.items():
        a = bands_for(g, rng, 12)
        b = a[:]
        for i in rng.sample(range(12), 4):
            b[i] = bands_for(g, rng, 1)[0]
        strains[(key, 'A')], strains[(key, 'B')] = a, b

    bc = load_barcodes(os.path.join(REPO, 'assets', 'onrep_barcodes.fasta'))
    reads = []
    for n, (key, strain, n_reads) in DESIGN.items():
        sp1, bc1, lk1 = bc[f'repBC{n:02d}_1st']
        sp2, bc2, lk2 = bc[f'repBC{n:02d}_2nd']
        head = sp1 + bc1 + lk1
        tail = rc(lk2) + bc2 + rc(sp2)
        # replicate noise: band weights vary +-30%, as between Rep-PCR wells
        bands = [(s, ln, w * rng.uniform(0.7, 1.3)) for s, ln, w in strains[(key, strain)]]
        weights = [w for _, _, w in bands]
        g = genomes[key]
        for i in range(n_reads):
            if rng.random() < 0.06:            # unspecific background ("smear")
                ln = int(rng.uniform(400, 3000))
                s = rng.randrange(0, len(g) - ln)
            else:
                s, ln, _ = rng.choices(bands, weights)[0]
            amp = g[s:s + ln]
            if rng.random() < 0.5:
                amp = rc(amp)
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(5, 25)))
            seq += mutate(head + amp + tail, rng, rng.uniform(0.03, 0.07))
            if rng.random() < 0.5:
                seq = rc(seq)
            q = ''.join(chr(33 + max(5, min(40, int(rng.gauss(17, 5))))) for _ in seq)
            reads.append((f'sim_{n:02d}_{i:05d} barcode=unclassified', seq, q))

    rng.shuffle(reads)
    os.makedirs(args.out, exist_ok=True)
    per = -(-len(reads) // args.files)
    for k in range(args.files):
        path = os.path.join(args.out, f'SIM00001_pass_rtrepseqtest_{k}.fastq.gz')
        with gzip.open(path, 'wt', compresslevel=9) as f:
            for name, seq, q in reads[k * per:(k + 1) * per]:
                f.write(f'@{name}\n{seq}\n+\n{q}\n')
    print(f'{len(reads)} reads in {args.files} files -> {args.out}')


if __name__ == '__main__':
    main()
