# Test data

`data/` holds a small **simulated** ON-rep-seq run (6 × fastq.gz, 13,100 reads, ~26 MB) to check an
installation. Reads are Rep-PCR-like amplicons cut from three public complete genomes, in the real
ON-rep-seq read structure (spacer – repBC barcode – (GTG)<sub>5</sub> … (GTG)<sub>5</sub> – barcode – spacer),
with ONT-like errors (3–7%) and quality strings. `simulate_test_data.py` recreates it
(`python test/simulate_test_data.py --download`).

| Genome (NCBI RefSeq) | |
|---|---|
| NC_004567.2 | *Lactiplantibacillus plantarum* WCFS1 |
| NC_000913.3 | *Escherichia coli* K-12 MG1655 |
| NC_002662.1 | *Lactococcus lactis* subsp. *lactis* IL1403 |

## Run the test

```bash
test/run_test.sh            # live: files arrive one by one, like during sequencing
test/run_test.sh --all-now  # all files at once
```

It creates a MinKNOW-style folder in `~/RT-Rep-Seq_test_run/` and starts `RT-Rep-Seq --min-reads 500`
there; the report appears in `~/Desktop/RT-Rep-Seq_SIM_test_<date_time>/latest/report.html`.
Most of the 5–15 min is Metabuli loading its database. Stop with `RT-Rep-Seq --stop`, clean up with
`rm -rf ~/RT-Rep-Seq_test_run ~/Desktop/RT-Rep-Seq_SIM_test_*`.

## Expected result (gtdb226 database, 90% strain similarity)

| Barcodes | Plate 1 | Reads | Species | Strain group |
|---|---|---|---|---|
| repBC01–03 | A1–C1 | ~1,100–1,300 each | *Lactiplantibacillus plantarum* | one group (3 replicates of strain A) |
| repBC04–05 | D1–E1 | ~1,050–1,250 | *Lactiplantibacillus plantarum* | second group (strain B) |
| repBC06–08 | F1–H1 | ~1,000–1,350 | *Escherichia coli* | one group (strain A) |
| repBC09 | A2 | ~1,150 | *Escherichia coli* | unique (strain B) |
| repBC10–11 | B2–C2 | ~1,100–1,250 | *Lactococcus lactis* | one group |
| repBC12 | D2 | 150 | *Lactococcus lactis* | below the 500-read threshold (red well) |

≥ 99% of reads are demultiplexed. *E. coli* reads are partly assigned to the GTDB sister clusters
*Escherichia coli_F* / *Escherichia* sp. — expected for single-read classification; the dominant
species is *E. coli*. Replicates of a simulated strain have identical band sets, so they group at
~100%; real technical replicates are typically 90–98% similar.
