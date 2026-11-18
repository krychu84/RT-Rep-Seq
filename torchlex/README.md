# torchlex (bundled)

FASTQ demultiplexer for Oxford Nanopore amplicons tagged with a spacer–barcode–linker structure.
torchlex is part of the **RubyRed** pipeline (https://github.com/reillyeo/RubyRed; method note:
https://nanoporetech.com/resource-centre/torchlex-method-real-time-demultiplexing-barcoded-oxford-nanopore-reads).

Contact: **Gjorgji Madjarov** — gjorgji.madjarov@finki.ukim.mk

This copy is adapted to the ON-rep-seq repBC barcode set (15-bp spacer, 24-bp barcode, (GTG)5
linker) and is imported directly by RT-Rep-Seq (no separate installation needed). Requires
biopython, polyleven and numpy. It keeps its authors' terms; the MIT licence of RT-Rep-Seq covers
the rest of this repository.

Algorithm: reads are filtered by mean Phred quality and length; in a window at each read end the
spacer and linker anchors are located, the 24-bp barcode between them is extracted and compared
with all reference barcodes by Levenshtein distance; the closest barcode is accepted if the
distance is below `max_barcode_distance` (default 5), otherwise the reverse complement is tried.
