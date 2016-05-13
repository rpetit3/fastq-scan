# fastq-stats

I needed a quick method to output summary statistics of an input FASTQ. I found
the question: "How To Efficiently Parse A Huge Fastq File?" on Biostars (link
below). From that question I used code from Pierre Lindenbaum's C++ solution as
the base for this program.

Biostars Link: https://www.biostars.org/p/10353/#10358


zcat YOUR_FASTQ.gz | fastq_stats GENOME_SIZE
