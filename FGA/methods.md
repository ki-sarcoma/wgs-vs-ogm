# Fraction of Genome Altered (FGA)

## Calculation

Let:

- $L_i$ = length (in base pairs) of the i-th altered genomic segment  
- $G$ = total size of the reference genome considered

Then:

$FGA = \frac{\sum_{i=1}^{n} L_i}{G}$

The total length of the reference genome was determined using chromosome size data from the UCSC Genome Browser. Only autosomal chromosomes (chr 1–22) were considered; sex chromosomes (chrX and chrY) and unplaced contigs were excluded.

- **WGS (hg19):** [hg19 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes), giving $G = 2,881,033,286$ bp
- **OGM (hg38):** [hg38 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes), giving $G = 2,875,001,522$ bp
