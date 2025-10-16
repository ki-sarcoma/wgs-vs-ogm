# Fraction of Genome Altered (FGA)

## Calculation

Let:

- $L_i$ = length (in base pairs) of the i-th altered genomic segment  
- $G$ = total size of the reference genome considered

Then:

$FGA = \frac{\sum_{i=1}^{n} L_i}{G}$

Only **autosomal chromosomes (chr1–chr22)** were considered; sex chromosomes (chrX and chrY) and unplaced contigs were excluded. The total length of the autosomes was determined using chromosome size data from the UCSC Genome Browser:

- **WGS (hg19):** [hg19 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes), giving $G = 2,881,033,286$ bp
- **OGM (hg38):** [hg38 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes), giving $G = 2,875,001,522$ bp

## Variant Filtering

## Optical Genome Mapping (OGM) Data

## Whole Genome Sequencing (WGS) Data
