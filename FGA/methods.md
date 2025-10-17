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

## Optical Genome Mapping (OGM) Data

## Variant Filtering

Discard chr > 22

Confidence > 0.5  (for example, for gains/losses). A lot of values have confidence 0 though
Mask_overlap_fract < 0.9 (avoid heavily masked regions)

Size > 5,000 bp

Confirm CNVs:
Loss: CopyNumber < 2
Gain: CopyNumber > 2
log2 ratio > +0.3 (gain) or < –0.3 (loss)
Adjust each CNV segment for tumor fraction

log2_ratio = log2(fractionalCopyNumber / 2)

num_overlap_DGV_calls < 5 (If many overlapping calls in DGV, likely benign)

## Whole Genome Sequencing (WGS) Data
