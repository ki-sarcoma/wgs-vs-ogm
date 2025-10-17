# Fraction of Genome Altered (FGA)

## Calculation

Let:

- $L_i$ = length (in base pairs) of the i-th altered genomic segment  
- $G$ = total size of the reference genome considered

Then:

$FGA = \frac{\sum_{i=1}^{n} L_i}{G}$

Onlyautosomal chromosomes (chr 1–22) were considered; sex chromosomes (chrX and chrY) and unplaced contigs were excluded. The total length of the autosomes was determined using chromosome size data from the UCSC Genome Browser:

- **WGS (hg19):** [hg19 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes), giving $G = 2,881,033,286$ bp
- **OGM (hg38):** [hg38 chrom.sizes](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes), giving $G = 2,875,001,522$ bp

## Optical Genome Mapping (OGM) Data

## Variant Filtering

### 1. Chromosomal Scope

- Include only autosomal chromosomes (chr 1–22).

### 2. Variant Size

- Include only variants with a size ≥ 50 kbp.

### 3. CNV Confirmation

- Include variants based on log<sub>2</sub> ratio:
  - **Gain:** log<sub>2</sub> ratio > +0.3
  - **Loss:** log<sub>2</sub> ratio < –0.3

The log<sub>2</sub> ratio is calculated using the fractional copy number:

$\log_2 \text{ ratio} = \log_2 \frac{fractionalCopyNumber}{2}$

> [!NOTE] Todo
> Adjust each CNV segment for tumor fraction.

### 4. Confidence Threshold

> [!NOTE] Todo  
>
> - Filter out variants with low confidence scores (`Confidence` column).
> - I've observed quite a lot of zeros, so this value may not be entirely reliable.

### 5. Masked Regions

> [!NOTE] Todo
>
> - Exclude variants overlapping heavily masked regions (`Mask_overlap_fract` column).

### 6. Overlap with DGV

> [!NOTE] Todo
>
> - Exclude CNVs with high overlap with the Database of Genomic Variants (DGV), as these are often benign (`num_overlap_DGV_calls` column).
> - We should include them only if the goal is to capture all genomic alterations.

## Whole Genome Sequencing (WGS) Data
