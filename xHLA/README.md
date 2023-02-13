# Requirements

## Download reference

This pipeline require you to download and index a reference genome **without alternative sequences**.

```bash
mkdir Reference
cd Reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
# Remove alternative sequences
samtools faidx -o hg38-noalt.fa hg38.fa `grep "^>" hg38.fa | grep -v "alt" | cut -c 2-`
# Index reference
bwa index hg38-noalt.fa
```

## Run the pipeline

```
nextflow main.nf --reference Reference/hg38-noalt.fa
```
