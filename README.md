# Understanding biological functions from a list of genes

Many curiosities or analysis boils down to a potential list of candidates. For
example, in cases of genome sequence or genome-wide assays, a set of analysis
leads to a list of genes to investigate. A natural next step is to first look
at in what set of biological functions they are involved. To do this, all we
need to do is to feed that list to existing R or Python packages and get a
visualization.

In this tiny project, we are going to simply do that. We will get a list of
genes after searching for a sequence of interest in promoter region of all the
genes and then look at their biological functions.

You will need the following datasets to do the project.

- [Human genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/), search for `hg38.fa.gz` and download

- Dreg from [EMBOSS](https://anaconda.org/bioconda/emboss)

- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) Bioconductor Package

- [bedtools](https://anaconda.org/channels/bioconda/packages/bedtools/overview) 

## Environment setup

If you have a Windows system, install Windows Subsystem for Linux, WSL. It will give you direct access to Linux system. Students with either a Linux or MacOS system will have easy time setting up the environment to execute the project. 

Setup [Miniforge](https://github.com/conda-forge/miniforge), and create a virtual environment with the following command:

```
mamba create -n go_enrichment python=3.12 
```

Once the above run successfully, you can do the following. 

```
mamba activate go_enrichment
```

```
mamba install bioconda::emboss
```

## Target sequences to search in the promoter region

As explained in the class, there about 1700 Transcription Factors in humans categorized in different families. You can know more about them [here](https://www.sciencedirect.com/science/article/pii/S0092867418301065) if you are interested. For the assignment, we will have one example sequence "GCGC..GCGC" (the dot represents any nucleotide). It is a recognition sequence for NRF1 Transcription Factor, a well known factor for regulating genes for Mitochondria Biogenesis, and a host of other biological functions. We are going to explore TSS upstream sequences and perform gene enrichment analysis. 

I hope you all have downloaded the human genome, and also have a working version of WSL/Linux/MacOS. Make sure the data is accessible in WSL (not in powershell). 

### Step 1

Convert human gene annotation file (`human_gene_annotation.tsv.gz`) into a bed file. Here, I am sharing a few lines of the bed file.

```
chr1    11869   11870   chr1@11869-11870|DDX11L2        .       +
chr1    12010   12011   chr1@12010-12011|DDX11L1        .       +
chr1    17436   17437   chr1@17436-17437|MIR6859-1      .       -
chr1    24886   24887   chr1@24886-24887|WASH7P .       -
```

Column description are below:

- Col1: Chromosome

- Col2: Transcription Start Site (TSS) 

- Col3: Col2 + 1

- Col4: Col1@Col2-Col3|<gene name>

- Col5: "." (place holder for score)

- Col6: Strand (+/-) 


Once you have the bed file - your goal is to extend the coordinates (Col2 or Col3) by 500 bases in a strand aware manner. Install `bedtools` and look at the `bedtools slop` command


In the `go_enrichment` environment, install `bedtools`. 

## Complete Workflow


### Common steps
DREG output is harder to convert into the BED file needed for gene intersection as it does'nt preserve the chromosome names, so I also used FIMO with the NRF1 PWM from the JASPAR database. Both approaches are included for completeness.

#### Step 1: Convert Gene Annotation to BED Format
**Input:** `1human_gene_annotation.tsv`  
**Output:** `3output.bed`  
**Tool:** Python script (`2bee.py`)  
**Purpose:** Convert TSV gene annotation into BED format with chromosomes, TSS coordinates, and gene names.

```bash
python 2bee.py
```

#### Step 2: Extend Promoter Regions (BED slop)
**Input:** `3output.bed`, `hg38.chrom.sizes`  
**Output:** `4bedtools.slop.output.bed`  
**Tool:** `bedtools slop`  
**Purpose:** Extend coordinates 500 bases upstream of TSS (strand-aware) to capture promoter regions where NRF1 may bind.

```bash
bedtools slop -s -i 3output.bed -g hg38.chrom.sizes -l 500 -r 0 > 4bedtools.slop.output.bed
```

#### Step 3: Extract Sequences from Promoter Regions
**Input:** `4bedtools.slop.output.bed`, `hg38.fa`  
**Output:** `5bedtools.getfasta.output.fa`  
**Tool:** `bedtools getfasta`  
**Purpose:** Extract DNA sequences from the extended promoter regions for motif scanning.

```bash
bedtools getfasta -s -fi hg38.fa -bed 4bedtools.slop.output.bed -fo 5bedtools.getfasta.output.fa
```

### Pipeline A: DREG (EMBOSS) motif search
Run these steps in `UsingDregForSequenceMatching/` after completing the common steps above.

```bash
cd UsingDregForSequenceMatching
```

#### Step 4A: Search for NRF1 binding sites with DREG
**Input:** `../5bedtools.getfasta.output.fa`  
**Output:** `6nrf1_hits.dreg`  
**Tool:** `dreg` (EMBOSS)  

```bash
dreg -sequence ../5bedtools.getfasta.output.fa -pattern GCGC..GCGC -outfile 6nrf1_hits.dreg
```

#### Step 5A: Convert DREG windows to BED
**Input:** `6nrf1_hits.dreg`, `../4bedtools.slop.output.bed`  
**Output:** `7nrf1_hits.bed`  
**Tool:** Python script (`7dreg_to_bed_for_intersect.py`)

```bash
ln -s ../4bedtools.slop.output.bed 4bedtools.slop.output.bed
python 7dreg_to_bed_for_intersect.py
```

#### Step 6A: Find genes with NRF1 hits (bedtools intersect)
**Input:** `../4bedtools.slop.output.bed`, `8nrf1_hits.bed`  
**Output:** `8nrf1_genes_intersected.bed`

```bash
bedtools intersect -a ../4bedtools.slop.output.bed -b 8nrf1_hits.bed -wa -u > 8nrf1_genes_intersected.bed
```

#### Step 7A: Extract gene names
**Input:** `9nrf1_genes_intersected.bed`  
**Output:** `10genes_with_nrf1.txt`

```bash
cut -f4 8.1nrf1_genes_intersected.bed | sed 's/.*|//' | sort -u > 10.1genes_with_nrf1.txt
```

#### Step 8A: Gene Ontology enrichment
**Input:** `10genes_with_nrf1.txt`  
**Output:** `11GO_annotation_results.csv`, `12nrf1_dotplot.pdf`

```bash
Rscript final_script.R
```

### Pipeline B: FIMO (MEME) motif search
Run these steps in `UsingFIMOforSequenceMatching/` after completing the common steps above.

```bash
cd UsingFIMOforSequenceMatching
```

#### Step 4B: Search for NRF1 binding sites with FIMO
**Input:** `../5bedtools.getfasta.output.fa`, `6nrf1.meme` (NRF1 PWM)  
**Output:** `7fimo_nrf1_results/fimo.txt`  
**Tool:** `fimo` (MEME Suite)

```bash
mamba install -c bioconda meme
# Download NRF1 PWM from JASPAR: https://jaspar.genereg.net/
fimo --oc 7fimo_nrf1_results 6nrf1.meme ../5bedtools.getfasta.output.fa
```

#### Step 5B: Convert FIMO results to BED
**Input:** `7fimo_nrf1_results/fimo.txt`  
**Output:** `8nrf1_absolute_hits.bed`

```bash
awk -F'\t' '!/^#/ && $1!="" {
    seq = $2;
    gsub(/\(.\)/, "", seq);
    split(seq, a, ":");
    chr = a[1];
    split(a[2], b, "-");
    offset = b[1];
    bed_start = offset + $3 - 1;
    bed_end = offset + $4;
    print chr "\t" bed_start "\t" bed_end "\t" $1 "\t" $6 "\t" $5
}' 7fimo_nrf1_results/fimo.txt > 8nrf1_absolute_hits.bed
```

#### Step 6B: Find genes with NRF1 binding sites (bedtools intersect)
**Input:** `../4bedtools.slop.output.bed`, `8nrf1_absolute_hits.bed`  
**Output:** `9nrf1_genes_intersected.bed`

```bash
bedtools intersect -a ../4bedtools.slop.output.bed -b 8nrf1_absolute_hits.bed -wa -u > 9nrf1_genes_intersected.bed
```

#### Step 7B: Extract gene names
**Input:** `9nrf1_genes_intersected.bed`  
**Output:** `10genes_with_nrf1.txt`

```bash
cut -f4 9nrf1_genes_intersected.bed | sed 's/.*|//' | sort -u > 10genes_with_nrf1.txt
```

#### Step 8B: Gene Ontology enrichment
**Input:** `10genes_with_nrf1.txt`  
**Output:** `11GO_annotation_results.csv`, `12nrf1_dotplot.pdf`

```bash
Rscript final_script.R
```

## Summary

The pipeline identifies genes regulated by NRF1 transcription factor and characterizes their biological roles through:
1. Genomic coordinate mapping (TSS → promoter regions)
2. Sequence extraction and motif detection (DREG or FIMO)
3. Gene-region mapping (bedtools intersect)
4. Functional enrichment (clusterProfiler in R)

This reveals that NRF1 preferentially regulates genes involved in mitochondrial biogenesis and energy metabolism, confirming known NRF1 biology.   
# GO_Enrichment
# GO_Enrichment
