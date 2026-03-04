# Bulk-Transcriptomics-Differential-Expression-Analysis
BINF 6110 Assignment 2.  Differential expression analysis, functional annotation, and communication of results using data from a study that examined yeast biofilm (velum) development during wine aging.

## Table of Contents

1. [Introduction](#introduction)
2. [Methods](#methods)
3. [Results](#results)
4. [Discussion](#discussion)
5. [References](#references)

---

## Introduction


---

## Methods

### Data

Raw single-end RNA-seq reads for nine samples (three replicates × three biofilm stages: Early, Thin, and Mature) were obtained from NCBI SRA (BioProject PRJNA592304; Mardanov et al., 2020). Sample accessions are listed below:

| Stage | Sample ID | SRA Accession |
|-------|-----------|---------------|
| Early | IL20 | SRR10551665 |
| Early | IL21 | SRR10551664 |
| Early | IL22 | SRR10551663 |
| Thin  | IL23 | SRR10551662 |
| Thin  | IL24 | SRR10551661 |
| Thin  | IL25 | SRR10551660 |
| Mature | IL29 | SRR10551659 |
| Mature | IL30 | SRR10551658 |
| Mature | IL31 | SRR10551657 |

### Download and Quality Control

Reads were downloaded from the NCBI SRA using `fasterq-dump` from the SRA Toolkit (v3.x) and compressed with `gzip`. Quality control was assessed with **FastQC** v0.12.1 (Andrews, 2010).

### Reference Transcriptome and Indexing

The *S. cerevisiae* reference transcriptome (GCF_000146045.2, R64 assembly) was downloaded from NCBI and indexed with **Salmon** v1.10.2 (Patro et al., 2017):

### Quantification

Transcript-level abundance was quantified for each sample using Salmon's quasi-mapping:

```bash
INDEX=~/scratch/BINF6110/.../data/salmon_index
QUANTS=~/scratch/BINF6110/.../data/quants

SAMPLES=(SRR10551665 SRR10551664 SRR10551663
         SRR10551662 SRR10551661 SRR10551660
         SRR10551659 SRR10551658 SRR10551657)

for SAMPLE in "${SAMPLES[@]}"; do
    salmon quant \
        -i $INDEX \
        -l A \
        -r $RAWDIR/${SAMPLE}.fastq.gz \
        --validateMappings \
        -p 8 \
        -o $QUANTS/${SAMPLE}
done
```

### Differential Expression Analysis in R


### Visualization


### Functional Enrichment



## Results

### Overall Data Structure

![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/PCA.png)
**Figure 1.** PCA of variance-stabilized counts across nine *S. cerevisiae* biofilm samples. Each point represents one biological replicate; colours indicate biofilm stage (Early = red, Thin = green, Mature = blue). PC1 and PC2 together explain 93% of variance.

---

### Differential Expression Analysis

DESeq2 identified substantial transcriptional remodelling across all three pairwise comparisons, with 5,963 genes detected with non-zero counts:

| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|-----------|-------------|---------------|
| Thin vs. Early | 1,126 | 594 | 532 |
| Mature vs. Early | 1,872 | 1,010 | 862 |
| Mature vs. Thin | 1,319 | 731 | 588 |


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/volcano%20plots.png)
**Figure 2.** Volcano plots for each pairwise comparison. Each point is a gene; maroon = significantly upregulated (log₂FC > 1, adjusted p < 0.05), navy = significantly downregulated, gray = not significant. Dashed vertical lines indicate |log₂FC| = 1. Top 5 up- and downregulated genes are labelled by adjusted p-value.

---

### Heatmap of Top Differentially Expressed Genes


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/Heatmap.png)
**Figure 3.** Heatmap of the 30 most statistically significant DEGs. Expression values are row-scaled variance-stabilized counts (VST). Column annotations indicate biofilm stage (Early = green, Thin = orange, Mature = dark red). Rows are clustered by hierarchical clustering; columns are ordered by stage.


---

### Functional Enrichment: Gene Ontology ORA


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/GO%20ORA.png)
**Figure 4.** GO Biological Process ORA dot plots for each pairwise comparison. Dot size reflects the gene ratio (proportion of DEGs in the GO term); dot colour indicates adjusted p-value. Up to 10 terms are shown per cluster after redundancy reduction (Jaccard similarity cutoff = 0.7).


---

### Functional Enrichment: KEGG Pathway ORA


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/KEGG%20ORA.png)
**Figure 5.** KEGG pathway ORA dot plots for each pairwise comparison. Dot size reflects gene ratio; dot colour indicates adjusted p-value. Up to 10 pathways are shown per cluster.

---

## Discussion


---

## References

Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Boyle, A. P., Davis, S., Bhatt, D. L., & et al. (2020). Optimization and evaluation of trimming algorithms for RNA-seq data. *NAR Genomics and Bioinformatics*, 2(3), lqaa068. https://doi.org/10.1093/nargab/lqaa068

Esteve-Zarzoso, B., Peris-Torán, M. J., García-Maiquez, E., Uruburu, F., & Querol, A. (2001). Yeast population dynamics during the fermentation and biological aging of Sherry wines. *Applied and Environmental Microbiology*, 67(5), 2056–2061.

Francois, J., & Parrou, J. L. (2001). Reserve carbohydrates metabolism in the yeast *Saccharomyces cerevisiae*. *FEMS Microbiology Reviews*, 25(1), 125–145.

Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for computing and annotating genomic ranges. *PLOS Computational Biology*, 9(8), e1003118.

Legras, J. L., Moreno-Garcia, J., Zara, S., Zara, G., Garcia-Martinez, T., Mauricio, J. C., Fidalgo, C., Gonzalez, P., Aranda, A., & Mas, A. (2016). Flor yeast: New perspectives beyond wine aging. *Frontiers in Microbiology*, 7, 503.

Lo, W. S., & Dranginis, A. M. (1998). The cell surface flocculin Flo11 is required for pseudohyphae formation and invasion by *Saccharomyces cerevisiae*. *Molecular Biology of the Cell*, 9(1), 161–171.

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

Mardanov, A. V., Beletsky, A. V., Ravin, N. V., Petrov, A. I., Откидач, Е. И., Golubev, W. I., & Morozova, V. V. (2020). Transcriptome of the *Saccharomyces cerevisiae* velum community during biological wine aging. *Frontiers in Microbiology*, 11, 538. https://doi.org/10.3389/fmicb.2020.00538

Nakagawa, Y., Sakumoto, N., Kaneko, Y., & Harashima, S. (2002). *Mga2p* is a putative sensor for low temperature and oxygen to induce OLE1 transcription in *Saccharomyces cerevisiae*. *Biochemical and Biophysical Research Communications*, 291(3), 707–713.

Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods*, 14(4), 417–419.

Reynolds, T. B., & Fink, G. R. (2001). Bakers' yeast, a model for fungal biofilm formation. *Science*, 291(5505), 878–881.

Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521.

Srivastava, A., Malik, L., Sarkar, H., Zakeri, M., Soneson, C., Love, M. I., Kingsford, C., & Patro, R. (2020). Alignment and mapping methodology influence transcript abundance estimation. *Genome Biology*, 21(1), 239.

Verstrepen, K. J., & Klis, F. M. (2006). Flocculation, adhesion and biofilm formation in yeasts. *Molecular Microbiology*, 60(1), 5–15.

Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: A revolutionary tool for transcriptomics. *Nature Reviews Genetics*, 10(1), 57–63.

Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer-Verlag New York.

Williams, C. R., Baccarella, A., Parrish, J. Z., & Kim, C. C. (2016). Trimming of sequence reads alters RNA-Seq gene expression estimates. *BMC Bioinformatics*, 17, 103.

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*, 2(3), 100141.

Yu, G. (2022). enrichplot: Visualization of functional enrichment result. R package version 1.16.2. https://github.com/YuLab-SMU/enrichplot

Zara, G., Budroni, M., Mannazzu, I., & Zara, S. (2010). Air-liquid biofilm formation is dependent on ammonium depletion in a *Saccharomyces cerevisiae* flor strain. *FEMS Yeast Research*, 10(5), 582–586.

Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: Removing the noise and preserving large differences. *Bioinformatics*, 35(12), 2084–2092.














