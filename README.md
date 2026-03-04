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

*Saccharomyces cerevisiae* is a yeast strain that produces a surface biofilm known as a velum during the biological aging of certain wines, particularly sherry wines. This biofilm develops on wine surfaces that are open to the air and play a crucial role in shaping the properties of the final product by consuming ethanol and producing acetaldehyde and acetals under oxidative, nutrient-deficient conditions (Mardanov et al., 2020; Legras et al., 2016).

The velum undergoes distinct morphological stages, **Early**, **Thin**, and **Mature**, each associated with different metabolic demands. Understanding the transcriptional basis of these transitions is valuable both for wine making and for research into *S. cerevisiae* stress adaptation, carbon source switching, and biofilm architecture. As glucose is progressively depleted from Early to Mature biofilm (Mardanov et al., 2020), cells are expected to shift from fermentative to oxidative metabolism, a transition controlled by large-scale changes in gene expression.

Bulk RNA sequencing (RNA-seq) is a method used for measuring differential gene expression across conditions. RNA-seq quantifies the abundance of all transcripts in a sample simultaneously, providing a genome-wide snapshot of the transcriptome (Wang et al., 2009). Here, **Salmon** was used for quasi-mapping-based quantification (Patro et al., 2017), which employs selective alignment against the reference transcriptome to produce accurate transcript-level counts without the computational cost of full alignment. This why it was selected over splice aware aligners like **STAR** or **minimap2** which take longer to run and require an extra step to get gene expression counts. Counts were then imported into R with **tximport** (Soneson et al., 2015) for transcript-to-gene aggregation, and differential expression was performed using **DESeq2** (Love et al., 2014).

DESeq2 was selected because it is the most widely adopted tool for bulk transcriptomics DE analysis, implements the median-of-ratios normalization method to account for differences in RNA composition across samples, and applies empirical Bayes shrinkage of log-fold change (LFC) estimates, which is particularly important for low-count genes common in studies with only three replicates per group (Love et al., 2014). LFC shrinkage was applied using **apeglm** (Zhu et al., 2019) to reduce noise in effect size estimates. Functional enrichment was performed using **clusterProfiler** (Wu et al., 2021) for both Gene Ontology (GO) Biological Process Over-Representation Analysis (ORA) and KEGG pathway ORA, allowing biologically interpretable categorization of differentially expressed gene (DEG) sets.


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

Reads were downloaded from the NCBI SRA using `fasterq-dump` from the SRA Toolkit v3.0.9 and compressed with `gzip`. Quality control was assessed with **FastQC** v0.12.1 (Andrews, 2010).

### Reference Transcriptome and Indexing

The *S. cerevisiae* reference transcriptome (GCF_000146045.2, R64 assembly) was downloaded from NCBI and indexed with **Salmon** v1.10.2 (Patro et al., 2017):

### Quantification

Transcript-level abundance was quantified for each sample using Salmon:

```bash
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

Salmon output was imported into R using **tximport** v1.38.2  (Soneson et al., 2015). A transcript-to-gene mapping was generated from the R64 GTF file using `makeTxDbFromGFF` from the **GenomicFeatures** v1.62.0 package (Lawrence et al., 2013). Samples were modelled with a single-factor design (`~stage`) in **DESeq2** v1.50.2 (Love et al., 2014), with `Early` as the reference level. To obtain the Mature vs. Thin contrast, the reference was relevelled to `Thin` and DESeq2 was re-run. Pairwise contrasts were extracted with `results()` at α = 0.05. Log-fold change shrinkage was applied with **apeglm** using `lfcShrink()` (Zhu et al., 2018). DEGs were defined as genes with |log₂FC| > 1 and adjusted p-value < 0.05 (Benjamini–Hochberg correction).

```r
dds <- DESeqDataSetFromTximport(txi, sample_info, ~stage)
dds <- DESeq(dds)

LFC_tve <- lfcShrink(dds, coef = "stage_Thin_vs_Early", type = "apeglm")
```

### Visualization

**PCA** was performed on variance-stabilized counts (VST) using `plotPCA` from DESeq2. **Volcano plots** were generated with **ggplot2** v4.0.2 (Wickham, 2016) and **ggrepel** v0.9.7, labelling the top 5 up- and down-regulated genes (by adjusted p-value) in each comparison. A **heatmap** of the 30 most significant genes (minimum adjusted p-value across all three comparisons) was produced with **pheatmap** v1.0.13  using row-scaled VST counts.

### Functional Enrichment

Gene Ontology Biological Process (GO BP) and KEGG pathway ORA were performed with `compareCluster` from **clusterProfiler** v4.18.4 (Wu et al., 2021), testing upregulated and downregulated DEGs separately against the background of all expressed genes. GO results were simplified with a Jaccard similarity cutoff of 0.7 using `clusterProfiler::simplify`. Dot plots showing the top 10 enriched terms per cluster were generated with `dotplot` from **enrichplot** v1.30.4 (Yu, 2022). KEGG queries used organism code `sce` (*S. cerevisiae*). Multiple testing correction was performed using the Benjamini–Hochberg method (p.adjust < 0.05, q-value < 0.2).

## Results

### Overall Data Structure

Principal component analysis of variance-stabilized read counts revealed clear separation of the three biofilm stages (Figure 1). PC1 (67% variance explained) separated the Early and Mature stages at opposite extremes, consistent with the large distance between the earliest and latest stages of biofilm development. Thin biofilm samples clustered in an intermediate position along PC1 but were clearly separated from both Early and Mature samples along PC2 (26% variance explained), suggesting Thin biofilm has a partially distinct transcriptional profile rather than just being a intermediate. Within-group replicates clustered tightly, indicating high reproducibility across biological triplicates and supporting the validity of the downstream differential expression analysis.

![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/PCA2.png)

**Figure 1.** PCA of variance-stabilized counts across nine *S. cerevisiae* biofilm samples. Each point represents one biological replicate; colours indicate biofilm stage (Early = red, Thin = green, Mature = blue). PC1 and PC2 together explain 93% of variance.

---

### Differential Expression Analysis

DESeq2 identified substantial transcriptional remodelling across all three pairwise comparisons, with 5,963 genes detected with non-zero counts. The number of DEGs (|log₂FC| > 1, adjusted p < 0.05) increased with developmental distance between the compared stages:

| Comparison | Total DEGs | Upregulated | Downregulated |
|------------|-----------|-------------|---------------|
| Thin vs. Early | 1,126 | 594 | 532 |
| Mature vs. Early | 1,872 | 1,010 | 862 |
| Mature vs. Thin | 1,319 | 731 | 588 |


The largest number of DEGs was observed in the Mature vs. Early comparison (1,872 genes; ~31% of the detected transcriptome), consistent with the two stages being sampled at the greatest temporal distance (38 vs. 109 days post-inoculation; Mardanov et al., 2020). Volcano plots confirmed both the magnitude and statistical significance of these changes (Figure 2). In the Mature vs. Early comparison, genes *YJL052W*, *YGL055W*, *YGR087C*, and *YHR094C* were among the most significantly downregulated (navy, left of dashed lines), while *YIR019C*, *YNR073C*, *YNR072W*, and *YNR071C* were among the most significantly upregulated (maroon, right of dashed lines). In the Thin vs. Early comparison, *YGR087C* and *YHR094C* showed high downregulation significance, while *YCR105W*, *YKR097W*, and *YNR073C* were top upregulated genes. The Mature vs. Thin comparison showed fewer extreme outliers, with *YJR152W* and *YGL055W* most significantly downregulated and *YBR117C*, *YPR127W*, *YDR085C*, and *YEL070W* most significantly upregulated.

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














