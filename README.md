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

*Saccharomyces cerevisiae* is a yeast strain that produces a surface biofilm known as a velum during the aging of certain wines, particularly sherry wines. This biofilm develops on wine surfaces that are open to the air and play a crucial role in shaping the properties of the final product by consuming ethanol and producing acetaldehyde and acetals under oxidative, nutrient-deficient conditions [^9][^5].

The velum undergoes distinct morphological stages, **Early**, **Thin**, and **Mature**, each associated with different metabolic demands. Understanding the transcriptional basis of these transitions is valuable both for wine making and for research into *S. cerevisiae* stress adaptation, carbon source switching, and biofilm architecture. As glucose is progressively depleted from Early to Mature biofilm[^9], cells are expected to shift from fermentative to oxidative metabolism, a transition controlled by large-scale changes in gene expression.

Bulk RNA sequencing (RNA-seq) is a method used for measuring differential gene expression across conditions. RNA-seq quantifies the abundance of all transcripts in a sample simultaneously, providing a genome-wide snapshot of the transcriptome.[^16] Here, **Salmon** was used for quasi-mapping-based quantification[^12], which employs selective alignment against the reference transcriptome to produce accurate transcript-level counts without the computational cost of full alignment. This why it was selected over splice aware aligners like **STAR** or **minimap2** which take longer to run and require an extra step to get gene expression counts. Counts were then imported into R with **tximport**[^7] for transcript-to-gene aggregation, and differential expression was performed using **DESeq2**.[^7]

DESeq2 was selected because it is the most widely adopted tool for bulk transcriptomics DE analysis, implements the median-of-ratios normalization method to account for differences in RNA composition across samples, and applies empirical Bayes shrinkage of log-fold change (LFC) estimates, which is particularly important for low-count genes common in studies with only three replicates per group.[^7] LFC shrinkage was applied using **apeglm**[^21] to reduce noise in effect size estimates. Functional enrichment was performed using **clusterProfiler**[^17] for both Gene Ontology (GO) Biological Process Over-Representation Analysis (ORA) and KEGG pathway ORA, allowing biologically interpretable categorization of differentially expressed gene (DEG) sets.


---

## Methods

### Data

Raw single-end RNA-seq reads for nine samples (three replicates × three biofilm stages: Early, Thin, and Mature) were obtained from NCBI SRA (BioProject PRJNA592304[^9]). Sample accessions are listed below:

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


## Additional Metadata from Mardanov et al. 2020

| Sample designation |Early biofilm | Thin biofilm | Mature biofilm |
|--------------------|--------------|--------------|----------------|
| Days from inoculation to sampling | 38 | 83 | 109 |
| Ethanol% (v/v) | 12.4 | 10.8 | 9.6 |
| Volatile acidity (g/l) | 0.2 | 0.2 | 0.1 |
| Total acidity (g/l) | 7.8 | 7.4 | 7.0 | 
| Aldehydes (mg/l) | 382.8 | 531.3 | 668.8 |
| Acetals (mg/l) | 147.2 | 253.7 | 280.3 |
| pH | 3.6 | 3.6 | 3.6 |
| Glucose (g/l) | 0.2 | 0.1 | <0.1 |
| Fructose (g/l) | <0.1 | <0.1 | <0.1 | 
| Oligosaccharides (g/l) | 0.3 | 0.2 | 0.2 |
| Glycerol (g/l) | 8.5 | 7.9 | 6.8 |

### Download and Quality Control

Reads were downloaded from the NCBI SRA using `fasterq-dump` from the SRA Toolkit v3.0.9 and compressed with `gzip`. Quality control was assessed with **FastQC** v0.12.1.[^1]

### Reference Transcriptome and Indexing

The *S. cerevisiae* reference transcriptome (GCF_000146045.2, R64 assembly) was downloaded from NCBI and indexed with **Salmon** v1.10.2[^12]:

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

Salmon output was imported into R using **tximport** v1.38.2.[^14][^7] A transcript-to-gene mapping was generated from the R64 GTF file using `makeTxDbFromGFF` from the **GenomicFeatures** v1.62.0 package.[^4] Samples were modelled with a single-factor design (`~stage`) in **DESeq2** v1.50.2 (Love et al., 2014), with `Early` as the reference level. To obtain the Mature vs. Thin contrast, the reference was relevelled to `Thin` and DESeq2 was re-run. Pairwise contrasts were extracted with `results()` at α = 0.05. Log-fold change shrinkage was applied with **apeglm** using `lfcShrink()` [^21]. DEGs were defined as genes with |log₂FC| > 1 and adjusted p-value < 0.05 (Benjamini–Hochberg correction).

```r
dds <- DESeqDataSetFromTximport(txi, sample_info, ~stage)
dds <- DESeq(dds)

LFC_tve <- lfcShrink(dds, coef = "stage_Thin_vs_Early", type = "apeglm")
```

### Visualization

**PCA** was performed on variance-stabilized counts (VST) using `plotPCA` from DESeq2. **Volcano plots** were generated with **ggplot2** v4.0.2[^17] and **ggrepel** v0.9.7, labelling the top 5 up- and down-regulated genes (by adjusted p-value) in each comparison. A **heatmap** of the 30 most significant genes (minimum adjusted p-value across all three comparisons) was produced with **pheatmap** v1.0.13  using row-scaled VST counts.

### Functional Enrichment

Gene Ontology Biological Process (GO BP) and KEGG pathway ORA were performed with `compareCluster` from **clusterProfiler** v4.18.4[^18], testing upregulated and downregulated DEGs separately against the background of all expressed genes. GO results were simplified with a Jaccard similarity cutoff of 0.7 using `clusterProfiler::simplify`. Dot plots showing the top 10 enriched terms per cluster were generated with `dotplot` from **enrichplot** v1.30.4.[^19] KEGG queries used organism code `sce` (*S. cerevisiae*). Multiple testing correction was performed using the Benjamini–Hochberg method (p.adjust < 0.05, q-value < 0.2).

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

Hierarchical clustering of the 30 most significant genes (ranked by minimum adjusted p-value across all comparisons) revealed two dominant expression modules (Figure 3). One cluster of genes, from *YOR247W* down to *YPL106C* was progressively upregulated from Early to Mature, reaching peak expression in the Mature stage. A second cluster, from *YJR152W* down to *YJL212C* showed the inverse pattern, with high expression in Early biofilm that declined through Thin and Mature stages. A third smaller cluster, including *YPR127W*, YBR117C*, and YDR085C* were specifically upregulated in the Mature stage relative to both Early and Thin. Finally, three genes, *YKR097W*, *YDL243C*, and especially *YCR105W*, were significantly upregulated in the thin stage.  *YCR105W* is especially interesting as it is only upregulated at the thin stage. 

![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/Heatmap.png)
**Figure 3.** Heatmap of the 30 most statistically significant DEGs. Expression values are row-scaled variance-stabilized counts (VST). Column annotations indicate biofilm stage (Early = green, Thin = orange, Mature = dark red). Rows are clustered by hierarchical clustering; columns are ordered by stage.


---

### Functional Enrichment: Gene Ontology ORA

In the **Thin vs. Early** comparison, upregulated genes were most significantly enriched for **cytoplasmic translation** (p.adjust < 1×e^-10), as well as **cellular respiration** and **proton transmembrane transport**. Downregulated genes were enriched for carbohydrate catabolism terms including **glycolytic process**, **monocarboxylic acid metabolic process**, and **small molecule catabolic process**. This suggests that as cells transition from Early to Thin biofilm, they upregulate mitochondrial oxidative capacity and ribosome biogenesis while reducing fermentation[^24].

In the **Mature vs. Early** comparison, upregulated genes were enriched for **mitochondrion organization**, **mitochondrial translation**, and **mitochondrial gene expression**, terms absent or less prominent in the Thin vs. Early comparison, suggesting a stronger commitment to oxidative metabolism by the Mature stage. Downregulated genes in this comparison showed particularly large gene ratios and low p.adj for **transmembrane transport**, **monocarboxylic acid metabolic process**, and **lipid metabolic process**, consistent with almost complete glucose exhaustion at 109 days and shutdown of fermentation pathways as evident in the metata provided by Mardanov et al. (2020)

In the **Mature vs. Thin** comparison, upregulated genes were enriched for **carbohydrate biosynthetic process**, **energy reserve metabolic process**, and **response to abiotic stimulus**, the latter consistent with oxidative stress adaptation at the biofilm surface.[^24] Downregulated genes showed enrichment for **monocarboxylic acid metabolic process**, **transmembrane transport**, **lipid biosynthetic process**, and a handful more.


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/GO%20ORA.png)
**Figure 4.** GO Biological Process ORA dot plots for each pairwise comparison. Dot size reflects the gene ratio (proportion of DEGs in the GO term); dot colour indicates adjusted p-value. Up to 10 terms are shown per cluster after redundancy reduction (cutoff = 0.7).


---

### Functional Enrichment: KEGG Pathway ORA

KEGG pathway ORA confirmed and extended the GO results (Figure 5). Across comparisons, the **Ribosome** pathway was the most significantly enriched among upregulated genes (p.adjust < 1×10⁻¹⁷ in Thin vs. Early; p.adjust < 1×10⁻⁹ in Mature vs. Early), consistent with GO results showing cytoplasmic translation enrichment. The **Citrate cycle (TCA cycle)** was enriched in the upregulated genes for both Thin vs. Early and Mature vs. Early comparisons, reflecting the shift toward aerobic respiration.

Among downregulated genes, the **Biosynthesis of secondary metabolites** pathway was consistently and strongly enriched across all comparisons. **Glycolysis / Gluconeogenesis**, **Carbon metabolism**, **Fatty acid metabolism**, and **Biosynthesis of amino acids** were enriched among downregulated genes in the Mature vs. Early comparison, suggesting the slowing of growth, sugar fermentation, and amino acid production. 

In the **Mature vs. Thin** comparison, **Starch and sucrose metabolism** emerged as the top upregulated pathway, while **Biosynthesis of secondary metabolites** remained the top downregulated pathway. The upregulation of carbohydrate storage pathways (including glycogen and trehalose biosynthesis from the GO analysis) at the Mature stage may represent a stress survival strategy under extremely low carbon availability [^3].


![Alt text for the image](https://github.com/mahnoor-nizz/Bulk-Transcriptomics-Differential-Expression-Analysis/blob/main/Figures/KEGG%20ORA.png)
**Figure 5.** KEGG pathway ORA dot plots for each pairwise comparison. Dot size reflects gene ratio; dot colour indicates adjusted p-value. Up to 10 pathways are shown per cluster.

---

## Discussion

This analysis reveals that the coordinated and progressive transcriptional reprogramming that occurs during *S. cerevisiae* velum biofilm development is consistent with a shift from fermentative to oxidative metabolism as glucose is depleted and ethanol accumulates over the 109-day aging period.

The most note worthy change across comparisons is the downregulation of fermentation and glycolysis pathways alongside upregulation of the TCA cycle, oxidative phosphorylation, and ribosomal biogenesis. The velum lifestyle is known to be oxidative: *S. cerevisiae* strains that are capable of velum formation preferentially oxidize ethanol via the glyoxylate cycle and gluconeogenesis under nutrient limiting conditions [^2] [^5]. The strong enrichment of **mitochondrial organization**, **mitochondrial translation** and **mitochondrial respiratory chain complex assembly** in upregulated genes of the Mature vs. Early comparison is consistent with published observations of increased mitochondrial content in biofilm-forming flor strains [^20].

Among the most consistently upregulated individual genes in thin and mature stages shown in the volcano plots and heatmap is *YIR019C* (*FLO11/MUC1*). *FLO11* is particularly notable as it encodes a GPI-anchored cell surface glycoprotein involved in flocculation, adhesion, making it a master regulator of yeast biofilm formation [^6][^13][^25]. Its consistent upregulation from Thin through Mature biofilm stages supports the role of Flo11 in maintaining biofilm structural integrity and yeast colony survival. The upregulation of *FLO11* is regulated in part by the cAMP-PKA and MAPK pathways in response to nutrient limitation, consistent with the glucose depletion observed across stages[^15].

Other notable genes are *YHR094C*(*HXT1/HOR4*), and *YPL095C*(*EEB1*), which the heatmap shows are significantly upregulated during the early stage. *HXT1* is a glucose transporter whose expression is induced in the presence of glucose and repressed when glucose is limiting, and *EEB1* is responsible for the major part of medium-chain fatty acid ethyl ester biosynthesis during fermentation[^24][^25].

Among the most significantly downregulated genes in mature vs early and mature vs thin are *YJL052W* (*TDH1*), and *YGL055W* (*OLE1*). *TDH1* encodes glyceraldehyde-3-phosphate dehydrogenase, a key glycolytic enzyme, whose downregulation is directly consistent with reduced glycolytic flux[^10][^25]. *OLE1* encodes a fatty acid desaturase and its downregulation in Mature biofilm may reflect altered membrane lipid composition under ethanol stress[^11][^25]. The strong and consistent enrichment of **Biosynthesis of secondary metabolites** among downregulated genes likely reflects coordinated reduction of anabolic pathways under carbon and energy limitation.

The upregulation of **trehalose metabolism** and **response to abiotic stimulus** in the Mature vs. Thin comparison is consistent with the known role of trehalose as a stress protectant in yeast, particularly under ethanol and oxidative stress[^3]. Together, these results suggest that Mature biofilm cells have substantially rewired their metabolism to survive a hostile environment characterized by low glucose, high ethanol (~9.6% v/v), elevated acetaldehyde (~669 mg/L), and oxidative conditions at the air-liquid interface[^9].

**Limitations:** This study used only three replicates per stage, which limits statistical power, particularly for low-abundance transcripts. No quality trimming was applied prior to quantification, which is consistent with current recommendations for gene expression studies [^8], but trimming would remain important for transcript assembly. Future studies could incorporate a Likelihood Ratio Test (LRT) to model the full developmental time course rather than pairwise contrasts, which would better capture continuously varying expression [^7].

---

## References

[^1]:	Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[^2]:	Esteve-Zarzoso, B., Peris-Torán, M. J., García-Maiquez, E., Uruburu, F., & Querol, A. (2001). Yeast population dynamics during the fermentation and biological aging of Sherry wines. Applied and Environmental Microbiology, 67(5), 2056–2061. https://doi.org/10.1128/AEM.67.5.2056-2061.2001

[^3]:	Francois, J., & Parrou, J. L. (2001). Reserve carbohydrates metabolism in the yeast Saccharomyces cerevisiae. FEMS Microbiology Reviews, 25(1), 125–145. https://doi.org/10.1111/j.1574-6976.2001.tb00574.x

[^4]:	Lawrence, M., Huber, W., Pagès, H., Aboyoun, P., Carlson, M., Gentleman, R., Morgan, M. T., & Carey, V. J. (2013). Software for computing and annotating genomic ranges. PLOS Computational Biology, 9(8), e1003118. https://doi.org/10.1371/journal.pcbi.1003118

[^5]:	Legras, J. L., Moreno-Garcia, J., Zara, S., Zara, G., Garcia-Martinez, T., Mauricio, J. C., Fidalgo, C., Gonzalez, P., Aranda, A., & Mas, A. (2016). Flor yeast: New perspectives beyond wine aging. Frontiers in Microbiology, 7, 503. https://doi.org/10.3389/fmicb.2016.00503

[^6]:	Lo, W. S., & Dranginis, A. M. (1998). The cell surface flocculin Flo11 is required for pseudohyphae formation and invasion by Saccharomyces cerevisiae. Molecular Biology of the Cell, 9(1), 161–171. https://doi.org/10.1091/mbc.9.1.161

[^7]:	Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

[^8]:	Liao, Y., & Shi, W. (2020). Read trimming is not required for mapping and quantification of RNA-seq reads at the gene level. NAR Genomics and Bioinformatics, 2(3), lqaa068. https://doi.org/10.1093/nargab/lqaa068

[^9]:	Mardanov, A. V., Beletsky, A. V., Ravin, N. V., Petrov, A. I., Golubev, W. I., & Morozova, V. V. (2020). Transcriptome of the Saccharomyces cerevisiae velum community during biological wine aging. Frontiers in Microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538

[^10]:	McAlister, L., & Holland, M. J. (1985). Isolation and characterization of yeast strains carrying mutations in the glyceraldehyde-3-phosphate dehydrogenase genes. Journal of Biological Chemistry, 260(28), 15013–15018.

[^11]:	Nakagawa, Y., Sakumoto, N., Kaneko, Y., & Harashima, S. (2002). Mga2p is a putative sensor for low temperature and oxygen to induce OLE1 transcription in Saccharomyces cerevisiae. Biochemical and Biophysical Research Communications, 291(3), 707–713. https://doi.org/10.1006/bbrc.2002.6507

[^12]:	Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

[^13]:	Reynolds, T. B., & Fink, G. R. (2001). Bakers' yeast, a model for fungal biofilm formation. Science, 291(5505), 878–881. https://doi.org/10.1126/science.291.5505.878

[^14]:	Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences. F1000Research, 4, 1521. https://doi.org/10.12688/f1000research.7563.2

[^15]:	Verstrepen, K. J., & Klis, F. M. (2006). Flocculation, adhesion and biofilm formation in yeasts. Molecular Microbiology, 60(1), 5–15. https://doi.org/10.1111/j.1365-2958.2006.05072.x

[^16]:	Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: A revolutionary tool for transcriptomics. Nature Reviews Genetics, 10(1), 57–63. https://doi.org/10.1038/nrg2484

[^17]:	Wickham, H. (2016). ggplot2: Elegant graphics for data analysis. Springer-Verlag.

[^18]:	Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141

[^19]:	Yu, G. (2022). enrichplot: Visualization of functional enrichment result (R package version 1.16.2). https://github.com/YuLab-SMU/enrichplot

[^20]:	Zara, G., Budroni, M., Mannazzu, I., & Zara, S. (2010). Air-liquid biofilm formation is dependent on ammonium depletion in a Saccharomyces cerevisiae flor strain. FEMS Yeast Research, 10(5), 582–586. https://doi.org/10.1111/j.1567-1364.2010.00634.x

[^21]:	Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: Removing the noise and preserving large differences. Bioinformatics, 35(12), 2084–2092. https://doi.org/10.1093/bioinformatics/bty895

[^22]:	Ozcan, S., & Johnston, M. (1999). Function and regulation of the Saccharomyces cerevisiae HXT family of sugar transporters. Microbiology and Molecular Biology Reviews, 63(3), 554–569. https://doi.org/10.1128/MMBR.63.3.554-569.1999  

[^23]:	Saerens, S. M. G., Verstrepen, K. J., Van Laere, S. D. M., Voet, A. R. D., Van Dijck, P., Delvaux, F. R., & Thevelein, J. M. (2006). The Saccharomyces cerevisiae EHT1 and EEB1 genes encode novel enzymes with medium-chain fatty acid ethyl ester synthesis and hydrolysis capacity. Journal of Biological Chemistry, 281(7), 4446–4456. https://doi.org/10.1074/jbc.M512028200

[^24]:  Moreno-García, J., Mauricio, J. C., Moreno, J., & García-Martínez, T. (2017). Differential Proteome Analysis of a Flor Yeast Strain under Biofilm Formation. International Journal of Molecular Sciences, 18(4), 720. https://doi.org/10.3390/ijms18040720

[^25]: https://www.yeastgenome.org/






