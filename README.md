## Summary 
This repository contains scripts for the analysis and integration of a comprehensive aging epigenomic atlas across 27 organs from young and aged mice, including six histone modifications, chromatin accessibility, DNA methylation, three-dimensional genome architecture, and gene expression.

To visualize and explore the data, please visit: 
https://zhangyx.lab.westlake.edu.cn/MusEpiAgingBulk/

## Notes on documentation and support
Some parts of the scripts are lightly documented. Many were written for **exploratory analyses**, and not all of them contributed to the final results. That said, you are welcome to reuse any pieces you find useful for your own work.

If you have any questions about the code, please feel free to open an [issue](https://github.com/silentFUSU/Mouse_aging_epigenomic_Atlas/issues) or email me.

## Project directory structure
The project directory is organized as follows. The `code/` directory corresponds to the contents of this GitHub repository.

```text  
project MusEpiAgingBulk
├── code/    ***[everything in this GitHub repository]***
├── data/
│   ├── rawdata/
│   ├── public_data/
│   └── samples/
│       ├── tissue1/
│       │   ├── H3K27me3/
│       │   ├── H3K9me3/
│       │   ├── H3K36me3/
│       │   ├── H3K4me1/
│       │   ├── H3K4me3/
│       │   ├── H3K27ac/
│       │   ├── ATAC/
│       │   ├── RNA/
│       │   ├── WGBS/
│       │   └── HiC/
│       ├── ...
│       └── tissueN/
│           └── (same as above)
└── result/
    ├── analysis_1/
    ├── ...
    └── analysis_N/
```

## Read alignment and initial processing
Read alignment and basic preprocessing were carried out using our in-house pipelines.

* **CUT&Tag / ATAC-seq preprocessing:**
Paired-end CUT&Tag and ATAC-seq reads were adapter-trimmed with TrimGalore (v0.6.10) and aligned to mm10 (mouse), dm6 (fly), or ce11 (worm) using Bowtie2 (v2.4.1; `-X 2000`). PCR duplicates were removed with Picard (v2.9.3; `REMOVE_DUPLICATES=true`). Peaks for ATAC-seq, H3K27ac, H3K4me1, and H3K4me3 were called by age group using MACS2 (v2.2.7.1; merged replicates; `-f BAMPE --nomodel -q 0.0001 --keep-dup all`). Broad domains for H3K27me3, H3K9me3, and H3K36me3 were identified with SICER (v1.0.3; window 5000, gap 10000). H3K27me3 age-domains were further called using EDD (v1.1.19; `--fdr 0.05 -n 50000 --gap-penalty 80`).

* **RNA-seq preprocessing:**
RNA-seq reads were trimmed with fastp (v0.23.2), aligned to mm10 with STAR (v2.7.10a), and quantified with featureCounts (v2.0.1). Differential expression between age groups was performed using edgeR (v3.40.2); DEGs were defined as FDR < 0.05 and |log2(old/young)| > 1.2.

* **WGBS preprocessing:**
WGBS reads were trimmed with TrimGalore (v0.6.10; `--length 20 --gzip --clip_R2 10`) and aligned to bisulfite-converted mm10 using Bismark (v0.24.1; Bowtie2 v2.4.1). Duplicates were removed and methylation was extracted using Bismark utilities. Differentially methylated regions were identified with DSS (v2.46.0; `minCG=3, delta=0.1, p.threshold=0.01`).

* **Hi-C preprocessing:**
Hi-C data were processed using HiC-Pro (v3.1.0) with default settings and mapped to mm10. Valid pairs were converted to `.hic` files (hicpro2juicebox) for visualization in Juicebox (v2.15) and downstream analysis (HOMER v5.0.1). Unless noted otherwise, contact matrices were analyzed at 10-kb resolution and ICE-normalized.

## Subfolder Information
* **Archive**: General-purpose utility functions.
* **call_peak_bin**: Peak calling and fixed-bin generation for CUT&Tag and ATAC-seq, including feature quantification at both the peak and bin levels.
* **chromHMM**: ChromHMM state definition and analyses of state changes (Figure 1).
* **Differential_analysis**: Differential analysis for CUT&Tag and ATAC-seq at both peak-level and bin-level resolutions.
* **Figures**: Code to generate selected integrative figures used in the manuscript.
* **GRN**: Gene regulatory network (GRN) analysis (Figure 6).
* **MEF**: Analyses related to mouse embryonic fibroblast (MEF) datasets.
* **public_data**: Exploratory mining and analysis of public datasets used during the data-mining phase of the study.
* **quality_control**: Data quality control (QC) scripts.
* **specific_modality**: Modality-specific analyses tailored to each data type.
* **specific_tissue**: Tissue-specific analyses tailored to individual organs and tissues.

## Data: 
All datasets generated in this study will be released via the [China National Center for Bioinformation (CNCB)](https://www.cncb.ac.cn/)