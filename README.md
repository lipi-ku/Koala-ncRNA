# Supplementary Data and Analysis Scripts for Koala RNA-seq Dataset PRJEB21505

## Overview

This repository contains supplementary data, RNA-seq analysis outputs, sample metadata, and R scripts associated with the manuscript entitled:

**“Noncoding RNA expression profiles associated with koala retrovirus subtype composition: an exploratory transcriptomic study.”**

The publicly available koala (*Phascolarctos cinereus*) RNA-seq dataset PRJEB21505 was independently analyzed to assess the reproducibility of candidate noncoding RNA (ncRNA) expression patterns associated with koala retrovirus (KoRV) subtype composition and to investigate expression associations between candidate ncRNAs and selected genomic-context genes.

## Data Source

The original RNA-seq data are publicly available through the European Nucleotide Archive (ENA).

- BioProject accession: PRJEB21505
- Organism: *Phascolarctos cinereus*
- Public dataset: https://www.ebi.ac.uk/ena/browser/view/PRJEB21505

## Sample Information

The metadata file contains information for all 29 sequencing accessions. Of these, 22 samples were assigned to five dataset-defined groups and included in the downstream comparative expression analyses.

The ABC-versus-A comparison was designated as the primary replicated subtype comparison. The AB and AC groups were interpreted descriptively because each contained only one biological sample.

## Repository Contents

This repository includes the following supplementary resources:

- **Sample metadata:** Sequencing accession numbers, available sample characteristics, and KoRV subtype classifications.
- **RNA-seq analysis outputs:** Gene and ncRNA expression data, count matrices, normalized expression values, and differential expression analysis results.
- **Candidate ncRNA analysis:** Expression profiles of the nine selected candidate ncRNAs across KoRV subtype groups.
- **Genomic-context gene analysis:** Expression profiles of eight selected genes associated with candidate ncRNA genomic regions.
- **Correlation analysis:** Spearman correlation coefficients, nominal p-values, and false discovery rate (FDR)-adjusted p-values for ncRNA–gene expression associations.
- **R scripts:** Scripts used for downstream statistical analyses, multiple-testing correction, data processing, and graphical visualization.
- **Figure outputs:** Visualizations supporting the independent RNA-seq analysis presented in the manuscript.

## Statistical Analysis and Reproducibility

The downstream analyses included differential expression assessment using DESeq2, counts-per-million (CPM)-based expression comparisons, Wilcoxon rank-sum tests, and Spearman rank correlation analyses. Multiple-testing correction was performed using the Benjamini–Hochberg procedure.

The accompanying data and R scripts are provided to facilitate independent verification and reproduction of the reported downstream statistical analyses and graphical outputs, including Figure 7 of the manuscript.

## Interpretation and Limitations

The independent RNA-seq analysis was conducted as an exploratory assessment of the reproducibility of candidate ncRNA expression patterns. The unequal sample sizes, limited biological replication in certain subtype groups, and potential differences in host and sample characteristics should be considered when interpreting the results.

The observed ncRNA–gene expression correlations represent statistical associations and do not establish direct molecular interactions, cis-regulation, or causal relationships.

## Data Availability

The original sequencing data remain publicly accessible through ENA under accession PRJEB21505. The processed expression data, sample annotations, statistical analysis outputs, and accompanying R scripts are provided in this repository to support research transparency and reproducibility.

## Citation

Users of these supplementary resources should cite the original PRJEB21505 dataset publication and the associated manuscript by Akter et al. in *npj Viruses*.
