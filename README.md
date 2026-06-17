# MIMUW Genomic Algorithms

Solutions, experiments and notes from the **Algorithms for Genomic Data Analysis** course at the Faculty of Mathematics, Informatics and Mechanics, University of Warsaw.

The repository focuses on algorithmic methods used in genomic data processing, especially problems related to high-throughput sequencing data. It includes implementations and experiments around read mapping, indexing, approximate matching, metagenomic classification, k-mer analysis, motif finding and genome assembly.

## Main assignments

### Read mapping

Implementation of a read mapping pipeline based on reference indexing, candidate generation and edit-distance verification.

The solution explores:

* suffix-array based indexing,
* exact and approximate shingle matching,
* candidate occurrence aggregation,
* edit-distance based verification,
* testing on sample and larger datasets,
* additional scripts for inspecting and visualizing matching behavior.

### Metagenomics

Implementation of a metagenomic classification assignment focused on identifying organisms from sequencing data.

The repository includes the main classifier implementation together with helper scripts and experiments used during solution development.

## Additional coursework

The lab folders contain smaller implementations and exercises covering selected topics from the course, including:

* text indexing,
* Burrows-Wheeler Transform,
* approximate pattern matching,
* high-error read mapping,
* k-mer spectra,
* randomized motif finding,
* genome assembly.

## Technologies

* Python
* Jupyter Notebook
* Shell scripts

## Notes

This is a coursework repository, not a production bioinformatics library. The goal was to implement, test and reason about algorithmic approaches to genomic data analysis, with emphasis on correctness, complexity and practical experimentation.
