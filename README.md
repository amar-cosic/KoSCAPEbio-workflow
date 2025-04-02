# Paper Workflow and Data Repository

This repository contains all data files, scripts, and high-quality figure PDFs necessary to reproduce the analyses and visualizations in the paper: **KoSCAPEbio differentiates species of the Klebsiella oxytoca Complex in 16S rRNA Amplicon Data**. Each folder corresponds to a specific figure in the paper, containing the data and scripts used to generate the results. Additionally, high-quality PDFs of the figures are provided use in the publication.

# Repository Organization

## Main Folders and Figures

- `Figure_1/`: Contains all data and scripts used to generate **Figure 1** and **Supplementary Figure S1** related to *Phylogeny* analysis.  
  Contents: `csv`, `fasta`, `nwk`, `txt`, `tsv`, `py`

- `Figure_2/`: Contains all data and scripts used to generate **Figure 2** and **Supplementary Figure S4**, both related to *Alignment*.  
  Contents: `fasta`, `nwk`, `txt`, `py`

- `Figure_3/`: Contains all data and scripts used to generate **Figure 3** related to the *Mock* dataset.  
  Contents: `fasta`, `qza`, `png`, `tsv`, `biom`, `aln`, `fastq`, `py`, `sh`

- `Figure_4/`: Contains all data and scripts used to generate **Figure 4** and **Supplementary Figure S5**, related to *Clinical Datasets*.  
  Contents: `fasta`, `qza`, `tsv`, `txt`, `png`, `biom`, `xlsx`

- `Figure_5/`: Contains all data and scripts used to generate **Figure 5** and **Supplementary Figure S5**, related to *Sim Clinical Dataset*.  
  Contents: `xlsx`, `R`

## Supplementary Figures without Dedicated Folder

- `No Folder`:  
    - **Supplementary Figure S2**: *Pipeline Database*
    - **Supplementary Figure S3**: *Pipeline Search*

## PDFs

- `PDFs/`: High-quality PDF versions of Figures 1–5 and Supplementary Figures S1–S5 used in the publication.

## Reproducing Figures

The folders for each figure contain the relevant data and scripts needed for analysis and visualization. Each script is labeled according to what it does.

1. **Install Required Software**: Ensure you have the necessary software installed, such as Python or R, along with any specific libraries used in the scripts. Requirements are noted within each script.

2. **Adjust Paths as Needed**:
   - Absolute file paths (e.g., working directory paths) have been omitted from the scripts for privacy reasons. Before running the scripts, please set the working directory or update file paths as needed to match the locations of the downloaded data on your system.
   
3. **Run Scripts**:
   - Navigate to the appropriate `Figure_*` folder for the figure you wish to replicate.
   - Run the scripts to process the data and generate the figure components.

4. **PDF Access**: High-resolution PDFs of each figure are available in the `PDFs/` folder.

## Citation

If you use this repository, please cite our paper:

**KoSCAPEbio differentiates species of the Klebsiella oxytoca Complex in 16S rRNA Amplicon Data**  
Authors: Amar Cosic, Bettina Halwachs, Kristina Schild, Ellen L. Zechner, & Sabine Kienesberger 
Journal: [Journal Name, Year]  
DOI: [DOI or link to paper]

## Contact
For help and support, please contact:
-  **Name**: Amar Cosic
-  **Email**: [amar.cosic995@gmail.com](mailto:amar.cosic995@gmail.com)