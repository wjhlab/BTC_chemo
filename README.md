# IMC Analysis for Manuscript 

This repository contains the analytical pipeline, configuration files, and custom R scripts used for the Imaging Mass Cytometry (IMC) data analysis presented in:

> **3D Multi-Omic Analysis Reveals Signatures of Response and Resistance to Chemotherapy in Pancreatic Cancer**

---

## Project Overview
This project performs high-dimensional single-cell analysis on IMC data. The workflow includes data normalization, clustering, spatial analysis, and visualization as described in the associated manuscript.

## Repository Contents
This repo is organized to ensure reproducibility of the findings:

* **cvsphenographpipeline_BTC.R**: main R script that loads global_data.RDS
* **`functions/`**: R scripts for data processing, visualization, and statistical testing
* **`Config/`**: Configuration files including metadata and panel

## Data Access
Due to GitHub's file size limits, global_data.RDS is hosted externally:

* **Processed RDS File:** https://zenodo.org/records/18818000  
    *This file is required to run the scripts provided in this repository.*
