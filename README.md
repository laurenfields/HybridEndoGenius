# HybridEndoGenius

HybridEndoGenius is a **hybrid neuropeptide identification workflow** designed to improve the discovery of endogenous peptides from mass spectrometry (MS/MS) data. It integrates **database-driven** and **de novo–assisted** strategies to enable sensitive and flexible neuropeptide identification across different species and experimental settings.

The workflow is designed to be **modular, reproducible, and scalable**, and is compatible with high-throughput computing environments such as **CHTC/HTCondor**.

---

## Key Features

- Hybrid peptide identification strategy combining:
  - Database search–based identification
  - De novo sequencing–assisted discovery
- Optimized for **neuropeptides**, which are typically:
  - Short
  - Low abundance
  - Poorly annotated in standard protein databases
- Flexible support for **species-specific configurations**
- Designed for **HTC execution** (HTCondor / DAGMan)

---

## Workflow Overview

At a high level, HybridEndoGenius performs the following steps:

1. **Input preparation**
   - MS/MS data (e.g. `.mzML`, `.ms2`, `.mgf`)
   - FASTA and CSV databases
2. **Database search**
   - Identification of known neuropeptides
3. **De novo sequencing**
   - Discovery of novel peptide candidates
4. **Filtering and post-processing**
   - Removal of peptides with low confidence
   - Identification of putative novel neuropeptides
5. **Result integration**
   - Database search using the FASTA generated from de novo sequencing

---
