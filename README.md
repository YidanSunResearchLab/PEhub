# PEhub

**PEhub** is a promoter-centric computational framework for identifying and characterizing **multi-way enhancer hubs** from HiChIP chromatin interaction data.
It resolves higher-order enhancer cooperation as **promoter-anchored regulatory units**, enabling quantitative analysis of hub architecture, stability, and statistical significance.

---

## Key Features

* **Promoter-centric modeling** of enhancer cooperation
* Detection of **multi-way enhancer hubs** beyond pairwise loops
* Flexible **weighting schemes** (distance-, signal-, and bin-aware)
* **Leiden clustering** on enhancer co-membership graphs
* **Bootstrap-based stability** assessment
* **Distance-aware null models** for global hub p-values
* Fully **CLI-driven**, reproducible end-to-end workflow
* Designed for **HiChIP / H3K27ac HiChIP / Hi-C–derived loops**

---

## Conceptual Overview

PEhub reframes enhancer hub detection as a **local, promoter-conditioned problem**.

Instead of clustering a global chromatin interaction graph, PEhub:

1. Anchors all interactions at individual promoters
2. Builds **weighted enhancer–enhancer co-membership matrices**
3. Partitions enhancers into hubs using **Leiden community detection**
4. Quantifies hub strength, density, stability, and statistical significance

> 📌 **Figure**:
> A schematic overview of the PEhub workflow is shown below.

```
docs/PEhub_overview.png
```

(Replace this path with your actual figure location.)

---

## Installation

PEhub is distributed as a single R script with a fully reproducible **Conda environment** that includes all required R, Bioconductor, and Python dependencies.

### 1. Clone the repository

```bash
git clone https://github.com/YidanSunResearchLab/PEhub
cd PEhub
```

---

### 2. Create Conda environment

PEhub provides a pre-defined Conda environment that installs:

* R (≥ 4.3)
* Required CRAN and Bioconductor packages
* Python + Leiden clustering backend

Create and activate the environment:

```bash
conda env create -f env.yml
conda activate pehub
```

This step only needs to be done once.

---

### 3. Verify installation

Check that R and Python are correctly configured:

```bash
Rscript PEhub.R --help
```

If no errors are shown, the environment is ready.

---

### 4. Python path for Leiden

When running PEhub, pass the Python interpreter from the active Conda environment:

```bash
which python
```

Example output:

```
/path/to/miniconda/envs/pehub/bin/python
```

Use this path with the `--python_path` argument.

---



## Input Data

### Required inputs

| Argument        | Description                                            |
| --------------- | ------------------------------------------------------ |
| `--tss_gtf`     | GTF file used to extract transcript TSS / promoters    |
| `--loop_sig`    | Significant interactions BEDPE (e.g. HicDCPlus output) |
| `--loop_all`    | All interactions BEDPE(.gz), used for null modeling    |
| `--outdir`      | Output directory                                       |
| `--sample`      | Sample identifier used in output file naming           |
| `--python_path` | Python interpreter path (for Leiden clustering)        |

---

### Optional but recommended parameters

| Argument            | Description                                                             |
| ------------------- | ----------------------------------------------------------------------- |
| `--method`          | Co-membership normalization (`log_minmax`, `log_zscore`, `log_maxnorm`) |
| `--weight_method`   | Edge weighting scheme (e.g. `bin_log_ratio_sig`)                        |
| `--promoter_window` | Promoter window around TSS (bp)                                         |
| `--k_min`           | Minimum enhancers per promoter                                          |
| `--resolution`      | Leiden resolution parameter                                             |
| `--quantile_cutoff` | Sparsification cutoff for co-membership                                 |
| `--B_stability`     | Bootstrap iterations                                                    |
| `--B_pvalue`        | Monte Carlo iterations for null p-values                                |
| `--null_mode`       | Null sampling strategy                                                  |
| `--workers`         | Parallel workers                                                        |

---

## Usage

### Basic example

```bash
Rscript PEhub.R \
  --tss_gtf example_data/genes.chr9.gtf.gz \
  --loop_sig example_data/hicdcplus.significant_interactions.bedpe \
  --loop_all example_data/hicdcplus.all_interactions.bedpe.gz \
  --outdir /path/to/pehub_test \
  --sample GM12878test \
  --python_path $(which python) \
  --method log_minmax \
  --weight_method bin_log_ratio_sig \
  --promoter_window 0 \
  --k_min 3 \
  --resolution 1.0 \
  --quantile_cutoff 0.2 \
  --pvalue_cutoff 0.05 \
  --stability_cutoff 0.5 \
  --B_stability 10 \
  --B_pvalue 1000 \
  --null_mode hist_matched \
  --workers 10
```

---

## Output Files

All outputs are written under `--outdir`.

### Intermediate

1. **Preprocessing**

```
multiple_result.<sample>.hub.all.preprocess.RData
```

2. **Observed hub discovery**

```
multiple_result.<sample>.hub.<method>.<weight_method>.preprocess.RData
```

---

### Final results

3. **Final hub statistics**

```
multiple_result.<sample>.hub.<method>.<weight_method>.RData
```

4. **Significant enhancer hubs**

* BED:

  ```
  multiple_result.<sample>.hub.<method>.<weight_method>.bed
  ```
* TXT (full annotation):

  ```
  multiple_result.<sample>.hub.<method>.<weight_method>.txt
  ```
* BEDPE (hub-associated interactions):

  ```
  multiple_result.<sample>.hub.<method>.<weight_method>.bedpe
  ```

5. **Filtered_out interactions**

```
multiple_result.<sample>.pairwise.<method>.<weight_method>.bed
multiple_result.<sample>.pairwise.<method>.<weight_method>.bedpe
```

---

## Interpreting the Results

Each enhancer hub is characterized by:

* **Hub size** (number of enhancers)
* **Weighted density / internal score**
* **Bootstrap stability metrics**

  * Jaccard stability score
  * Reproducibility rate
  * Existence rate
* **Global p-value and FDR**

  * Distance-aware null models
  * Distance-matched resampling

Significant hubs represent **cooperative regulatory units** rather than isolated enhancer–promoter loops.

---

## Reproducibility

* All random steps use fixed seeds
* Parameters are fully CLI-controlled
* Intermediate RData files allow stepwise re-analysis
* Designed for HPC and parallel execution

---

## Citation

If you use PEhub in your work, please cite:

> **Tan J and Sun Y**
> *PEhub resolves the hierarchical regulatory architecture of multi-way enhancer hubs in the human brain.*
> *(Manuscript in submission)*


---

## License

This project is released under the **MIT License**.
See the `LICENSE` file for details.

---

