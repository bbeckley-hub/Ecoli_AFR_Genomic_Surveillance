
# Ecoli_AFR_Genomic_Surveillance

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue.svg)](https://www.python.org/)
[![R 4.3+](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)

Comparative genomic analysis of *Escherichia coli* sequence types **ST10**,
**ST131**, and **ST167** across 14 African countries and Jordan.

---

## Overview

This repository contains the analysis code, intermediate data tables, and
supplementary data supporting the manuscript:

> **Beckley B., et al.** *Comparative genomic analysis of Escherichia coli ST10,
> ST131, and ST167 in Africa and Jordan reveals distinct resistance,
> virulence, and surface-antigen architectures.* (Manuscript under review.)

**Dataset:** 444 high-quality *E. coli* genomes
(ST167 = 189, ST131 = 151, ST10 = 104) selected from 1,816 initial assemblies
across 28 NCBI BioProjects.

**Analysis layers:**

| Layer | Tool(s) | Script(s) |
|---|---|---|
| Species confirmation | FastANI | EcoliTyper wrapper |
| Assembly QC | BUSCO 6.1.0 | `busco_command.txt`, `generate_busco_summary.sh` |
| MLST / serotype / CH type / phylogroup | EcoliTyper 1.3.0 | `ecolityper_command.txt` |
| AMR genes | AMRFinderPlus, ResFinder, CARD, MEGARes, ARG-ANNOT (via ABRicate) | `amr_analysis.py`, `amr_burden.py` |
| Chromosomal mutations | AMRFinderPlus | `mutation_analysis.py`, `mutation_burden.py` |
| Virulence genes | VFDB (via ABRicate) | `virulence_analysis.py` |
| Plasmid replicons | PlasmidFinder (via ABRicate) | `plasmid_analysis.py` |
| Capsule/O-antigen loci | Kaptive 3.3.2 | `capsule_summary.py`, `surface_antigen_summary.py` |
| Annotation | Prokka 1.14.6 | `run_prokka.sh` |
| Pangenome | Roary 3.13.0 | `run_roary.sh`, `pangenome_analysis.py`, `pangenome_summary_plot.py` |
| Core phylogeny | IQ-TREE 3.0.1, ClonalFrameML 1.2.0, Gubbins 3.4.3 | Phylogeny.txt & Recombination.sh |
| Statistics | R 4.3.3, Python 3.12.2 | `statistical_analysis.py`, `run_final_statistics.py` |

---

## Repository structure

```
Ecoli_AFR_Genomic_Surveillance/
├── LICENSE
├── README.md
├── .gitignore
├── example_data/                  # Small example outputs (format reference)
├── Intermediate_data_tables/      # Contingency tables, associations, stats
├── metadata/                      # Sample metadata & ST-by-country summary
├── Scripts/                       # All analysis scripts (see below)
└── Supplementary_Data/            # Supplementary Files 1–6
```

---

## Requirements

### System tools

| Tool | Version used | Purpose |
|---|---|---|
| [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | v18.29.1 | Genome download |
| [FastANI](https://github.com/ParBLiSS/FastANI) | 1.34 | Species confirmation |
| [BUSCO](https://busco.ezlab.org/) | 6.1.0 | Assembly completeness |
| [EcoliTyper](https://github.com/bbeckley-hub/EcoliTyper) | 1.3.0 | Pipeline for MLST, serotyping, AMR, etc. |
| [Prokka](https://github.com/tseemann/prokka) | 1.14.6 | Genome annotation |
| [Roary](https://sanger-pathogens.github.io/Roary/) | 3.13.0 | Pangenome |
| [IQ-TREE](http://www.iqtree.org/) | 3.0.1 | Phylogeny |
| [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML) | 1.2.0 | Recombination |
| [Gubbins](https://github.com/nickjcroucher/gubbins) | 3.4.3 | Recombination |
| [Kaptive](https://github.com/klebgenomics/Kaptive) | 3.3.2 | Capsule/O-locus typing |
| [EC-K-typing](https://github.com/rgladstone/EC-K-typing) | 4.0.0 | K-locus typing |


### Python packages

Install via `requirements.txt` (see below). Core:

- `pandas`, `numpy`, `scipy`, `statsmodels`
- `matplotlib` (for `pangenome_summary_plot.py`, `st_enriched_genes.py`)

---

## Installation

### 1. Clone the repo

```bash
git clone git@github.com:bbeckley-hub/Ecoli_AFR_Genomic_Surveillance.git
cd Ecoli_AFR_Genomic_Surveillance
```

### 2. Set up Python environment

```bash
Conda install or pip install [package]
```

### 3. Install external tools

The simplest way is via conda/mamba (see `environment.yml`):

```bash
conda env create -f environment.yml
conda activate ecoli_afr
```

### 4. Download the BUSCO lineage dataset (once)

`busco_command.txt` uses `--offline`, so you must fetch the lineage first:

```bash
busco --download enterobacterales_odb12
```

---

## Reproducing the analysis

> **Working directory:** All commands below assume you are inside `Scripts/`
> unless stated otherwise. The Python scripts read inputs from and write
> outputs to the current working directory.

```bash
cd Scripts
```

### Step 1 — Download genomes

Reads `ecoli_clean.tsv` (columns: `acc`, `biosample`, `taxonomy`,
`bioproject`, `country`) and writes one FASTA per accession to `genomes/`.

```bash
bash download_genomes.sh
```

- Retries each accession up to 3 times (10 s delay).
- Failed accessions are written to `failed_accessions.txt`.
- Full log: `download.log`.

### Step 2 — Assembly quality control (BUSCO)

```bash
# Paste the loop from busco_command.txt, or run:
bash busco_command.txt
```

Outputs `busco_results/<accession>_busco/short_summary.*.txt`. Runs two
genomes in parallel at 4 CPUs each.

### Step 3 — Summarise BUSCO results

⚠️ **Edit first:** `generate_busco_summary.sh` contains a hardcoded path
(`cd ~/Africa/genomes/Sts/busco_results/`). Change it to your `busco_results/`
directory before running.

```bash
bash generate_busco_summary.sh
```

Outputs `summary_busco.tsv`.

### Step 4 — Genotyping (EcoliTyper)

```bash
ecolityper -i "*.fna" -o results \
    --amr-min-identity 0.95 --amr-min-coverage 0.9 \
    --abricate-minid 90 --abricate-mincov 85
```

Produces per-genome `amr_genes.csv`, `virulence_genes.csv`,
`plasmid_markers.csv`, `mutations.csv`, `sample_overview.csv`, and an
interactive HTML report (see `example_data/` for the format).

### Step 5 — Annotate genomes (Prokka)

```bash
bash run_prokka.sh
```

Writes `prokka_annotations/<accession>/` (`.gff`, `.gbk`, `.faa`, …).
Summary in `prokka_summary.tsv`. Uses 14 CPUs.

### Step 6 — Build pangenome (Roary)

```bash
bash run_roary.sh
```

Equivalent to:

```bash
roary -e --mafft -p 14 -i 90 -cd 99 -s \
    -f roary_output prokka_annotations/*/*.gff
```

Outputs `roary_output/gene_presence_absence.csv` and `.Rtab`.

### Step 7 — Core-genome phylogeny

```bash
iqtree3 -s roary_output/core_gene_alignment.aln \
    -m MFP -B 1000 -T AUTO --prefix core_tree
clonalframeml -s roary_output/core_gene_alignment.aln \
    -em -if core_tree.treefile
run_gubbins.py -p gubbins_out roary_output/core_gene_alignment.aln
```

(Filter recombinant sites with Gubbins before the final IQ-TREE run.)

### Step 8 — Downstream analyses

Each script reads from and writes to the current directory.

#### 8a. AMR

```bash
python amr_analysis.py          # → amr_stat_results.csv
python amr_burden.py            # → stdout (median/IQR + KW + Dunn)
```

#### 8b. Chromosomal resistance mutations

```bash
python mutation_analysis.py     # → stdout (Fisher + FDR)
python mutation_burden.py       # → stdout
```

#### 8c. Virulence

```bash
python virulence_analysis.py    # → virulence_associations.csv
```

#### 8d. Plasmids

```bash
python plasmid_analysis.py      # → plasmid_associations.csv
```

#### 8e. Pangenome enrichment & ordination

```bash
python pangenome_analysis.py    # → st_enriched_genes_fast.csv
python pangenome_summary_plot.py  # → pangenome_summary_bar.{png,pdf}
```

`st_enriched_genes.py` plots MDS coordinates from
`pangenome_mds_coordinates.csv`:

```bash
python st_enriched_genes.py     # → pangenome_mds_plot.png
```

#### 8f. Capsule / surface-antigen loci (Kaptive outputs)

```bash
python capsule_summary.py          # → stdout (E. coli Group 2/3 only)
python surface_antigen_summary.py  # → stdout (E. coli + Klebsiella K + O)
```

#### 8g. Population structure & chi-square tests

Two related scripts:

```bash
python compute_geo.py             # → stdout (crosstabs, assembly metrics)
python statistical_analysis.py    # → ecoli_comprehensive_harmonized.tsv
python run_final_statistics.py    # → statistical_summary.tsv + per-test tables
```

`run_final_statistics.py` produces, for each contingency table
(`st_country`, `st_phylogroup`, `st_serotype`, `st_chtype`, `st_otype`):

- `<key>_observed.tsv`
- `<key>_expected.tsv`
- `<key>_residuals.tsv`

plus a combined `statistical_summary.tsv`.

---

## Data dictionary

### `metadata/`

| File | Description |
|---|---|
| `ecoli_comprehensive.csv` | **Tab-separated** despite `.csv` extension. Columns: `Assembly`, `Country`, `MLST`, `Phylogroup`, `Serotype`, `CH_Type`, `O_Type`, `N50`, `GC_Content`, `Total_Length`, … |
| `st_by_country.csv` | Long-format summary of ST counts per country |

### `Intermediate_data_tables/`

| File | Produced by | Description |
|---|---|---|
| `amr_final_results.csv` | EcoliTyper aggregation | Merged AMR gene presence across all databases |
| `amr_gene_summary.csv` | EcoliTyper aggregation | Canonical gene × genome presence |
| `amr_raw_gene_summary.csv` | EcoliTyper aggregation | Pre-canonicalisation |
| `amr_stat_results.csv` | `amr_analysis.py` | Fisher's exact + BH-FDR per gene × ST |
| `plasmid_associations.csv` | `plasmid_analysis.py` | Fisher's exact + BH-FDR per replicon × ST |
| `virulence_associations.csv` | `virulence_analysis.py` | Fisher's exact + BH-FDR per gene × ST |
| `st_country_{observed,expected,residuals}.tsv` | `run_final_statistics.py` | Contingency table + chi-square components |
| `st_chtype_*.tsv` | `run_final_statistics.py` | Same, for CH types |
| `st_otype_*.tsv` | `run_final_statistics.py` | Same, for O types |
| `st_phylogroup_*.tsv` | `run_final_statistics.py` | Same, for phylogroups |
| `st_serotype_*.tsv` | `run_final_statistics.py` | Same, for serotypes |

### `example_data/`

| File | Description |
|---|---|
| `sample_overview.csv` | Per-genome ST / country / phylogroup summary |
| `amr_genes.csv` | EcoliTyper AMR output format |
| `virulence_genes.csv` | EcoliTyper virulence output format |
| `plasmid_markers.csv` | EcoliTyper plasmid output format |
| `mutations.csv` | EcoliTyper chromosomal-mutation output |
| `bacmet_genes.csv` | BacMet biocidal-metal resistance hits |
| `high_risk_combinations.csv` | Predefined high-risk AMR combinations |
| `genius_ecoli_ultimate_gene_centric_report.html` | Full EcoliTyper HTML report |

### `Supplementary_Data/`

See manuscript for full descriptions.

- **Supplementary File 1** — BioProjects & accession lists (28 files + search strategy)
- **Supplementary File 2** — QC metrics (`FASTA_QC_summary.tsv/html`, `summary_busco.tsv`)
- **Supplementary File 3** — Curated 57-gene virulence list (XLSX)
- **Supplementary File 4** — Interactive EcoliTyper HTML report
- **Supplementary File 5** — Kaptive output (E. coli Group 2/3, Klebsiella K & O)
- **Supplementary File 6** — EcoliTyper results, IQ-TREE tree, Roary matrices
- **Supplementary Figure S1** — Pangenome composition bar chart (PNG/PDF)

---

## Outputs at a glance

| Figure / Table | Source | Script |
|---|---|---|
| Fig. 1 — Geographic distribution map | MapChart | (manual) |
| Fig. 2 — PCoA of accessory genome | `pangenome_*` | `pangenome_summary_plot.py` / Python ordination |
| Table 1 — ST × country | `st_country_*` | `run_final_statistics.py` |
| Table 2 — Phylogroup / serotype / CH summary | `metadata/` | `statistical_analysis.py` |
| Table 3 — AMR genes | `amr_stat_results.csv` | `amr_analysis.py` |
| Table 4 — Virulence genes | `virulence_associations.csv` | `virulence_analysis.py` |
| Table 6 — Chromosomal mutations | stdout | `mutation_analysis.py` |
| Supp. Fig. S1 — Pangenome bar | `pangenome_summary_bar.*` | `pangenome_summary_plot.py` |

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{Beckley2026EcoliAFR,
  author  = {Beckley, Brown},
  title   = {Comparative genomic analysis of {Escherichia} coli {ST10},
             {ST131}, and {ST167} in {Africa} and {Jordan} reveals distinct
             resistance, virulence, and surface-antigen architectures},
  journal = {[Manuscript under review]},
  year    = {2026}
}
```

A machine-readable `CITATION.cff` is included.

---

## License

Code is released under the [MIT License](LICENSE).

Derived data tables in `Intermediate_data_tables/`, `metadata/`, and
`Supplementary_Data/` are provided for convenience under the same terms.
Underlying genome assemblies are publicly available from NCBI/ENA and
remain subject to their respective terms of use.

---

## Contact

**Brown Beckley**
Department of Medical Biochemistry, University of Ghana Medical School, Accra
Department of Biochemistry and Biotechnology, KNUST, Kumasi
brownbeckley94@gmail.com
[@bbeckley-hub](https://github.com/bbeckley-hub)
