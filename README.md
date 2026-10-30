# Snippy example: clonal structure within host 22M21

Example [Snippy](https://github.com/tseemann/snippy) (v4.6.0) workflow for within-host clonal structure of *Enterococcus faecalis* from host **22M21** (11 isolates).

Pipeline: RefSeq reference → static mask (repeats / phage / IS) → Snippy on QC-ed FASTQ → keep ALIGNED/LENGTH ≥ 0.8 → masked core SNPs → FastTree + snp-dists → heatmap & tree. Optional Gubbins is commented in §8.

## Contents

| Path | Description |
|------|-------------|
| `run_pipeline.sh` | End-to-end script (§0–7 main; §8 optional Gubbins) |
| `metadata.txt` | Sample table (`Id`, `Species`, `Completeness`, …) |
| `script/plot_heatmap.R` | Pairwise SNP distance heatmap |
| `script/plot_tree.R` | Core SNP tree (tip color by cluster; ★ = support ≥ 0.7) |
| `ref/refseq_reference_genomes.tsv` | RefSeq reference list (Species → download URL) |
| `result/mask_static.bed` | Example static mask |
| `result/core.tree` / `core.snpdists.tsv` | Example core tree & pairwise distances |
| `result/core_gubbins.tree` / `core_gubbins.snpdists.tsv` | Example after optional Gubbins |
| `figure/` | Example plots: `heatmap` / `tree` (pdf + png) |

Raw reads and large intermediates (`*.full.aln`, etc.) are **not** included. Contact the author for reads.

## Requirements

- conda / mamba
- Env `snippy` (versions in `run_pipeline.sh` §0):  
  `snippy=4.6.0` `fastp=1.4.0` `snp-dists=0.8.2` `fasttree=2.2.0` `minimap2=2.28` `bedtools=2.31.1` `samtools=1.20` `phastaf=0.1.0` `isescan=1.7.3`
- Optional env `gubbins`: `gubbins=3.4.3` `python=3.10` (separate; no bedtools)
- Plots: R with `ggplot2`, `scales`; tree also needs `ggtree`, `treeio`, `ape`  
  (conda: `r-ape r-patchwork bioconductor-ggtree bioconductor-treeio`)

```bash
conda create -y -n snippy -c bioconda -c conda-forge \
  snippy=4.6.0 fastp=1.4.0 snp-dists=0.8.2 fasttree=2.2.0 \
  minimap2=2.28 bedtools=2.31.1 samtools=1.20 phastaf=0.1.0 isescan=1.7.3
conda create -y -n gubbins -c bioconda -c conda-forge gubbins=3.4.3 python=3.10
conda activate snippy
```

## Inputs

Run `bash run_pipeline.sh` from the repo root. The script does `mkdir -p snippy && cd snippy`, so prepare:

```text
# repo root
run_pipeline.sh
script/plot_heatmap.R
script/plot_tree.R

# after cd snippy/
snippy/
  metadata.txt
  ref/refseq_reference_genomes.tsv
  seq/{id}_1.fq.gz
  seq/{id}_2.fq.gz
```

`{id}` = `Id` in `metadata.txt`. Prefer QC-ed paired-end FASTQ (`--R1/--R2`). `snippy --ctgs` is shown commented in §4 (pseudo-reads; no real depth/quality).

## Quick start

```bash
conda activate snippy
bash run_pipeline.sh

# re-plot from snippy/ workdir (paths as written by the pipeline):
# Rscript ../script/plot_heatmap.R result/core/core.snpdists.tsv figure
# Rscript ../script/plot_tree.R result/core/core.tree result/core/core.snpdists.tsv figure
```

Steps in the script:

1. **QC** — fastp  
2. **Reference** — Species → URL in `ref/refseq_reference_genomes.tsv` → `ref/ref.fna`  
3. **Static mask** — repeats + phage + IS → `temp/mask/mask_static.bed`  
4. **Snippy** — per-sample calling  
5. **Filter** — ALIGNED/LENGTH ≥ 0.8 → `result/core.highAlign.txt`  
6. **Core SNPs** — `snippy-core --mask` → `result/core/`  
7. **Tree / distances / plots** — FastTree + snp-dists → heatmap & tree in `figure/`  
8. **Optional** — Gubbins (§8; plots go to `figure/gubbins`)

## Outputs

When you run the pipeline, main products land under `snippy/` as:

| Path | Description |
|------|-------------|
| `temp/fastp/` | QC FASTQ + reports |
| `temp/mask/mask_static.bed` | Merged static mask |
| `temp/core_all/` | Unfiltered snippy-core (for the ALIGNED filter) |
| `result/snippy/{id}/` | Per-sample Snippy dirs |
| `result/snippy_pass/` | Symlinks to samples that pass the filter |
| `result/core.highAlign.txt` | Pass list |
| `result/core/core.aln` / `core.full.aln` | Masked core alignments |
| `result/core/core.tree` | FastTree Newick |
| `result/core/core.snpdists.tsv` | Pairwise SNP distances |
| `result/gubbins/` | Optional Gubbins outputs (§8) |
| `figure/heatmap.pdf` / `.png` | Distance heatmap |
| `figure/tree.pdf` / `.png` | Phylogenetic tree |

This repo ships a **small example slice** under `result/` and `figure/` (tree, distances, mask, plots) for inspection without re-running.

## Notes

- Mask BED contig names must match the reference (see commented check in §6).
- Soft links under `result/snippy_pass/` use absolute paths; moving the project breaks them.
- Threshold `0.8` is set in §5; or edit `result/core.highAlign.txt` and re-run masked `snippy-core`.
- Cluster cut for tip colors matches the heatmap (pairwise SNP ≤ 100).
- Large regenerable files (`*.full.aln`) are omitted from the example package.
- §7 (static mask) is the default endpoint; uncomment §8 only if you need recombination-aware distances.

## Citation

- Snippy: Seemann T. https://github.com/tseemann/snippy  
- This example: Zhu Zhihao (contact for reads and questions)

## License

Pipeline scripts and example outputs in this repository: MIT. Snippy itself is GPL-2.0.
