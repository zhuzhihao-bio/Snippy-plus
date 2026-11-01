# Snippy example: clonal structure within host 22M21

Example [Snippy](https://github.com/tseemann/snippy) (v4.6.0) workflow for within-host clonal structure of *Enterococcus faecalis* from host **22M21** (11 isolates).

Pipeline: optional GSA FASTQ download → QC → RefSeq reference → static mask (repeats / phage / IS) → Snippy → keep ALIGNED/LENGTH ≥ 0.8 → masked core SNPs → FastTree + snp-dists → heatmap & tree. Optional Gubbins is commented at the end.

## Contents

| Path | Description |
|------|-------------|
| `run_pipeline.sh` | End-to-end script (download commented; §0–7 main; Gubbins optional) |
| `metadata.txt` | Sample table (`Id`, `Species`, `Bioproject`, `Run`=CRR*, …) |
| `script/plot_heatmap.R` | Pairwise SNP distance heatmap |
| `script/plot_tree.R` | Core SNP tree (tip color by cluster; ★ = support ≥ 0.7) |
| `ref/refseq_reference_genomes.tsv` | RefSeq reference list (Species → download URL) |
| `result/mask_static.bed` | Example static mask |
| `result/core.tree` / `core.snpdists.tsv` | Example core tree & pairwise distances |
| `result/core_gubbins.tree` / `core_gubbins.snpdists.tsv` | Example after optional Gubbins |
| `figure/` | Example plots: `heatmap` / `tree` (pdf + png) |

Raw reads are **not** shipped. Use the commented GSA download in `run_pipeline.sh`, or place your own FASTQs.

## Requirements

- conda / mamba
- Env `snippy` (see `run_pipeline.sh` §0):  
  `snippy=4.6.0` `fastp=1.4.0` `snp-dists=0.8.2` `fasttree=2.2.0` `minimap2=2.28` `bedtools=2.31.1` `samtools=1.20` `phastaf=0.1.0` `isescan=1.7.3`
- Optional download: `iseq` (`conda install -y -c bioconda -c conda-forge iseq`; Aspera optional for `-a`)
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
  metadata.txt          # Run column = GSA CRR accessions
  ref/refseq_reference_genomes.tsv
  seq/{id}_1.fq.gz      # from download step, or your own
  seq/{id}_2.fq.gz
```

`{id}` = `Id` in `metadata.txt`. Prefer paired-end FASTQ (`--R1/--R2`). `snippy --ctgs` is shown commented in §4 (pseudo-reads; no real depth/quality).

## Quick start

```bash
conda activate snippy
# optional: uncomment the GSA download block in run_pipeline.sh, then:
#   conda install -y -c bioconda -c conda-forge iseq
bash run_pipeline.sh

# re-plot from snippy/ workdir:
# Rscript ../script/plot_heatmap.R result/core/core.snpdists.tsv figure
# Rscript ../script/plot_tree.R result/core/core.tree result/core/core.snpdists.tsv figure
```

Steps in the script:

1. **Download (commented)** — iSeq by `Run` (CRR*) → `seq/{id}_1/2.fq.gz`  
2. **QC** — fastp  
3. **Reference** — Species → URL in `ref/refseq_reference_genomes.tsv` → `ref/ref.fna`  
4. **Static mask** — repeats + phage + IS → `temp/mask/mask_static.bed`  
5. **Snippy** — per-sample calling  
6. **Filter** — ALIGNED/LENGTH ≥ 0.8 → `result/core.highAlign.txt`  
7. **Core SNPs** — `snippy-core --mask` → `result/core/`  
8. **Tree / distances / plots** — FastTree + snp-dists → `figure/`  
9. **Optional** — Gubbins (commented; plots → `figure/gubbins`)

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
| `result/gubbins/` | Optional Gubbins outputs |
| `figure/heatmap.pdf` / `.png` | Distance heatmap |
| `figure/tree.pdf` / `.png` | Phylogenetic tree |

This repo ships a **small example slice** under `result/` and `figure/` (tree, distances, mask, plots) for inspection without re-running. Large regenerable files (`*.full.aln`) are omitted.

## Notes

- Mask BED contig names must match the reference (see commented check before masked `snippy-core`).
- Soft links under `result/snippy_pass/` use absolute paths; moving the project breaks them.
- Threshold `0.8` is set in the filter step; or edit `result/core.highAlign.txt` and re-run masked `snippy-core`.
- Cluster cut for tip colors matches the heatmap (pairwise SNP ≤ 100).
- Uncomment Gubbins only if you need recombination-aware distances.

## Citation

- Snippy: Seemann T. https://github.com/tseemann/snippy  
- This example: Zhu Zhihao (contact for reads and questions)

## License

Pipeline scripts and example outputs in this repository: MIT. Snippy itself is GPL-2.0.
