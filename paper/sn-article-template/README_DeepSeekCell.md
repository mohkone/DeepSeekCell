# DeepSeekCell BMC Bioinformatics Manuscript Files

Use these adapted files for the BMC Bioinformatics submission draft:

- `deepseekcell_bmc_bioinformatics.tex`: main manuscript in Springer Nature/BMC `sn-article` format.
- `deepseekcell_bmc_refs.bib`: manuscript-specific bibliography.

The original template files (`sn-article.tex`, `sn-article.pdf`, `sn-jnl.cls`, and `bst/`) are unchanged.

## Compile Command

From this folder, run:

```powershell
pdflatex deepseekcell_bmc_bioinformatics.tex
bibtex deepseekcell_bmc_bioinformatics
pdflatex deepseekcell_bmc_bioinformatics.tex
pdflatex deepseekcell_bmc_bioinformatics.tex
```

The manuscript uses:

```tex
\documentclass[referee,lineno,pdflatex,sn-vancouver-num]{sn-jnl}
```

This gives double-line spacing, line numbering, and Vancouver-style numbered references suitable for BMC-style review drafts.

## Before Submission

Current manuscript status:

- Author names, funding statement, competing interests statement, figures, tables, and captions have been inserted in the BMC template manuscript.
- The GitHub-Zenodo archive DOI has been inserted in the Code availability and Data availability text: https://doi.org/10.5281/zenodo.20680434.
- Add ORCID IDs if the final author list requires them.

## Figure Status

Embedded figure files:

- `deepseekcell_workflow.pdf`
- `deepseekcell_shiny_interface.png`
- `benchmark_macroF1.pdf`
- `benchmark_accuracy.pdf`
- `benchmark_clade_accuracy.pdf`
- `benchmark_runtime.pdf`
- `global_ARI.pdf`
- `global_MacroF1.pdf`
- `global_Accuracy.pdf`
- `global_Runtime.pdf`
