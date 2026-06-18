# ORA Manuscript Draft

This directory contains the LaTeX manuscript draft and verified BibTeX references.

## Files

- `main.tex`: full manuscript draft with main sections, methods, limitations, figure callouts, and bibliography.
- `references.bib`: BibTeX entries with DOI-backed references.

## Build

The draft references manuscript figures from `../results/figures/manuscript_figure*.pdf`.

```bash
PYTHON=.venv/bin/python make manuscript-figures
make manuscript
```

The current local machine does not have `latexmk`, `pdflatex`, or `bibtex` on `PATH`, so `make manuscript` will report the missing TeX engine until one is installed.
