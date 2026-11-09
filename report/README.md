# Report Guide

Technical report for the Falcon 9 simulation, written in LaTeX with the TUM `tumarticle` class.

## Structure

- `main.tex` — master document: title, authors, packages, and the ordered list of chapter includes.
- `content/` — one file per chapter (e.g. `introduction.tex`), each starting with `\section{...}`.
- `content/acronyms.tex` — all acronyms, declared with `\DeclareAcronym{KEY}{short = ..., long = ...}`; use `\ac{KEY}` in text.
- `content/sample_content.tex` — reference for table, figure, and heading syntax.
- `images/` — all figures (PDF preferred), referenced as `images/<file>.pdf`.
- `tum/` — TUM LaTeX classes and styles; found automatically via `latexmkrc`.

## Adding a chapter

1. Create `content/<name>.tex`, starting with `\section{Chapter Title}`.
2. Register it in `main.tex` below the `% +++++ Main part` marker:

   ```latex
   \input{content/introduction}
   %\input{content/next_section}
   ```

   The order of `\input` lines defines the chapter order in the PDF.

## Generating the PDF

From this folder (`report/`), run:

```bash
latexmk -pdf main.tex
```

The result is `main.pdf`. `latexmk` automatically runs `pdflatex` and `bibtex` as many times as needed (for references and the table of contents), so a single call is enough.

Clean up the auxiliary files afterwards with:

```bash
latexmk -c
```

If you use the VS Code task: `` Ctrl+Shift+B `` runs the build from anywhere in the repo.

## Conventions

- Headings: `\section` / `\subsection` / `\subsubsection` — no manual numbering.
- Figures: `\begin{figure}[htbp]` with `\centering`, `\includegraphics`, `\caption`, and `\label{fig:<name>}`; reference with `Figure~\ref{fig:<name>}`.
- Tables: `booktabs` style (`\toprule`, `\midrule`, `\bottomrule`), see `content/sample_content.tex`.
- Math: inline `$...$`; units in math mode with a thin space, e.g. `$8\,\mathrm{MHz}$`.
