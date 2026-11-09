# Report Guide (for AI agents)

This report is written in LaTeX using the TUM `tumarticle` class. It is compiled from `report/main.tex` into `report/main.pdf`.

## Structure

- `report/main.tex` — master document. Metadata (title, subtitle, authors), package loading, and the ordered list of chapter includes.
- `report/content/` — one file per chapter (e.g. `introduction.tex`). Each file starts with `\section{...}` and holds the full chapter including its subsections.
- `report/images/` — all figures (PDF preferred). Referenced from chapter files as `images/<file>.pdf`.
- `report/content/acronyms.tex` — all acronyms, declared with `\DeclareAcronym{KEY}{short = ..., long = ...}`.
- `report/content/sample_content.tex` — reference for tables, figures, and heading syntax.

## How to add a chapter

1. Create the chapter file: `report/content/<name>.tex`, starting with `\section{Chapter Title}`.
2. Register it in `report/main.tex` below the `% +++++ Main part` marker:

   ```latex
   \input{content/introduction}
   %\input{content/next_section}
   ```

   The order of `\input` lines defines the order of chapters in the PDF. Use the `.tex` extension implicitly — `\input{content/<name>}` resolves to `content/<name>.tex`.

3. Rebuild with `latexmk` (a `latexmkrc` exists in `report/`).

## Conventions

- Headings: `\section` / `\subsection` / `\subsubsection` — no manual numbering.
- Acronyms: use `\ac{KEY}` in text; define new keys only in `content/acronyms.tex`.
- Figures (place in `report/images/`, then include in the chapter file):

  ```latex
  \begin{figure}[htbp]
    \centering
    \includegraphics[width=\linewidth]{images/<name>.pdf}
    \caption{Caption below the figure.}
    \label{fig:<name>}
  \end{figure}
  ```

  Always set both `caption` and `label`; reference with `Figure~\ref{fig:<name>}`.
- Tables: use the `booktabs` style from `sample_content.tex` (`\toprule`, `\midrule`, `\bottomrule`).
- Math: inline with `$...$`; units in math mode with a thin space, e.g. `$8\,\mathrm{MHz}$`.

## Writing style

The report documents the simulator (the whole project). It is a technical report, not an academic essay: describe the system, its models, and its implementation.

- Professional, precise, direct English. Present tense for how the simulator works; the system is the subject, not the authors.
- Be concrete: name the actual component, model, or parameter being described, give values with units, and show the governing equations.
- State assumptions, model limitations, and simplifications explicitly where the affected model is described.
- Organize each chapter around one part of the simulator (e.g. environment, aerodynamics, propulsion); keep motivational text brief.
- Support every design decision with an equation, figure, table, or reference; keep notation and units consistent across chapters.
- Report validation against reference data where available (e.g. published Falcon 9 / Falcon Heavy performance figures).
- No filler or emphasis words ("very", "robust", "obviously"). State what the system does, quantify it, and where relevant explain why it was built that way.
- Short paragraphs, one topic each, topic sentence first.
