# VS Code Markdown Setup for PhD Research

This note outlines a highly optimized VS Code extension stack tailored for academic writing, reference management, literature synthesis, and manuscript compilation.

## 📚 Citation & Reference Management

* **[Zotero Better BibTeX Citations](https://github.com)**: Connects your local Zotero library directly to your Markdown files. Pressing `Alt + Shift + Z` opens a citation picker UI to search your library and automatically insert formatted Pandoc-style citation keys (e.g., `[@author2026]`). *(Requires the Zotero desktop app and the Better BibTeX plugin)*.
* **[Zotex](https://visualstudio.com)**: Instantly pulls and autocompletes citation keys from your `.bib` files. Supports Pandoc citation styles, traditional LaTeX `\cite{}` commands, and Markdown footnotes.

## 🗂️ Personal Knowledge Management & Outlining

* **[Foam](https://visualstudio.com)**: Transforms VS Code into a "second brain" identical to Obsidian or Roam Research. It links concepts using `[[WikiLinks]]`, auto-generates interactive graph visualizations of your notes, and manages back-links to track literature reviews and cross-reference themes.
* **[Dendron](https://visualstudio.com)**: Specialized for deep, hierarchical note management. It organizes highly structured data (e.g., naming notes `methods.lab-protocol.cell-culture`) to keep massive quantities of field notes or archive files searchable.

## 🔬 Scientific Drafting Tools

* **[Rxiv-Maker](https://henriqueslab.org)**: Provides a targeted framework for writing scientific manuscripts. Offers specialized cross-referencing capabilities (for equations, figures, and tables), autocompletes math definitions, and ensures drafting meets formal preprint formatting benchmarks.
* **[Science Writer](https://visualstudio.com)**: Offers pre-built boilerplate snippets and structural phrasing shortcuts. Acts as a useful blueprinting assistant when structuring abstract, methodology, or discussion sections.

## 🖨️ Thesis Compiling & Formatted Exporting

* **[Quarto](https://visualstudio.com)**: The gold standard for modern academic publishing. Combines plain Markdown text with live executable code blocks (Python, R, or Julia) to generate dynamic, publication-ready PDFs, Word documents, or HTML templates with rendered bibliographies.
* **[Markdown Preview Enhanced](https://visualstudio.com)**: Essential for reviewing intensive typesetting. Flawlessly renders complete LaTeX math blocks ($$E=mc^2$$), footnotes, and internal document anchors directly alongside your workspace.

---

### Recommended Academic Workflow
1. **Gather & Manage**: Zotero (Desktop) + Better BibTeX
2. **Brainstorm & Synthesise**: VS Code + Foam (`[[WikiLinks]]` for literature notes)
3. **Draft & Compile**: Quarto (for seamless embedding of data analysis and bibliography rendering)
