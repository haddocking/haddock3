# Mathematical and Symbol Formatting Rules

Always format mathematical expressions, scientific notation, algorithmic complexities, and chemical/physical symbols as clean plain text, code spans, or standard Unicode. **Never use LaTeX dollar math syntax (`$...$` or `$$...$$`).**

## Formatting Standards

| Category | Do NOT Use (LaTeX) | REQUIRED Format | Examples |
| :--- | :--- | :--- | :--- |
| **Algorithmic Complexity** | `$O(N^2)$`, `$O(N \cdot M)$` | `O(N^2)` or `O(N²)` | `O(N^2) scaling`, `O(N * M)` |
| **Scientific Notation & Powers** | `$10^{-5}$`, `$10^4$` | `10^-5`, `1e-5`, or `10⁻⁵` | `1e-5 relative tolerance` |
| **Matrix Operations** | `$A \cdot A^T$`, `$A \times B$` | `A * A^T`, `A · A^T`, or `A x B` | `A * A^T dot product` |
| **Approximations & Ranges** | `$\approx 50,000,000$` | `~50,000,000` or `approx. 50,000,000` | `~50 million pairs` |
| **Metrics & Subscripts** | `$f_{nat}$`, `$i\text{-RMSD}$`, `$l\text{-RMSD}$` | `Fnat` or `f_nat`, `i-RMSD`, `l-RMSD` | `Fnat, i-RMSD, and l-RMSD` |
| **Units & Dimensions** | `$< 5.0\text{Å}$`, `$(N, M, 3)$` | `< 5.0 Å` or `(N, M, 3)` | `distance cutoff < 5.0 Å` |
| **Formulas / Fractions** | `$\frac{N(N-1)}{2}$` | `N * (N - 1) / 2` | `N * (N - 1) / 2 pairs` |

## Rationale
LaTeX math delimiters (`$...$`, `$$...$$`) often fail to render or clutter markdown readers, terminal outputs, and text diffs with backslashes and braces. Clean markdown formatting keeps documentation human-readable in plain text, code reviews, and rich previews.
