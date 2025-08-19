# Contributing

Thanks for considering contributing to **cleanRfish**!

## How to contribute
1. Fork the repo and create a feature branch.
2. Install dev dependencies in R: `install.packages(c("devtools","testthat","roxygen2"))`
3. Run checks: `devtools::check()`
4. Add tests under `tests/testthat/`.
5. Submit a pull request with a clear description.

## Reporting bugs
Please open an issue with a minimal **reprex** and session info.

## Style
- Document with **roxygen2** (markdown = TRUE).
- Prefer tidyverse style guides.
