# pigCRCvirome

Code and processed data for **Gut virome restructuring in a porcine model of Western diet-associated colorectal cancer**. Experimental background is described in the parent study by [Osswald, Wortmann et al.](https://doi.org/10.1136/gutjnl-2024-332243).

The study includes 73 virome libraries from 16 APC-mutant pigs across three diet groups: control feed (CTRL, n = 5), red meat and lard (RL, n = 5), and RL plus colestyramine (COL, n = 6). Raw sequence accessions are available in ENA under [PRJEB121175](https://www.ebi.ac.uk/ena/browser/view/PRJEB121175) and [PRJEB123405](https://www.ebi.ac.uk/ena/browser/view/PRJEB123405).

## Reproduction

The analysis environment uses R, [Quarto](https://quarto.org/), and [qproj](https://github.com/rujinlong/qproj), with dependencies recorded in [renv.lock](renv.lock).

1. Restore package dependencies into the project library:

```sh
Rscript -e 'install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(library = renv::paths$library(), prompt = FALSE)'
```

2. Execute all analyses and verification checks:

```sh
Rscript -e 'renv::load(); source(here::here("scripts", "reproduce.R"))'
```

`scripts/reproduce.R` renders the Quarto reports, writes output tables and publication figures to `analyses/data/`, and runs [scripts/verify.R](scripts/verify.R) to validate statistics and figures against the submitted source data.

## Analyses

Analyses are implemented in Quarto notebooks under `analyses/`:

| Module | Analysis |
|---|---|
| [01-community](analyses/01-community.qmd) | Bray-Curtis PERMANOVA, dispersion, ordination, paired change, and litter-restricted tests |
| [02-lifestyle](analyses/02-lifestyle.qmd) | Temperate phage enrichment and LinDA differential abundance models |
| [03-host-associations](analyses/03-host-associations.qmd) | Candidate phages (bile acid metabolism), host-family correlations, and sensitivity tests |
| [04-figures](analyses/04-figures.qmd) | Figures 1 and 2 |

Input matrices, genome annotations, and clinical tables are documented in [analyses/data/00-raw/README.md](analyses/data/00-raw/README.md).

## License

GPL-3.0-or-later; see [LICENSE](LICENSE).
