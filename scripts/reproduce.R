# Run from the repository root: Rscript scripts/reproduce.R
here::i_am("scripts/reproduce.R")
modules <- c("01-community", "02-lifestyle", "03-host-associations", "04-figures")
for (module in modules) {
  status <- system2("quarto", c("render", shQuote(here::here("analyses", paste0(module, ".qmd"))),
                                "--execute", "--no-cache"))
  if (status != 0L) stop("Rendering failed: ", module)
}
source(here::here("scripts", "verify.R"))
