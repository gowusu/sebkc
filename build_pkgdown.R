## build_pkgdown.R
## Builds the pkgdown documentation website for sebkc and (optionally) deploys it
## to the gh-pages branch. The two manuals in manual/*.md are the single source of
## truth; this script regenerates them as pkgdown articles under vignettes/articles/
## (a build artifact, git-ignored) so nothing is duplicated on master and no HTML
## ever lands on the default branch.
##
##   Rscript build_pkgdown.R          # build only, into docs/ (for a local check)
##   Rscript build_pkgdown.R deploy   # build AND push to the gh-pages branch
##
## Requires the pkgdown package and a git executable.
deploy <- "deploy" %in% commandArgs(trailingOnly = TRUE)

## make the GitHub Desktop git visible to R (it is not on PATH here)
gitdir <- "C:/Users/PC/AppData/Local/GitHubDesktop/app-3.6.1/resources/app/git/cmd"
if (dir.exists(gitdir) && Sys.which("git") == "")
  Sys.setenv(PATH = paste(gitdir, Sys.getenv("PATH"), sep = .Platform$path.sep))

## point R/rmarkdown at the standalone Pandoc (Rscript does not bundle it)
pandocdir <- "C:/Users/PC/AppData/Local/Pandoc"
if (dir.exists(pandocdir)) {
  Sys.setenv(RSTUDIO_PANDOC = pandocdir)
  if (!nzchar(Sys.which("pandoc")))
    Sys.setenv(PATH = paste(pandocdir, Sys.getenv("PATH"), sep = .Platform$path.sep))
}

pkg <- "C:/Users/PC/Documents/GitHub/sebkc"
art <- file.path(pkg, "vignettes", "articles")
dir.create(file.path(art, "figures"), recursive = TRUE, showWarnings = FALSE)

## turn a manual/<name>.md into a pkgdown article vignettes/articles/<name>.Rmd
make_article <- function(mdfile, rmdname) {
  x <- readLines(mdfile, warn = FALSE, encoding = "UTF-8")
  h1 <- grep("^# ", x)[1]
  title <- sub("^# ", "", x[h1])
  body <- x[-h1]                                   # drop H1 (YAML supplies the title)
  ## repoint cross-links between the two manuals to their article pages
  body <- gsub("manual/crop_water_needs.md", "crop_water_needs.html", body, fixed = TRUE)
  body <- gsub("manual/sebkc_manual.md",     "sebkc_manual.html",     body, fixed = TRUE)
  yaml <- c("---",
            paste0('title: "', gsub('"', "'", title), '"'),
            "---", "")
  writeLines(c(yaml, body), file.path(art, rmdname), useBytes = TRUE)
}

make_article(file.path(pkg, "manual", "crop_water_needs.md"), "crop_water_needs.Rmd")
make_article(file.path(pkg, "manual", "sebkc_manual.md"),     "sebkc_manual.Rmd")

## figures live next to the articles so the `figures/...` links resolve
figs <- list.files(file.path(pkg, "manual", "figures"), full.names = TRUE)
file.copy(figs, file.path(art, "figures"), overwrite = TRUE)

if (deploy) {
  message("Building and deploying to gh-pages ...")
  pkgdown::deploy_to_branch(pkg, branch = "gh-pages", new_process = FALSE)
} else {
  message("Building site into docs/ (local check only) ...")
  pkgdown::build_site(pkg, preview = FALSE, install = FALSE, new_process = FALSE)
}
message("done.")
