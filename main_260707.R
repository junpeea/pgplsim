getwd()
usethis::proj_get()
usethis::use_git_config(
  user.name = "Yoonbae Jun",
  user.email = "junpeea@snu.ac.kr"
)
usethis::use_git()
usethis::use_github()

gitcreds::gitcreds_set()
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

Sys.which("git")

gitcreds::gitcreds_set()

gh::gh_whoami()

usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)

usethis::use_github(private = TRUE)

devtools::document()
devtools::load_all()
devtools::test()

usethis::use_testthat(edition = 3)

devtools::test(
  filter = "predict"
)

# I recommend this package-development workflow
#
# Every time you make changes:
#
#   cd C:\Users\yoonbaej\Documents\LuGPLSIM
#
# R
#
# Then in R:
#
#   devtools::document()
# devtools::check()
#
# Back in Command Prompt:
#
#   git status
# git add .
# git commit -m "renamed exported functions"
# git pull origin master
# git push

# This sequence keeps your local repository and GitHub repository synchronized while ensuring the package still passes R CMD check before each push.
