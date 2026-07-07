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

devtools::check()
