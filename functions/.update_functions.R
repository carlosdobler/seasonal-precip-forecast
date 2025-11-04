# MODULES
rts <- c("general_tools")

# *****

if (!fs::dir_exists("functions")) {
  fs::dir_create("functions")
}

purrr::walk(rts, \(rt) {
  download.file(
    stringr::str_glue(
      "https://raw.github.com/carlosdobler/spatial-routines/master/{rt}.R"
    ),
    stringr::str_glue("functions/{rt}.R"),
    quiet = T
  )
})
