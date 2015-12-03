.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is vegan ",
                          utils::packageDescription("vegan",
                                                    fields="Version"),
                          appendLF = TRUE)
}
