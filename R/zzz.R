.onAttach <- function(lib, pkg)  {
    packageStartupMessage("This is vegan ",
                          utils::packageDescription("vegan",
                                                    field="Version"),
                          appendLF = TRUE)
}
