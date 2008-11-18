.First.lib <- function(lib, pkg)  {
    library.dynam("vegan", pkg, lib)
    packageStartupMessage("This is vegan ",
                          utils::packageDescription("vegan", field="Version"),
                          appendLF = TRUE)
}
