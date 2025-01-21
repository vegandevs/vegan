.onAttach <- function(lib, pkg)  {
    if(interactive())
        packageStartupMessage(
            "This is vegan ",
            utils::packageDescription("vegan", fields="Version"),
            appendLF = TRUE)
}
