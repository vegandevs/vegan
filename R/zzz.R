.onAttach <- function(lib, pkg)  {
    if(interactive())
        packageStartupMessage(
            "This is vegan devel version ",
            utils::packageDescription("vegan", fields="Version"),
            appendLF = TRUE)
}
