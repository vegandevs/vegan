"ripley.subsets" <-
    function(n, r, s = 1:n)
{
    if(mode(n) != "numeric" || length(n) != 1
       || n < 1 || (n %% 1) != 0) stop("bad value of n")
    if(mode(r) != "numeric" || length(r) != 1
       || r < 1 || (r %% 1) != 0) stop("bad value of r")
    if(!is.atomic(s) || length(s) < n)
        stop("s is either non-atomic or too short")
    fun <- function(n, r, s)
        if(r <= 0) vector(mode(s), 0) else if(r >= n) s[1:n] else
    rbind(cbind(s[1], Recall(n - 1, r - 1, s[-1])),
          Recall(n - 1, r, s[-1]))
    fun(n, r, s)
}
