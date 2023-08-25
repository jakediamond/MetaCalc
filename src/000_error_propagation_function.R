# 
# Purpose: A function for propagating error
# Author: Jake Diamond via Lee Pang (https://www.r-bloggers.com/2015/01/easy-error-propagation-in-r/)
# Date: 17 October 2022
# 

mutate_with_error = function(.data, f) {
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    # expression to compute new variable errors
    sapply(all.vars(f[[3]]), function(v) {
      dfdp = deparse(D(f[[3]], v))
      sprintf('(d%s*(%s))^2', v, dfdp)
    }) %>%
      paste(collapse='+') %>%
      sprintf('sqrt(%s)', .)
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  .data %>%
    # the standard evaluation alternative of mutate()
    mutate_(.dots=exprs)
}
