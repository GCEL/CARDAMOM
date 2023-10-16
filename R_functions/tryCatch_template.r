
# Example error handling using tryCatch from base R
# Copied from https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
# thank you the R community

#urls <- c(
#    "http://stat.ethz.ch/R-manual/R-devel/library/base/html/connections.html",
#    "http://en.wikipedia.org/wiki/Xz",
#    "xxxxx"
#)
#readUrl <- function(url) {
#    out <- tryCatch(
#        {
#            # Just to highlight: if you want to use more than one
#            # R expression in the "try" part then you'll have to
#            # use curly brackets.
#            # 'tryCatch()' will return the last evaluated expression
#            # in case the "try" part was completed successfully
#
#            message("This is the 'try' part")
#
#            readLines(con=url, warn=FALSE)
#            # The return value of `readLines()` is the actual value
#            # that will be returned in case there is no condition
#            # (e.g. warning or error).
#            # You don't need to state the return value via `return()` as code
#            # in the "try" part is not wrapped inside a function (unlike that
#            # for the condition handlers for warnings and error below)
#        },
#        error=function(cond) {
#            message(paste("URL does not seem to exist:", url))
#            message("Here's the original error message:")
#            message(cond)
#            # Choose a return value in case of error
#            return(NA)
#        },
#        warning=function(cond) {
#            message(paste("URL caused a warning:", url))
#            message("Here's the original warning message:")
#            message(cond)
#            # Choose a return value in case of warning
#            return(NULL)
#        },
#        finally={
#        # NOTE:
#        # Here goes everything that should be executed at the end,
#        # regardless of success or error.
#        # If you want more than one expression to be executed, then you
#        # need to wrap them in curly brackets ({...}); otherwise you could
#        # just have written 'finally=<expression>'
#            message(paste("Processed URL:", url))
#            message("Some other message at the end")
#        }
#    )
#    return(out)
#}
#
#tryCatch has a slightly complex syntax structure. However, once we understand the 4 parts which constitute a complete tryCatch call as shown below, it becomes easy to remember:
#
#expr: [Required] R code(s) to be evaluated
#
#error : [Optional] What should run if an error occured while evaluating the codes in expr
#
#warning : [Optional] What should run if a warning occured while evaluating the codes in expr
#
#finally : [Optional] What should run just before quitting the tryCatch call, irrespective of if expr ran successfully, with an error, or with a warning
#
#tryCatch(
#    expr = {
#        # Your code...
#        # goes here...
#        # ...
#    },
#    error = function(e){
#        # (Optional)
#        # Do this if an error is caught...
#    },
#    warning = function(w){
#        # (Optional)
#        # Do this if an warning is caught...
#    },
#    finally = {
#        # (Optional)
#        # Do this at the end before quitting the tryCatch structure...
#    }
#)
#Thus, a toy example, to calculate the log of a value might look like:
#
#log_calculator <- function(x){
#    tryCatch(
#        expr = {
#            message(log(x))
#            message("Successfully executed the log(x) call.")
#        },
#        error = function(e){
#            message('Caught an error!')
#            print(e)
#        },
#        warning = function(w){
#            message('Caught an warning!')
#            print(w)
#        },
#        finally = {
#            message('All done, quitting.')
#        }
#    )
#}
#Now, running three cases:
#
#A valid case
#
#log_calculator(10)
## 2.30258509299405
## Successfully executed the log(x) call.
## All done, quitting.
#A "warning" case
#
#log_calculator(-10)
## Caught an warning!
## <simpleWarning in log(x): NaNs produced>
## All done, quitting.
#An "error" case
#
#log_calculator("log_me")
## Caught an error!
## <simpleError in log(x): non-numeric argument to mathematical function>
## All done, quitting.
