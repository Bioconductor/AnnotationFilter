#' @include AnnotationFilter.R

## Functionality to translate a query condition to an AnnotationFilter.

#' Adapted from GenomicDataCommons.
#' @importFrom methods is validObject
#' @noRd
.binary_op <- function(sep) {
    force(sep)
    function(e1, e2) {
        ## First create the class. Throws an error if not possible i.e. no
        ## class for the field available.
        field <- as.character(substitute(e1))
        clName <- .fieldToClass(field)
        tryCatch(
            newCl <- new(clName, condition = sep, field = field),
            error = function(e) {
                stop("No AnnotationFilter class ", clName, " (for field '",
                     field, "') defined!")
            }
        )
        ## Fill with values.
        force(e2)
        if (is(newCl, "CharacterFilter"))
            e2 <- as.character(e2)
        else if (is(newCl, "IntegerFilter"))
            e2 <- as.integer(e2)
        newCl@value <- e2
        if (validObject(newCl))
            return(newCl)
    }
}

#' Combine filters into a AnnotationFilterList combbined with \code{sep}
#' @noRd
.combine_op <- function(sep) {
    force(sep)
}

#' The \code{.LOG_OP_REG} is an \code{environment} providing functions for common
#' logical operations to translate expressions into AnnotationFilter objects.
#' 
#' @noRd
filtReg <- new.env()
## Assign conditions.
filtReg$`==` <- .binary_op("==")
filtReg$`%in%` <- .binary_op("==")
filtReg$`!=` <- .binary_op("!=")
filtReg$`>` <- .binary_op(">")
filtReg$`<` <- .binary_op("<")
filtReg$`>=` <- .binary_op(">=")
filtReg$`<=` <- .binary_op("<=")
assign(".LOG_OP_REG", filtReg, envir = topenv())

#' That's to test the behaviour.
#' @param expr Is expected to be an expression.
#' @noRd
.expressionToAnnotationFilter <- function(expr) {
    x <- substitute(expr)
    eval(x, envir = .LOG_OP_REG)
}

