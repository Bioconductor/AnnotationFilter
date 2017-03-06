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
    function(e1, e2) {
        ## Avoid implicit nesting of AnnotationFilterList - should be done
        ## eventually by ( later.
        if (is(e1, "AnnotationFilterList")) {
            sep <- c(logOp(e1), sep)
            e1 <- value(e1)
        } else
            e1 <- list(e1)
        if (is(e2, "AnnotationFilterList")) {
            sep <- c(logOp(e2), sep)
            e2 <- value(e2)
        } else
            e2 <- list(e2)
        ## Don't use the constructor here.
        afl <- new("AnnotationFilterList", c(e1, e2), logOp = sep)
        return(afl)
    }
}

.combine_op2 <- function(sep) {
    force(sep)
    function(e1, e2) {
        cat("2class e1: ", class(e1), "\n")
        cat("2class e2: ", class(e2), "\n")
        return(AnnotationFilterList(e1, e2, logOp = sep))
    }
}

## .group_op <- function(sep) {
##     force(sep)
##     cat("calling (\n")
##     ## do.call(sep, list(...))
##     function(...) {
##         do.call("(", list(...))
##     }
## }

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
## combine filters
filtReg$`&` <- .combine_op("&")
filtReg$`|` <- .combine_op("|")
## group filters
## filtReg$`(` <- .group_op("(")
assign(".LOG_OP_REG", filtReg, envir = topenv())


#' @rdname translate-utils
#'
#' @title Converting filter expressions into AnnotationFilters
#' 
#' @description \code{convertFilterExpression} \emph{translates} a logical
#' expression such as \code{gene_id == "BCL2"} into a filter object extending the
#' \code{\link{AnnotationFilter}} class (in the example a
#' \code{\link{GeneIdFilter}} object) or an \code{\link{AnnotationFilterList}} if
#' the expression contains multiple conditions.
#'
#' @note No nesting of filter expressions (with \code{(}) is supported yet.
#' 
#' @param expr an expression describing the filter rules to be applied. See
#' examples below.
#'
#' @examples
#'
#' ## Convert a filter expression based on a gene ID to a GeneIdFilter
#' gnf <- convertFilterExpression(gene_id == "BCL2")
#' gnf
#'
#' ## Same conversion but for two gene IDs.
#' gnf <- convertFilterExpression(gene_id %in% c("BCL2", "BCL2L11"))
#' gnf
#'
#' ## Converting an expression that combines multiple filters. As a result we
#' ## get an AnnotationFilterList containing the corresponding filters.
#' ## Be aware that nesting of expressions/filters does not work.
#' flt <- convertFilterExpression(gene_id %in% c("BCL2", "BCL2L11") &
#'                                tx_biotype == "nonsense_mediated_decay" |
#'                                seq_name == "Y")
#' flt
#'
#' @export
convertFilterExpression <- function(expr) {
    x <- substitute(expr)
    eval(x, envir = .LOG_OP_REG)
}

