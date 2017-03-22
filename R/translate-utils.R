#' @include AnnotationFilter.R

## Functionality to translate a query condition to an AnnotationFilter.

#' Adapted from GenomicDataCommons.
#'
#' @importFrom methods is validObject initialize
#'
#' @noRd
.binary_op <- function(sep) {
    force(sep)
    function(e1, e2) {
        ## First create the class. Throws an error if not possible i.e. no
        ## class for the field available.
        field <- as.character(substitute(e1))
        class <- .fieldToClass(field)
        filter <- tryCatch({
            new(class, condition = sep, field = field)
        }, error = function(e) {
            stop("No AnnotationFilter class '", class, "' for field '",
                field, "' defined")
        })
        ## Fill with values.
        force(e2)
        if (is(filter, "CharacterFilter")) {
            e2 <- as.character(e2)
        } else if (is(filter, "IntegerFilter")) {
            e2 <- as.integer(e2)
        }
        initialize(filter, value = e2)
    }
}

#' Combine filters into a AnnotationFilterList combbined with \code{sep}
#'
#' @noRd
.combine_op <- function(sep) {
    force(sep)
    function(e1, e2) {
        ## Avoid implicit nesting of AnnotationFilterList - should be done
        ## eventually
        if (is(e1, "AnnotationFilterList")) {
            sep <- c(.logOp(e1), sep)
            e1 <- .aflvalue(e1)
        } else
            e1 <- list(e1)
        if (is(e2, "AnnotationFilterList")) {
            sep <- c(.logOp(e2), sep)
            e2 <- .aflvalue(e2)
        } else
            e2 <- list(e2)
        ## Don't use the constructor here.
        new("AnnotationFilterList", c(e1, e2), logOp = sep)
    }
}

#' The \code{.LOG_OP_REG} is an \code{environment} providing functions for
#' common logical operations to translate expressions into AnnotationFilter
#' objects.
#'
#' @noRd
.LOG_OP_REG <- local({
    filtReg <- new.env() # parent=emptyenv() breaks look-up

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

    filtReg
})

#' @rdname translate-utils
#'
#' @title Converting filter expressions into AnnotationFilters
#'
#' @description \code{convertFilterExpression} \emph{translates} a logical
#'     expression such as \code{gene_id == "BCL2"} into a filter object
#'     extending the \code{\link{AnnotationFilter}} class (in the example a
#'     \code{\link{GeneIdFilter}} object) or an
#'     \code{\link{AnnotationFilterList}} if the expression contains multiple
#'     conditions.
#'
#' @note No nesting of filter expressions (with \code{(}) is supported yet.
#'
#' @param expr an expression describing the filter rules to be applied. See
#'     examples below.
#'
#' @return \code{convertFilterExpression} and
#'     \code{convertFilterExpressionQuoted} return an
#'     \code{\link{AnnotationFilter}} or an \code{\link{AnnotationFilterList}}.
#'
#' @examples
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
#' @export
convertFilterExpression <-
    function(expr)
{
    x <- substitute(expr)
    convertFilterExpressionQuoted(x)
}

#' @rdname translate-utils
#'
#' @description \code{convertFilterExpressionQuoted} takes a \emph{quoted}
#'     filter expression (e.g. using \code{substitute}) and, as
#'     \code{convertFilterExpression}, translates it into an
#'     \code{\link{AnnotationFilter}} or \code{\link{AnnotationFilterList}}
#'     object.
#'
#' @details The \code{convertFilterExpression} function is designed to be used
#'     interactively, while the \code{convertFilterExpressionQuoted} is usually
#'     being called by other functions.
#'
#' @examples
#'
#' ## Define a simple function that calls the convertFilterExpressionQuoted
#' ## function
#' testFun <- function(x)
#'     convertFilterExpressionQuoted(substitute(x))
#'
#' ## Now we can use this function to translate a filter expression.
#' testFun(gene_id == 100)
#'
#' ## Alternatively we can call convertFilterExpressionQuoted passing the
#' ## quoted expression
#' filter_expr <- substitute(gene_id == 100)
#' convertFilterExpressionQuoted(filter_expr)
#' 
#' @export
convertFilterExpressionQuoted <-
    function(expr)
{
    eval(expr, envir = .LOG_OP_REG)
}
