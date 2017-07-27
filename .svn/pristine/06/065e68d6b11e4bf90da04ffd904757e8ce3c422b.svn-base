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
            sep <- c(logicOp(e1), sep)
            e1 <- .aflvalue(e1)
        } else
            e1 <- list(e1)
        if (is(e2, "AnnotationFilterList")) {
            sep <- c(logicOp(e2), sep)
            e2 <- .aflvalue(e2)
        } else
            e2 <- list(e2)
        ## Don't use the constructor here.
        new("AnnotationFilterList", c(e1, e2), logOp = sep)
    }
}

#' The \code{.LOG_OP_REG} is a \code{list} providing functions for
#' common logical operations to translate expressions into AnnotationFilter
#' objects.
#'
#' @noRd
.LOG_OP_REG <- list()
## Assign conditions.
.LOG_OP_REG$`==` <- .binary_op("==")
.LOG_OP_REG$`%in%` <- .binary_op("==")
.LOG_OP_REG$`!=` <- .binary_op("!=")
.LOG_OP_REG$`>` <- .binary_op(">")
.LOG_OP_REG$`<` <- .binary_op("<")
.LOG_OP_REG$`>=` <- .binary_op(">=")
.LOG_OP_REG$`<=` <- .binary_op("<=")
## combine filters
.LOG_OP_REG$`&` <- .combine_op("&")
.LOG_OP_REG$`|` <- .combine_op("|")

#' @rdname AnnotationFilter
#'
#' @description \code{AnnotationFilter} \emph{translates} a filter
#'     expression such as \code{~ gene_id == "BCL2"} into a filter object
#'     extending the \code{\link{AnnotationFilter}} class (in the example a
#'     \code{\link{GeneIdFilter}} object) or an
#'     \code{\link{AnnotationFilterList}} if the expression contains multiple
#'     conditions (see examples below). Filter expressions have to be written
#'     in the form \code{~ <field> <condition> <value>}, with \code{<field>}
#'     being the default field of the filter class (use the
#'     \code{supportedFilter} function to list all fields and filter classes),
#'     \code{<condition>} the logical expression and \code{<value>} the value
#'     for the filter.
#'
#' @details Filter expressions for the \code{AnnotationFilter} class have to be
#'     written as formulas, i.e. starting with a \code{~}.
#'
#' @note Translation of nested filter expressions using the
#'     \code{AnnotationFilter} function is not yet supported.
#' 
#' @param expr A filter expression, written as a \code{formula}, to be
#'     converted to an \code{AnnotationFilter} or \code{AnnotationFilterList}
#'     class. See below for examples.
#'
#' @return \code{AnnotationFilter} returns an
#'     \code{\link{AnnotationFilter}} or an \code{\link{AnnotationFilterList}}.
#' 
#' @importFrom lazyeval f_eval
#'
#' @examples
#' 
#' ## Convert a filter expression based on a gene ID to a GeneIdFilter
#' gnf <- AnnotationFilter(~ gene_id == "BCL2")
#' gnf
#'
#' ## Same conversion but for two gene IDs.
#' gnf <- AnnotationFilter(~ gene_id %in% c("BCL2", "BCL2L11"))
#' gnf
#'
#' ## Converting an expression that combines multiple filters. As a result we
#' ## get an AnnotationFilterList containing the corresponding filters.
#' ## Be aware that nesting of expressions/filters does not work.
#' flt <- AnnotationFilter(~ gene_id %in% c("BCL2", "BCL2L11") &
#'                         tx_biotype == "nonsense_mediated_decay" |
#'                         seq_name == "Y")
#' flt
#' 
#' @export
AnnotationFilter <- function(expr) {
    f_eval(expr, data = .LOG_OP_REG)
}
