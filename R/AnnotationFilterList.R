#' @include AnnotationFilter.R

#' @name AnnotationFilterList
#'
#' @title Combining annotation filters
#'
#' @aliases AnnotationFilterList-class
#' 
#' @description The \code{AnnotationFilterList} allows to combine filter objects
#' extending the \code{\link{AnnotationFilter}} class to construct more complex
#' queries. Consecutive filter objects in the \code{AnnotationFilterList} can be
#' combined by a logical \emph{and} (\code{&}) or \emph{or} (\code{|}). The
#' \code{AnnotationFilterList} extends \code{list}, individual elements can thus
#' be accessed with \code{[[}.
#' 
#' @exportClass AnnotationFilterList
#' @rdname AnnotationFilterList
#' @name AnnotationFilterList
NULL

.AnnotationFilterList <- setClass(
    "AnnotationFilterList",
    contains = "list",
    slots = c(
        logOp = "character"
    ),
    prototype = prototype(
        logOp = character()
    )
)

.LOG_OPS <- c("&", "|")

setValidity("AnnotationFilterList", function(object) {
    txt <- character()
    vals <- object@.Data
    logOp <- object@logOp
    if (length(vals) == 0) {
        if (length(logOp))
            txt <- c(txt, "'logOp' can not have length > 0 if the object is empty")
    } else {
        ## Note: we allow length of vals being 1, but then logOp has to be empty.
        ## Check content:
        if (!all(unlist(lapply(vals, function(z) {
            is(z, "AnnotationFilter") | is(z, "AnnotationFilterList")
        }), use.names = FALSE)))
            txt <- c(txt, paste0("the object should contain only ",
                                 "'AnnotationFilter' or 'AnnotationFilterList' ",
                                 "objects"))
        ## Check that logOp has length object -1
        if (length(logOp) != length(vals) - 1)
            txt <- c(txt, "length of 'logOp' has to be length of the object -1")
        ## Check content of logOp.
        if (!all(logOp %in% .LOG_OPS))
            txt <- c(txt, "'logOp' can only contain '&' and '|'")
        if (length(txt)) txt else TRUE
    }
})

## AnnotationFilterList constructor function.
#' @rdname AnnotationFilterList
#' @name AnnotationFilterList
#' 
#' @param ... individual \code{\link{AnnotationFilter}} objects or a mixture of
#' \code{AnnotationFilter} and \code{AnnotationFilterList} objects.
#' 
#' @param logOp \code{character} of length being equal to the numner of submitted
#' \code{AnnotationFilter} objects -1. Each value representing the logical
#' operation to combine consecutive filters, i.e. the first element being the
#' logical operation to combine the first and second \code{AnnotationFilter}, the
#' second element being the logical operation to combine the second and third
#' \code{AnnotationFilter} and so on. Allowed values are \code{"&"} and
#' \code{"|"}. The function assumes a logical \emph{and} between all elements by
#' default.
#'
#' @seealso \code{\link{supportedFilters}} for available
#' \code{\link{AnnotationFilter}} objects
#' 
#' @examples
#'
#' ## Create some AnnotationFilters
#' gf <- GenenameFilter(c("BCL2", "BCL2L11"))
#' tbtf <- TxBiotypeFilter("protein_coding", condition = "!=")
#'
#' ## Combine both to an AnnotationFilterList. By default elements are combined
#' ## using a logical "and" operator. The filter list represents thus a query
#' ## like: get all features where the gene name is either ("BCL2" or "BCL2L11")
#' ## and the transcript biotype is not "protein_coding".
#' afl <- AnnotationFilterList(gf, tbtf)
#' afl
#'
#' ## Get the logical operator combining the filters
#' logOp(afl)
#'
#' ## and the list with the filter objects
#' value(afl)
#'
#' ## Access individual filters.
#' afl[[1]]
#'
#' ## Create a filter in the form of: get all features where the gene name is
#' ## either ("BCL2" or "BCL2L11") and the transcript biotype is not
#' ## "protein_coding" or the seq_name is "Y". Hence, this will get all feature
#' ## also found by the previous AnnotationFilterList and returns also all
#' ## features on chromosome Y.
#' afl <- AnnotationFilterList(gf, tbtf, SeqNameFilter("Y"), logOp = c("&", "|"))
#' afl
#' 
#' @export
AnnotationFilterList <- function(..., logOp = character()) {
    vals <- list(...)
    ## By default we're assuming & between elements.
    if (length(vals) > 1 & length(logOp) == 0)
        logOp <- rep("&", (length(vals) - 1))
    return(.AnnotationFilterList(vals, logOp = logOp))
}

#' @rdname AnnotationFilterList
#' @export
setMethod("show", "AnnotationFilterList", function(object) {
    cat("class:", class(object))
    if (length(object) == 0) {
        cat("empty object")
    } else {
        cat("\nfilters:\n")
        if (length(object) == 1) {
            show(object[[1]])  
        } else {
            cat("\n")
            show(object[[1]])
            for (i in seq_along(object@logOp)) {
                cat("", object@logOp[i], "\n")
                show(object[[i + 1]])
            }
        }
    }
})

#' @aliases logOp
#' @rdname AnnotationFilterList
#' 
#' @description \code{logOp()} get the logical operators that combine the
#' filters. Returns a \code{character} vector of length
#' \code{length(object) - 1}.
#'
#' @param object An \code{AnnotationFilterList}.
#' @export
setMethod("logOp", "AnnotationFilterList", function(object) {
    object@logOp
})

#' @rdname AnnotationFilterList
#'
#' @description \code{value()} get a \code{list} with the \code{AnnotationFilter}
#' objects. Use \code{[[} to access individual filters.
#' @export
setMethod("value", "AnnotationFilterList", function(object) {
    object@.Data
})
