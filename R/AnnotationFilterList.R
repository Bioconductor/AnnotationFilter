#' @include AnnotationFilter.R

#' @rdname AnnotationFilterList
#'
#' @name AnnotationFilterList
#'
#' @title Combining annotation filters
#'
#' @aliases AnnotationFilterList-class
#'
#' @description The \code{AnnotationFilterList} allows to combine
#'     filter objects extending the \code{\link{AnnotationFilter}}
#'     class to construct more complex queries. Consecutive filter
#'     objects in the \code{AnnotationFilterList} can be combined by a
#'     logical \emph{and} (\code{&}) or \emph{or} (\code{|}). The
#'     \code{AnnotationFilterList} extends \code{list}, individual
#'     elements can thus be accessed with \code{[[}.
#'
#' @note The \code{AnnotationFilterList} does not support containing empty
#'     elements, hence all elements of \code{length == 0} are removed in
#'     the constructor function.
#'
#' @exportClass AnnotationFilterList
NULL

.AnnotationFilterList <- setClass(
    "AnnotationFilterList",
    contains = "list",
    slots = c(logOp = "character",
              not = "logical",
              .groupingFlag = "logical")
)

.LOG_OPS <- c("&", "|")

setValidity("AnnotationFilterList",
    function(object)
{
    txt <- character()
    filters <- .aflvalue(object)
    logOp <- .logOp(object)
    not <- .not(object)
    if (length(filters) == 0 && length(logOp)) {
        txt <- c(
            txt, "'logicOp' can not have length > 0 if the object is empty"
        )
    } else if (length(filters) != 0) {
        ## Note: we allow length of filters being 1, but then logOp has
        ## to be empty.  Check content:
        fun <- function(z)
            is(z, "AnnotationFilter") || is(z, "AnnotationFilterList")
        test <- vapply(filters, fun, logical(1))
        if (!all(test)){
            txt <- c(
                txt, "only 'AnnotationFilter' or 'AnnotationFilterList' allowed"
            )
        }
        # Check that all elements are non-empty (issue #17). Doing this
        ## separately from the check above to ensure we get a different error
        ## message.
        if (!all(lengths(filters) > 0))
            txt <- c(txt, "Lengths of all elements have to be > 0")
        ## Check that logOp has length object -1
        if (length(logOp) != length(filters) - 1)
            txt <- c(txt, "length of 'logicOp' has to be length of the object -1")
        ## Check content of logOp.
        if (!all(logOp %in% .LOG_OPS))
            txt <- c(txt, "'logicOp' can only contain '&' and '|'")
    }

    if (length(txt)) txt else TRUE
})

## AnnotationFilterList constructor function.
#' @rdname AnnotationFilterList
#'
#' @name AnnotationFilterList
#'
#' @param ... individual \code{\link{AnnotationFilter}} objects or a
#'     mixture of \code{AnnotationFilter} and
#'     \code{AnnotationFilterList} objects.
#'
#' @param logicOp \code{character} of length equal to the number
#'     of submitted \code{AnnotationFilter} objects - 1. Each value
#'     representing the logical operation to combine consecutive
#'     filters, i.e. the first element being the logical operation to
#'     combine the first and second \code{AnnotationFilter}, the
#'     second element being the logical operation to combine the
#'     second and third \code{AnnotationFilter} and so on. Allowed
#'     values are \code{"&"} and \code{"|"}. The function assumes a
#'     logical \emph{and} between all elements by default.
#'
#' @param logOp Deprecated; use \code{logicOp=}.
#'
#' @param .groupingFlag Flag desginated for internal use only.
#'
#' @param not \code{logical} of length one. Indicates whether the grouping
#'      of \code{AnnotationFilters} are to be negated.
#'
#' @seealso \code{\link{supportedFilters}} for available
#'     \code{\link{AnnotationFilter}} objects
#'
#' @return \code{AnnotationFilterList} returns an \code{AnnotationFilterList}.
#' 
#' @examples
#' ## Create some AnnotationFilters
#' gf <- GeneNameFilter(c("BCL2", "BCL2L11"))
#' tbtf <- TxBiotypeFilter("protein_coding", condition = "!=")
#'
#' ## Combine both to an AnnotationFilterList. By default elements are combined
#' ## using a logical "and" operator. The filter list represents thus a query
#' ## like: get all features where the gene name is either ("BCL2" or "BCL2L11")
#' ## and the transcript biotype is not "protein_coding".
#' afl <- AnnotationFilterList(gf, tbtf)
#' afl
#'
#' ## Access individual filters.
#' afl[[1]]
#'
#' ## Create a filter in the form of: get all features where the gene name is
#' ## either ("BCL2" or "BCL2L11") and the transcript biotype is not
#' ## "protein_coding" or the seq_name is "Y". Hence, this will get all feature
#' ## also found by the previous AnnotationFilterList and returns also all
#' ## features on chromosome Y.
#' afl <- AnnotationFilterList(gf, tbtf, SeqNameFilter("Y"),
#'                             logicOp = c("&", "|"))
#' afl
#'
#' @export
AnnotationFilterList <-
    function(..., logicOp = character(), logOp = character(), not = FALSE,
        .groupingFlag=FALSE)
{
    if (!missing(logOp) && missing(logicOp)) {
        logicOp <- logOp
        .Deprecated(msg = "'logOp' deprecated, use 'logicOp'")
    }
    filters <- list(...)

    ## Remove empty nested lists and AnnotationFilterLists
    removal <- lengths(filters) != 0
    filters <- filters[removal]

    if (length(filters) > 1 & length(logicOp) == 0)
        ## By default we're assuming & between elements.
        logicOp <- rep("&", (length(filters) - 1))
    .AnnotationFilterList(filters, logOp = logicOp, not = not,
        .groupingFlag=.groupingFlag)
}

.logOp <- function(object) object@logOp

.aflvalue <- function(object) object@.Data

.not <- function(object) object@not

#' @rdname AnnotationFilterList
#'
#' @description \code{value()} get a \code{list} with the
#'     \code{AnnotationFilter} objects. Use \code{[[} to access
#'     individual filters.
#'
#' @return \code{value()} returns a \code{list} with \code{AnnotationFilter}
#'     objects.
#' 
#' @export
setMethod("value", "AnnotationFilterList", .aflvalue)

#' @rdname AnnotationFilterList
#'
#' @aliases logicOp
#'
#' @description \code{logicOp()} gets the logical operators separating
#'     successive \code{AnnotationFilter}.
#'
#' @return \code{logicOp()} returns a \code{character()} vector of
#'     \dQuote{&} or \dQuote{|} symbols.
#'
#' @export logicOp
setMethod("logicOp", "AnnotationFilterList", .logOp)

#' @rdname AnnotationFilterList
#'
#' @aliases not
#'
#' @description \code{not()} gets the logical operators separating
#'     successive \code{AnnotationFilter}.
#'
#' @return \code{not()} returns a \code{character()} vector of
#'     \dQuote{&} or \dQuote{|} symbols.
#'
#' @export not
setMethod("not", "AnnotationFilterList", .not)

.distributeNegation <- function(object, .prior_negation=FALSE)
{
    if(.not(object))
        .prior_negation <- ifelse(.prior_negation, FALSE, TRUE)
    filters <- lapply(object, function(x){
        if(is(x, "AnnotationFilterList"))
            distributeNegation(x, .prior_negation)   
        else{
            if(.prior_negation) x@not <- ifelse(x@not, FALSE, TRUE)
            x
        }
    })
    ops <- vapply(logicOp(object), function(x) {
        if(.prior_negation){
            if(x == '&')
                '|'
            else
                '&'
        }
        else
            x
    }
        ,character(1)
    )
    ops <- unname(ops)
    filters[['logicOp']] <- ops
    do.call("AnnotationFilterList", filters)
}

#' @rdname AnnotationFilterList
#'
#' @aliases distributeNegation
#'
#' @description
#'
#' @param .prior_negation \code{logical(1)} unused argument.
#'
#' @return \code{AnnotationFilterList} object with DeMorgan's law applied to
#'      it such that it is equal to the original \code{AnnotationFilterList}
#'      object but all \code{!}'s are distributed out of the
#'      \code{AnnotationFilterList} object and to the nested
#'      \code{AnnotationFilter} objects.
#'
#' @examples
#' afl <- AnnotationFilter(~!(symbol == 'ADA' | symbol %startsWith% 'SNORD'))
#' afl <- distributeNegation(afl)
#' afl
#' @export
setMethod("distributeNegation", "AnnotationFilterList", .distributeNegation)

#' @importFrom utils head
#'
#' @noRd
.convertFilterList <- function(object, show, granges=list(), nested=FALSE)
{
    filters <- value(object)
    result <- character(length(filters))
    for (i in seq_len(length(filters))) {
        if (is(filters[[i]], "AnnotationFilterList")) {
            res <- .convertFilterList(filters[[i]], show=show, granges=granges,
                nested=TRUE)
            granges <- c(granges, res[[2]])
            result[i] <- res[[1]]
        }
        else if (field(filters[[i]]) == "granges") {
            if(!show)
                result[i] <- .convertFilter(filters[[i]])
            else {
                nam <- paste0("GRangesFilter_", length(granges) + 1)
                granges <- c(granges, list(filters[[i]]))
                result[i] <- nam
            }
        }
        else
            result[i] <- .convertFilter(filters[[i]])
    }

    result_last <- tail(result, 1)
    result <- head(result, -1)
    result <- c(rbind(result, logicOp(object)))
    result <- c(result, result_last)
    result <- paste(result, collapse=" ")
    if(nested || object@not)
        result <- paste0("(", result, ")")
    if(object@not)
        result <- paste0("!", result)
        
    list(result, granges)
}

#' @rdname AnnotationFilterList
#'
#' @aliases convertFilter
#'
#' @description Converts an \code{AnnotationFilterList} object to a
#'      \code{character(1)} giving an equation that can be used as input to
#'      a \code{dplyr} filter.
#'
#' @return \code{character(1)} that can be used as input to a \code{dplyr}
#'      filter.
#'
#' @examples
#' afl <- AnnotationFilter(~symbol=="ADA" & tx_start > "400000")
#' result <- convertFilter(afl)
#' result
#' @export
setMethod("convertFilter", signature(object = "AnnotationFilterList",
                                     db = "missing") , function(object)
{
    result <- .convertFilterList(object, show=FALSE)
    result[[1]]
})

#' @rdname AnnotationFilterList
#'
#' @param object An object of class \code{AnnotationFilterList}.
#'
#' @importFrom utils tail
#' @export
setMethod("show", "AnnotationFilterList", function(object)
{
    result <- .convertFilterList(object, show=TRUE)
    granges <- result[[2]]
    result <- result[[1]]
    cat("AnnotationFilterList of length", length(object), "\n")
    cat(result)
    cat("\n")
    for(i in seq_len(length(granges))) {
        cat("\n")
        cat("Symbol: GRangesFilter_", i, "\n", sep="")
        show(granges[[1]])
        cat("\n")
    }
})
