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
			  not = "logical")
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
        ## Check that all elements are non-empty (issue #17). Doing this
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
#' @seealso \code{\link{supportedFilters}} for available
#'     \code{\link{AnnotationFilter}} objects
#'
#' @return \code{AnnotationFilterList} returns an \code{AnnotationFilterList}.
#' 
#' @examples
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
    function(..., logicOp = character(), logOp = character(), not = FALSE)
{
    if (!missing(logOp) && missing(logicOp)) {
        logicOp <- logOp
        .Deprecated(msg = "'logOp' deprecated, use 'logicOp'")
    }
    filters <- list(...)
	## Remove empty nested AnnotationFilterLists
	for (i in seq_len(length(filters))) {
		if(is(filters[[i]], "AnnotationFilterList") && length(filters[[i]]) == 0)
			filters <- filters[-i]
	}
    if (length(filters) > 1 & length(logicOp) == 0)
    ## By default we're assuming & between elements.
        logicOp <- rep("&", (length(filters) - 1))
    .AnnotationFilterList(filters, logOp = logicOp, not = not)
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

#' @rdname AnnotationFilterList
#'
#' @aliases subset
#'
#' @description \code{[} gets the logical operators separating
#'     successive \code{AnnotationFilter}.
#'
#' @return \code{[} returns a \code{character()} vector of
#'     \dQuote{&} or \dQuote{|} symbols.
#'
#' @export
setMethod('[', 'AnnotationFilterList', function(x, i, j, ..., drop=TRUE) {
	res <- callNextMethod()
	i <- seq(length(x)) %in% i
	ops <- .logicOp_subset(logicOp(x), i)
	do.call(AnnotationFilterList, c(res, logicOp=ops, not=T))
})

.logicOp_subset <- function(op, fields_subset) {
	keepOp <- rep(TRUE, length(op))

    first <- seq_along(fields_subset) - 1L
    second <- seq_along(fields_subset)

    isFirst <- first == 0L
    keepOp[ which(!fields_subset & isFirst)  ] <- FALSE

    isLast <- second == length(fields_subset)
    keepOp[ which(!fields_subset & isLast) - 1 ] <- FALSE

    isOther <- !(isFirst | isLast)
    isDifferent <- c(FALSE, op[first] != op[second[-length(second)]])

    keepOp[ second[isOther & !fields_subset & c(FALSE, op[first] == "&")] ] <- FALSE
    keepOp[  first[isOther & !fields_subset & c(FALSE, op[first] != "&")] ] <- FALSE
    keepOp[ second[isOther & isDifferent] ] <- FALSE

    if(!any(fields_subset[-1]) || !any(head(fields_subset, -1)))   ## Catch an edge case
		keepOp <- rep(FALSE, length(op))

    if(length(fields_subset) >= 4 &&
               !any(fields_subset[2:(length(fields_subset)-1)]) &&
               fields_subset[1] == TRUE && fields_subset[length(fields_subset)]) {
        if(any(op[c(1, length(op))] == c('|', '|')))
            return('|')
        else
            return('&')
    }

    op[keepOp]
}

#' @rdname AnnotationFilterList
#'
#' @param object An object of class \code{AnnotationFilterList}.
#'
#' @importFrom utils tail
#' @export
setMethod("show", "AnnotationFilterList",
    function(object)
{
    cat(
        "class: ", class(object), "\n",
        "length: ", length(object), "\n",
        sep = ""
    )
	if(not(object))
		cat("NOT\n")
    if (length(object)) {
        cat("filters:\n\n")
        show(object[[1]])
        for (i in tail(seq_along(object), -1L)) {
            cat("\n", logicOp(object)[i - 1L], "\n\n")
            show(object[[i]])
        }
    }
})

