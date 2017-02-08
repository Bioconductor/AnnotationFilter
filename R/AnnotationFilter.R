.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CHAR_FIELDS <- c("symbol")

.INT_FIELDS <- character()

#' @rdname AnnotationFilter
#' @export
.AnnotationFilter <- setClass(
    "AnnotationFilter",
    representation(
        "VIRTUAL",
        field="character",
        condition="character",
        value="ANY",
        .valueIsCharacter="logical"
    ),
    prototype=list(
        condition= "==",
        value=character(),
        .valueIsCharacter=TRUE
    )
)

setValidity("AnnotationFilter", function(object) {
    value <- .value(object)
    condition <- .condition(object)
    isCharacter <- .isCharacter(object)
    txt <- character()
    condNa <- any(is.na(condition))
    condOne <- length(condition) == 1L
    ## Check condition
    if (!condOne)
        txt <- c(txt, "'condition' must be length 1")
    if (condNa)
        txt <- c(txt, "'condition' can not be NA")
    if (!condNa && condOne && !condition %in% .OPS)
        txt <- c(txt,
                 sprintf("'condition' must be one of %s",
                         paste("'", .OPS, "'", collapse=", ")))
    if (!condNa && condOne && condition  %in%
        c("startsWith", "endsWith", ">", "<", ">=", "<=") && length(value) > 1L)
        txt <- c(txt,
                 paste0("'value' must be length 1 when condition is '",
                        condition, "'"))
    if (!condNa && condOne && condition  %in% c("startsWith", "endsWith") &&
        !isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with character value"))
    if (!condNa && condOne && condition  %in% c(">", "<", ">=", "<=") &&
        isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with integer value"))
    ## Check value
    if (any(is.na(value)))
        txt <- c(txt, "'value' can not be NA")
    if (isCharacter && !is.character(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take character value"))
    if (!isCharacter && !is.integer(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take integer value"))
    if (length(txt)) txt else TRUE
})

## slot accessors

.field <- function(x) x@field

.condition <- function(x) x@condition

.value <- function(x) x@value

.isCharacter <- function(x) x@.valueIsCharacter

## helper functions

.fieldToClass <- function(field) {
    class <- sub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, if (length(class)) "Filter" else character(0))
}


#' @importFrom methods new
.filterFactory <- function(field, class, .valueIsCharacter) {
    force(field); force(class)          # watch for lazy evaluation
    as.value <-
        if (.valueIsCharacter) {
            as.character
        } else {
            function(x) {
                stopifnot(is.numeric(x))
                as.integer(x)
            }
        }
    function(value, condition = "==") {
        value <- as.value(value)
        new(class, field=field, condition = as.character(condition),
            value=value, .valueIsCharacter=.valueIsCharacter)
    }
}

## Install-time class creation

local({
    field <- c(.CHAR_FIELDS, .INT_FIELDS)
    class <- .fieldToClass(field)
    for (i in seq_along(field)) {
        setClass(class[[i]], contains="AnnotationFilter", where=topenv())
        assign(class[[i]],
               .filterFactory(
                   field[[i]], class[[i]], field[[i]] %in% .CHAR_FIELDS
               ),
               envir=topenv())
    }
})

#' @aliases SymbolFilter SymbolFilter-class
#'
#' @title Filters for annotation objects
#'
#' These functions are used to create filters for annotation
#' resources.
#'
#' \code{supportedFilters()} lists all defined filters; filters are
#' only available for tables containing the field on which the filter
#' acts.
#'
#' @rdname AnnotationFilter
#'
#' @examples
#' supportedFilters()
#'
#' ## Create a SymbolFilter to filter on a gene's symbol.
#' sf <- SymbolFilter("BCL2")
#'
#' @export SymbolFilter
#' @exportClass SymbolFilter
#' @export
supportedFilters <- function() {
    .fieldToClass(c(.CHAR_FIELDS, .INT_FIELDS))
}

#' @param object An \code{AnnotationFilter} object
#'
#' @importFrom methods show
#' @rdname AnnotationFilter
#' @export
setMethod("show", "AnnotationFilter", function(object){
    cat("class:", class(object),
        "\ncondition:", .condition(object),
        "\nvalue:", .value(object), "\n")
})

############################################################
## Methods for the filter classes
## 

#' @aliases condition
#' @description \code{condition,condition<-} get or set the \code{condition}
#' value for the filter \code{object}.
#' 
#' @rdname AnnotationFilter
#' @export
setMethod("condition", "AnnotationFilter", function(object) {
    .condition(object)
})
#' @aliases condition<-
#'
#' @param value The value for the object.
#' 
#' @importFrom methods validObject
#' @rdname AnnotationFilter
#' @export
setReplaceMethod("condition", "AnnotationFilter", function(object, value) {
    object@condition <- as.character(value)
    if (validObject(object))
        object
})

#' @aliases value
#' @description \code{value,value<-} get or set the \code{value} for the filter
#' \code{object}.
#' 
#' @rdname AnnotationFilter
#' @export
setMethod("value", "AnnotationFilter", function(object) {
    .value(object)
})
#' @aliases value<-
#' @rdname AnnotationFilter
#' @export
setReplaceMethod("value", "AnnotationFilter", function(object, value) {
    if (.isCharacter(object))
        value <- as.character(value)
    object@value <- value
    if (validObject(object))
        object
})
