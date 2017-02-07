.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CHAR_FIELDS <- character()

.INT_FIELDS <- character()

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
    if (length(condition) != 1L)
        txt <- c(txt, "'condition' must be length 1")
    if (!condition %in% .OPS)
        txt <- c(txt,
                 sprintf("'condition' must be one of %s",
                         paste("'", .OPS, "'", collapse=", ")))
    if (isCharacter && !is.character(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take character value"))
    if (!isCharacter && (!is.integer(value)) || is.na(value))
        txt <- c(txt,
                 paste0("'", class(object),
                        "' can only take integer value"))
    if (condition  %in% c("startsWith", "endsWith", ">", "<", ">=", "<=") &&
        length(value) > 1L)
        txt <- c(txt,
                 paste0("'value' must be length 1 when condition is '",
                        condition, "'"))
    if (condition  %in% c("startsWith", "endsWith") && !isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with character value"))
    if (condition  %in% c(">", "<", ">=", "<=") && isCharacter)
        txt <- c(txt,
                 paste0("'", condition,
                        "' can only work with integer value"))
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
        new(class, field=field, condition=condition,
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

#' Filters for annotation objects
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
