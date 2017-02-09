.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CHAR_FIELDS <- c("exon_id", "exon_name", "gene_id", "genename", "gene_biotype",
                  "entrez", "symbol", "tx_id", "tx_name", "tx_biotype",
                  "protein_id", "uniprot", "seq_name", "seq_strand")

.INT_FIELDS <- c("cds_start", "cds_end", "exon_start", "exon_rank", "exon_end",
                 "gene_start", "gene_end", "tx_start", "tx_end")

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

#' @aliases CdsStartFilter CdsEndFilter ExonIdFilter ExonNameFilter
#' ExonStartFilter ExonEndFilter ExonRankFilter GeneIdFilter GenenameFilter
#' GeneBiotypeFilter GeneStartFilter GeneEndFilter EntrezFilter SymbolFilter
#' TxIdFilter TxNameFilter TxBiotypeFilter TxStartFilter TxEndFilter
#' ProteinIdFilter UniprotFilter SeqNameFilter SeqStrandFilter
#' AnnotationFilter-class CdsStartFilter-class CdsEndFilter-class
#' ExonIdFilter-class ExonNameFilter-class ExonStartFilter-class
#' ExonEndFilter-class ExonRankFilter-class GeneIdFilter-class
#' GenenameFilter-class GeneBiotypeFilter-class GeneStartFilter-class
#' GeneEndFilter-class EntrezFilter-class SymbolFilter-class TxIdFilter-class
#' TxNameFilter-class TxBiotypeFilter-class TxStartFilter-class TxEndFilter-class
#' ProteinIdFilter-class UniprotFilter-class SeqNameFilter-class
#' SeqStrandFilter-class
#'
#' 
#' @title Filters for annotation objects
#'
#' @description These functions are used to create filters for annotation
#' resources.
#'
#' \code{supportedFilters()} lists all defined filters; by default filters are
#' only available for tables containing the field on which the filter
#' acts. See the vignette for a description to use filters for databases in which
#' the database table column name differs from the default `field` of the filter.
#'
#' @usage CdsStartFilter(value, condition = "==")
#' CdsEndFilter(value, condition = "==")
#' ExonIdFilter(value, condition = "==")
#' ExonNameFilter(value, condition = "==")
#' ExonRankFilter(value, condition = "==")
#' ExonStartFilter(value, condition = "==")
#' ExonEndFilter(value, condition = "==")
#' GeneIdFilter(value, condition = "==")
#' GenenameFilter(value, condition = "==")
#' GeneBiotypeFilter(value, condition = "==")
#' GeneStartFilter(value, condition = "==")
#' GeneEndFilter(value, condition = "==")
#' EntrezFilter(value, condition = "==")
#' SymbolFilter(value, condition = "==")
#' TxIdFilter(value, condition = "==")
#' TxNameFilter(value, condition = "==")
#' TxBiotypeFilter(value, condition = "==")
#' TxStartFilter(value, condition = "==")
#' TxEndFilter(value, condition = "==")
#' ProteinIdFilter(value, condition = "==")
#' UniprotFilter(value, condition = "==")
#' SeqNameFilter(value, condition = "==")
#' SeqStrandFilter(value, condition = "==")
#'
#' @param value Value for the filter.
#'
#' @param condition character(1) defining the condition to be used in the filter.
#' For numeric/integer filters one of \code{"=="}, \code{"!="}, \code{">"},
#' \code{"<"}, \code{">="} and \code{"<="}. For character filter/values
#' \code{"=="}, \code{"!="}, \code{"startsWith"} and \code{"endsWith"} are
#' allowed. Default condition is "==".
#' 
#' @rdname AnnotationFilter
#'
#' @examples
#' supportedFilters()
#'
#' ## Create a SymbolFilter to filter on a gene's symbol.
#' sf <- SymbolFilter("BCL2")
#' sf
#'
#' ## Create a GeneStartFilter to filter based on the genes' chromosomal start
#' ## coordinates
#' gsf <- GeneStartFilter(10000, condition = ">")
#' gsf
#' 
#' @export CdsStartFilter CdsEndFilter ExonIdFilter ExonNameFilter
#' @export ExonStartFilter ExonEndFilter ExonRankFilter GeneIdFilter
#' @export GenenameFilter GeneBiotypeFilter GeneStartFilter GeneEndFilter
#' @export EntrezFilter SymbolFilter TxIdFilter TxNameFilter TxBiotypeFilter
#' @export TxStartFilter TxEndFilter ProteinIdFilter UniprotFilter SeqNameFilter
#' @export SeqStrandFilter
#' 
#' @exportClass AnnotationFilter CdsStartFilter CdsEndFilter ExonIdFilter
#' ExonNameFilter ExonStartFilter ExonEndFilter ExonRankFilter GeneIdFilter
#' GenenameFilter GeneBiotypeFilter GeneStartFilter GeneEndFilter EntrezFilter
#' SymbolFilter TxIdFilter TxNameFilter TxBiotypeFilter TxStartFilter TxEndFilter
#' ProteinIdFilter UniprotFilter SeqNameFilter SeqStrandFilter
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

############################################################
## Methods for the filter classes
## 

#' @aliases condition
#' @description \code{condition} get the \code{condition} value for the filter
#' \code{object}.
#' 
#' @rdname AnnotationFilter
#' @export
setMethod("condition", "AnnotationFilter", function(object) {
    .condition(object)
})

#' @aliases value
#' @description \code{value} get or set the \code{value} for the filter
#' \code{object}.
#' 
#' @rdname AnnotationFilter
#' @export
setMethod("value", "AnnotationFilter", function(object) {
    .value(object)
})

