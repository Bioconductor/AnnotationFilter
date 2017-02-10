#' @name AnnotationFilter
#'
#' @title Filters for annotation objects
#'
#' @aliases CdsStartFilter CdsEndFilter ExonIdFilter ExonNameFilter
#'     ExonStartFilter ExonEndFilter ExonRankFilter GeneIdFilter
#'     GenenameFilter GeneBiotypeFilter GeneStartFilter GeneEndFilter
#'     EntrezFilter SymbolFilter TxIdFilter TxNameFilter
#'     TxBiotypeFilter TxStartFilter TxEndFilter ProteinIdFilter
#'     UniprotFilter SeqNameFilter SeqStrandFilter
#'     AnnotationFilter-class CharacterFilter-class
#'     IntegerFilter-class CdsStartFilter-class CdsEndFilter-class
#'     ExonIdFilter-class ExonNameFilter-class ExonStartFilter-class
#'     ExonEndFilter-class ExonRankFilter-class GeneIdFilter-class
#'     GenenameFilter-class GeneBiotypeFilter-class
#'     GeneStartFilter-class GeneEndFilter-class EntrezFilter-class
#'     SymbolFilter-class TxIdFilter-class TxNameFilter-class
#'     TxBiotypeFilter-class TxStartFilter-class TxEndFilter-class
#'     ProteinIdFilter-class UniprotFilter-class SeqNameFilter-class
#'     SeqStrandFilter-class supportedFilters
#'     show,AnnotationFilter-method show,CharacterFilter-method
#'     show,IntegerFilter-method show,GRangesFilter-method
#'
#' @description
#'
#' The filters extending the base \code{AnnotationFilter} class
#' represent a simple filtering concept for annotation resources.
#' Each filter object is thought to filter on a single (database)
#' table column based on the provided values and the defined
#' condition.
#'
#' Filter instances should be created using the constructor functions (e.g.
#' \code{GeneIdFilter}) and not by calls to \code{new}.
#'
#' \code{supportedFilters()} lists all defined filters. Packages using
#' \code{AnnotationFilters} should implement the \code{supportedFilters} for
#' their annotation resource object (e.g. for \code{object = "EnsDb"} in the
#' \code{ensembldb} package) to list all supported filters for the specific
#' resource.
#'
#' @details
#'
#' By default filters are only available for tables containing the
#' field on which the filter acts (i.e. that contain a column with the
#' name matching the value of the \code{field} slot of the
#' object). See the vignette for a description to use filters for
#' databases in which the database table column name differs from the
#' default \code{field} of the filter.
#'
#' @usage
#'
#' CdsStartFilter(value, condition = "==")
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
#' @param value Value for the filter. For \code{GRangesFilter} a
#'     \code{\link[GenomicRanges]{GRanges}} object.
#'
#' @param condition character(1) defining the condition to be used in
#'     the filter.  For numeric/integer filters one of \code{"=="},
#'     \code{"!="}, \code{">"}, \code{"<"}, \code{">="} and
#'     \code{"<="}. For character filter/values \code{"=="},
#'     \code{"!="}, \code{"startsWith"} and \code{"endsWith"} are
#'     allowed. Default condition is "==". For \code{GRangesFilter} it
#'     can be \code{"within"} (for the feature to be completely within
#'     the range) or \code{"overlapping"}, for the feature to be
#'     (partially) overlapping with the range.
NULL

.OPS <- c("==", "!=", "startsWith", "endsWith", ">", "<", ">=", "<=")

.CHAR_FIELDS <- c("exon_id", "exon_name", "gene_id", "genename", "gene_biotype",
                  "entrez", "symbol", "tx_id", "tx_name", "tx_biotype",
                  "protein_id", "uniprot", "seq_name", "seq_strand")

.INT_FIELDS <- c("cds_start", "cds_end", "exon_start", "exon_rank", "exon_end",
                 "gene_start", "gene_end", "tx_start", "tx_end")

############################################################
## AnnotationFilter
##

#' @exportClass AnnotationFilter
.AnnotationFilter <- setClass(
    "AnnotationFilter",
    contains = "VIRTUAL",
    slots = c(
        field="character",
        condition="character",
        value="ANY"
    ),
    prototype=list(
        condition= "=="
    )
)

setValidity("AnnotationFilter", function(object) {
    txt <- character()

    value <- .value(object)
    condition <- .condition(object)
    isCharacter <- is(value, "character")
    condNa <- any(is.na(condition))
    condOne <- length(condition) == 1L

    if (!condOne)
        txt <- c(txt, "'condition' must be length 1")
    if (condNa)
        txt <- c(txt, "'condition' can not be NA")

    test0 <- condition  %in% c("startsWith", "endsWith", ">", "<", ">=", "<=")
    if (!condNa && condOne && test0 && length(value) > 1L)
        txt <- c(txt, paste0("'", condition, "' requires length 1 'value'"))

    if (any(is.na(value)))
        txt <- c(txt, "'value' can not be NA")

    if (length(txt)) txt else TRUE
})

## slot accessors

.field <- function(object) object@field

.condition <- function(object) object@condition

.value <- function(object) object@value

#' @aliases condition
#' @description \code{condition()} get the \code{condition} value for
#'     the filter \code{object}.
#'
#' @param object An \code{AnnotationFilter} object.
#' @rdname AnnotationFilter
#' @export
setMethod("condition", "AnnotationFilter", .condition)

#' @aliases value
#' @description \code{value()} get the \code{value} for the filter
#'     \code{object}.
#'
#' @rdname AnnotationFilter
#' @export
setMethod("value", "AnnotationFilter", .value)

#' @importFrom methods show
#' @export
setMethod("show", "AnnotationFilter", function(object){
    cat("class:", class(object),
        "\ncondition:", .condition(object), "\n")
})

############################################################
## CharacterFilter, IntegerFilter
##

#' @exportClass CharacterFilter
.CharacterFilter <- setClass(
    "CharacterFilter",
    contains = c("VIRTUAL", "AnnotationFilter"),
    slots = c(value = "character"),
    prototype = list(
        value = character()
    )
)

setValidity("CharacterFilter", function(object) {
    txt <- character()

    condition <- .condition(object)
    test0 <-length(condition) == 1L && !is.na(condition)

    test <- test0 && condition %in% c("==", "!=", "startsWith", "endsWith")
    if (!test)
        txt <- c(txt, paste("'", condition, "' not allowed"))
})

#' @importFrom methods show callNextMethod
#' @export
setMethod("show", "CharacterFilter", function(object) {
    callNextMethod()
    cat("value:", .value(object), "\n")
})

#' @exportClass IntegerFilter
.IntegerFilter <- setClass(
    "IntegerFilter",
    contains = c("VIRTUAL", "AnnotationFilter"),
    slots = c(value = "integer"),
    prototype = list(
        value = integer()
    )
)

setValidity("IntegerFilter", function(object) {
    txt <- character()

    condition <- .condition(object)
    test0 <-length(condition) == 1L && !is.na(condition)

    test <-  test0 && condition %in% setdiff(.OPS, c("startsWith", "endsWith"))
    if (!test)
        txt <- c(txt, paste("'", condition, "' not allowed"))
})

#' @export
setMethod("show", "IntegerFilter", function(object) {
    callNextMethod()
    cat("value:", .value(object), "\n")
})

#' @rdname AnnotationFilter
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @exportClass GRangesFilter
.GRangesFilter <- setClass(
    "GRangesFilter",
    contains = "AnnotationFilter",
    slots = c(
        value = "GRanges",
        feature = "character"
    ),
    prototype = list(
        value  = GRanges(),
        condition = "overlapping",
        field = "granges",
        feature = "gene"
    )
)

#' @rdname AnnotationFilter
#' @examples
#' ## filter by GRanges
#' GRangesFilter(GenomicRanges::GRanges("chr10:87869000-87876000"))
#' @export
GRangesFilter <-
    function(value, condition = "overlapping", feature = "gene")
{
    .GRangesFilter(
        field = "granges",
        value = value,
        condition = condition,
        feature = feature)
}

.feature <- function(object) object@feature

setValidity("GRangesFilter", function(object) {
    condition <- .condition(object)
    txt <- character()
    if (!(condition %in% c("within", "overlapping")))
        txt <- c(txt, "'condition' must be \"within\" or \"overlapping\"")
    if (length(txt)) txt else TRUE
})

#' @importFrom GenomicRanges show
#' @export
setMethod("show", "GRangesFilter", function(object) {
    callNextMethod()
    cat("feature:", .feature(object),
        "\nvalue:\n")
    show(value(object))
})


############################################################
## Create install-time classes
##

#' @rdname AnnotationFilter
#' @name AnnotationFilter
#'
#' @param feature \code{character(1)} defining on what feature the
#'     \code{GRangesFilter} should be applied. Choices could be
#'     \code{"gene"}, \code{"tx"} or \code{"exon"}.
#'
#' @examples
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
#' @export GenenameFilter GeneBiotypeFilter GeneStartFilter
#' @export GeneEndFilter EntrezFilter SymbolFilter TxIdFilter
#' @export TxNameFilter TxBiotypeFilter TxStartFilter TxEndFilter
#' @export ProteinIdFilter UniprotFilter SeqNameFilter SeqStrandFilter
#'
#' @importFrom methods new
#'
#' @exportClass CdsStartFilter CdsEndFilter ExonIdFilter
#'     ExonNameFilter ExonStartFilter ExonEndFilter ExonRankFilter
#'     GeneIdFilter GenenameFilter GeneBiotypeFilter GeneStartFilter
#'     GeneEndFilter EntrezFilter SymbolFilter TxIdFilter TxNameFilter
#'     TxBiotypeFilter TxStartFilter TxEndFilter ProteinIdFilter
#'     UniprotFilter SeqNameFilter SeqStrandFilter
NULL

.fieldToClass <- function(field) {
    class <- sub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, if (length(class)) "Filter" else character(0))
}

.filterFactory <- function(field, class) {
    force(field); force(class)          # watch for lazy evaluation
    as.value <-
        if (field %in% .CHAR_FIELDS) {
            as.character
        } else {
            function(x) {
                stopifnot(is.numeric(x))
                as.integer(x)
            }
        }
    function(value, condition = "==") {
        value <- as.value(value)
        condition <- as.character(condition)
        new(class, field=field, condition = condition, value=value)
    }
}

local({
    makeClass <- function(field, contains) {
        class <- .fieldToClass(field)
        for (i in seq_along(field)) {
            setClass(class[[i]], contains=contains, where=topenv())
            assign(
                class[[i]],
                .filterFactory(field[[i]], class[[i]]),
                envir=topenv()
            )
        }
    }
    makeClass(.CHAR_FIELDS, "CharacterFilter")
    makeClass(.INT_FIELDS, "IntegerFilter")
})

############################################################
## Utilities - supportedFilters
##

.supportedFilters <- function() {
    sort(c(.fieldToClass(c(.CHAR_FIELDS, .INT_FIELDS)), "GRangesFilter"))
}

#' @rdname AnnotationFilter
#' @examples
#' supportedFilters()
#' @export
setMethod("supportedFilters", "missing", function(object) {
    .supportedFilters()
})
