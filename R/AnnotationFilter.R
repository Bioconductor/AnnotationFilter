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
#' table column using the provided values and the defined condition.
#'
#' Filter instances created using the constructor functions (e.g.
#' \code{GeneIdFilter}).
#'
#' \code{supportedFilters()} lists all defined filters. It returns a two column
#' \code{data.frame} with the filter class name and its default field.
#' Packages using \code{AnnotationFilter} should implement the
#' \code{supportedFilters} for their annotation resource object (e.g. for
#' \code{object = "EnsDb"} in the \code{ensembldb} package) to list all
#' supported filters for the specific resource.
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
#' CdsStartFilter(value, condition = "==", not = FALSE)
#' CdsEndFilter(value, condition = "==", not = FALSE)
#' ExonIdFilter(value, condition = "==", not = FALSE)
#' ExonNameFilter(value, condition = "==", not = FALSE)
#' ExonRankFilter(value, condition = "==", not = FALSE)
#' ExonStartFilter(value, condition = "==", not = FALSE)
#' ExonEndFilter(value, condition = "==", not = FALSE)
#' GeneIdFilter(value, condition = "==", not = FALSE)
#' GenenameFilter(value, condition = "==", not = FALSE)
#' GeneBiotypeFilter(value, condition = "==", not = FALSE)
#' GeneStartFilter(value, condition = "==", not = FALSE)
#' GeneEndFilter(value, condition = "==", not = FALSE)
#' EntrezFilter(value, condition = "==", not = FALSE)
#' SymbolFilter(value, condition = "==", not = FALSE)
#' TxIdFilter(value, condition = "==", not = FALSE)
#' TxNameFilter(value, condition = "==", not = FALSE)
#' TxBiotypeFilter(value, condition = "==", not = FALSE)
#' TxStartFilter(value, condition = "==", not = FALSE)
#' TxEndFilter(value, condition = "==", not = FALSE)
#' ProteinIdFilter(value, condition = "==", not = FALSE)
#' UniprotFilter(value, condition = "==", not = FALSE)
#' SeqNameFilter(value, condition = "==", not = FALSE)
#' SeqStrandFilter(value, condition = "==", not = FALSE)
#'
#' @param value \code{character()}, \code{integer()}, or
#'     \code{GRanges()} value for the filter
#'
#' @param condition \code{character(1)} defining the condition to be
#'     used in the filter. For \code{IntegerFilter}, one of
#'     \code{"=="}, \code{"!="}, \code{">"}, \code{"<"}, \code{">="}
#'     or \code{"<="}. For \code{CharacterFilter}, one of \code{"=="},
#'     \code{"!="}, \code{"startsWith"}, \code{"endsWith"} or \code{"contains"}.
#'     Default condition is \code{"=="}.
#'
#' @param not \code{logical(1)} whether the \code{AnnotationFilter} is negated.
#'     \code{TRUE} indicates is negated (!). \code{FALSE} indicates not
#'     negated. Default not is \code{FALSE}.
#'
#' @return The constructor function return an object extending
#'     \code{AnnotationFilter}. For the return value of the other methods see
#'     the methods' descriptions.
#' 
#' @seealso \code{\link{AnnotationFilterList}} for combining
#'     \code{AnnotationFilter} objects.
NULL

.CONDITION <- list(
    IntegerFilter = c("==", "!=", ">", "<", ">=", "<="),
    CharacterFilter =  c("==", "!=", "startsWith", "endsWith", "contains"),
    GRangesFilter = c("any", "start", "end", "within", "equal")
)

.FIELD <- list(
    CharacterFilter = c(
        "exon_id", "exon_name", "gene_id", "genename", "gene_biotype",
        "entrez", "symbol", "tx_name", "tx_biotype",
        "protein_id", "uniprot", "seq_name", "seq_strand"),
    IntegerFilter = c(
        "cds_start", "cds_end", "exon_start", "exon_rank", "exon_end",
        "gene_start", "gene_end", "tx_id", "tx_start", "tx_end")
)

.valid_condition <- function(condition, class) {
    txt <- character()

    test0 <- length(condition) == 1L
    if (!test0)
        txt <- c(txt, "'condition' must be length 1")

    test1 <- test0 && (condition %in% .CONDITION[[class]])
    if (!test1) {
        value <- paste(sQuote(.CONDITION[[class]]), collapse=" ")
        txt <- c(txt, paste0("'", condition, "' must be in ", value))
    }

    if (length(txt)) txt else TRUE
}

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
        value="ANY",
        not="logical"
    ),
    prototype=list(
        condition= "==",
        not= FALSE
    )
)

setValidity("AnnotationFilter", function(object) {
    txt <- character()

    value <- .value(object)
    condition <- .condition(object)
    not <- .not(object)
    test_len <- length(condition) == 1L
    test_NA <- !any(is.na(condition))

    if (test_len && !test_NA)
        txt <- c(txt, "'condition' can not be NA")
    test0 <- test_len && test_NA

    test1 <- condition  %in% c("startsWith", "endsWith", "contains", ">",
                               "<", ">=", "<=")
    if (test0 && test1 && length(value) > 1L)
        txt <- c(txt, paste0("'", condition, "' requires length 1 'value'"))

    if(length(not) != 1)
        txt <- c(txt, '"not" value must be of length 1.')

    if (any(is.na(value)))
        txt <- c(txt, "'value' can not be NA")

    if (length(txt)) txt else TRUE
})

.field <- function(object) object@field

.condition <- function(object) object@condition

.value <- function(object) object@value

.not <- function(object) object@not

#' @rdname AnnotationFilter
#'
#' @aliases condition
#'
#' @description \code{condition()} get the \code{condition} value for
#'     the filter \code{object}.
#'
#' @param object An \code{AnnotationFilter} object.
#' 
#' @export
setMethod("condition", "AnnotationFilter", .condition)

#' @rdname AnnotationFilter
#'
#' @aliases value
#'
#' @description \code{value()} get the \code{value} for the filter
#'     \code{object}.
#'
#' @export
setMethod("value", "AnnotationFilter", .value)

#' @rdname AnnotationFilter
#'
#' @aliases field
#'
#' @description \code{field()} get the \code{field} for the filter
#'     \code{object}.
#'
#' @export
setMethod("field", "AnnotationFilter", .field)

#' @rdname AnnotationFilter
#'
#' @description \code{not()} get the \code{not} for the filter \code{object}.
#'
#' @export
setMethod("not", "AnnotationFilter", .not)

#' @importFrom methods show
#'
#' @export
setMethod("show", "AnnotationFilter", function(object){
    if(.not(object)) cat("NOT\n")
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
    .valid_condition(.condition(object), "CharacterFilter")
})

#' @importFrom methods show callNextMethod
#'
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
    .valid_condition(.condition(object), "IntegerFilter")
})

#' @export
setMethod("show", "IntegerFilter", function(object) {
    callNextMethod()
    cat("value:", .value(object), "\n")
})

#' @rdname AnnotationFilter
#'
#' @importFrom GenomicRanges GRanges
#'
#' @importClassesFrom GenomicRanges GRanges
#'
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
        condition = "any",
        field = "granges",
        feature = "gene"
    )
)

setValidity("GRangesFilter", function(object) {
    .valid_condition(.condition(object), "GRangesFilter")
})

.feature <- function(object) object@feature

#' @rdname AnnotationFilter
#'
#' @param type \code{character(1)} indicating how overlaps are to be
#'     filtered. See \code{findOverlaps} in the IRanges package for a
#'     description of this argument.
#'
#' @examples
#' ## filter by GRanges
#' GRangesFilter(GenomicRanges::GRanges("chr10:87869000-87876000"))
#' @export
GRangesFilter <-
    function(value, feature = "gene",
             type = c("any", "start", "end", "within", "equal"))
{
    condition <- match.arg(type)
    .GRangesFilter(
        field = "granges",
        value = value,
        condition = condition,
        feature = feature)
}

.feature <- function(object) object@feature

#' @aliases feature
#'
#' @description \code{feature()} get the \code{feature} for the
#'     \code{GRangesFilter} \code{object}.
#'
#' @rdname AnnotationFilter
#'
#' @export
feature <- .feature

#' @importFrom GenomicRanges show
#'
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
#'
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
    class <- gsub("_([[:alpha:]])", "\\U\\1", field, perl=TRUE)
    class <- sub("^([[:alpha:]])", "\\U\\1", class, perl=TRUE)
    paste0(class, if (length(class)) "Filter" else character(0))
}

.filterFactory <- function(field, class) {
    force(field); force(class)          # watch for lazy evaluation
    as.value <-
        if (field %in% .FIELD[["CharacterFilter"]]) {
            function(x) {
#               if(!is.character(x))
#                  stop("Input to a ", field,
#                       "filter must be a character vector.")
                as.character(x)
            }
        } else {
            function(x) {
                if(!is.numeric(x))
                    stop("Input to a ", field,
                         "filter must be a numeric vector.")
                as.integer(x)
            }
        }

    function(value, condition = "==", not = FALSE) {
        value <- as.value(value)
        condition <- as.character(condition)
        not <- as.logical(not)
        new(class, field=field, condition = condition, value=value, not=not)
    }
}

local({
    makeClass <- function(contains) {
        fields <- .FIELD[[contains]]
        classes <- .fieldToClass(fields)
        for (i in seq_along(fields)) {
            setClass(classes[[i]], contains=contains, where=topenv())
            assign(
                classes[[i]],
                .filterFactory(fields[[i]], classes[[i]]),
                envir=topenv()
            )
        }
    }
    for (contains in names(.FIELD))
        makeClass(contains)
})

############################################################
## Utilities 
##

.convertFilter <- function(object) {
    field <- field(object)
    if (field == "granges")
        stop("GRangesFilter cannot be converted using convertFilter().")
    value <- value(object)
    condition <- condition(object)
    not <- not(object)

    op <- switch(
        condition,
        "==" = if (length(value) == 1) "==" else "%in%",
        "!=" = if (length(value) == 1) "!=" else "%in%",
        "startsWith" = "%like%",
        "endsWith" = "%like%",
        "contains" = "%like%"
    )

    not_val <- ifelse(not, '!', '')

    if (condition %in% c("==", "!="))
        value <- paste0("'", value, "'", collapse=", ")

    if (!is.null(op) && op %in% c("==", "!="))
        sprintf("%s%s %s %s", not_val, field, op, value)
    else if ((condition == "==") && op == "%in%")
        sprintf("%s%s %s c(%s)", not_val, field, op, value)
    else if ((condition == "!=") && op == "%in%")
        if(not) sprintf("%s %s c(%s)", field, op, value)
        else sprintf("!%s%s %s c(%s)", not_val, field, op, value)
    else if (condition == "startsWith")
        sprintf("%s%s %s '%s%%'", not_val, field, op, value)
    else if (condition == "endsWith")
        sprintf("%s%s %s '%%%s'", not_val, field, op, value)
    else if (condition == "contains")
        sprintf("%s%s %s '%s'", not_val, field, op, value)
    else if (condition %in% c(">", "<", ">=", "<=")) {
        sprintf("%s%s %s %s", not_val, field, condition, as.integer(value))
    }
}

#' @rdname AnnotationFilter
#'
#' @description Converts an \code{AnnotationFilter} object to a 
#'      \code{character(1)} giving an equation that can be used as input to
#'      a \code{dplyr} filter.
#'
#' @return \code{character(1)} that can be used as input to a \code{dplyr} 
#'      fitlter.
#'
#' @examples
#' filter <- SymbolFilter("ADA", "==")
#' result <- convertFilter(filter)
#' result
#' @export
setMethod("convertFilter", "AnnotationFilter", .convertFilter)

.FILTERS_WO_FIELD <- c("GRangesFilter")

.supportedFilters <- function() {
    fields <- unlist(.FIELD, use.names=FALSE)
    filters <- .fieldToClass(fields)
    d <- data.frame(
      filter=c(filters, .FILTERS_WO_FIELD),
      field=c(fields, "granges") #rep(NA, length(.FILTERS_WO_FIELD)))
    )
    d[order(d$filter),]
}

#' @rdname AnnotationFilter
#'
#' @examples
#' supportedFilters()
#' @export
setMethod("supportedFilters", "missing", function(object) {
    .supportedFilters()
})
