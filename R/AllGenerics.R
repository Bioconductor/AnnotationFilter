## Generic methods.
setGeneric("condition", function(object, ...) standardGeneric("condition"))

setGeneric("field", function(object, ...) standardGeneric("field"))

setGeneric("value", function(object, ...) standardGeneric("value"))

setGeneric("logicOp", function(object, ...) standardGeneric("logicOp"))

setGeneric("not", function(object, ...) standardGeneric("not"))

setGeneric("simplify", function(object, ...) standardGeneric("simplify"))

setGeneric("convertFilter", function(object, ...)
    standardGeneric("convertFilter"))

setGeneric("distributeNegation", function(object, ...)
    standardGeneric("distributeNegation"))

setGeneric("supportedFilters", function(object, ...)
    standardGeneric("supportedFilters"))
