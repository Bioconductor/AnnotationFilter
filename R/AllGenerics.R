## Generic methods.
setGeneric("condition", function(object, ...) standardGeneric("condition"))

setGeneric("field", function(object, ...) standardGeneric("field"))

setGeneric("value", function(object, ...) standardGeneric("value"))

setGeneric("logicOp", function(object, ...) standardGeneric("logicOp"))

setGeneric("not", function(object, ...) standardGeneric("not"))

setGeneric("supportedFilters", function(object, ...)
    standardGeneric("supportedFilters"))
