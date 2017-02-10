context("AnnotationFilter")

test_that("supportedFilters() works", {
    expect_true(inherits(supportedFilters(), "character"))
    expect_identical(
        length(supportedFilters()),
        length(c(
            AnnotationFilters:::.CHAR_FIELDS,
            AnnotationFilters:::.INT_FIELDS
        ))
    )
})

test_that("SymbolFilter as representative for character filters", {
    expect_true(validObject(new("SymbolFilter")))
    expect_error(SymbolFilter())
    expect_error(SymbolFilter(1, ">"))
    expect_error(SymbolFilter(1, "foo"))
    expect_error(SymbolFilter(c("foo","bar"), "startsWith"))
    ## Getter / setter
    fl <- SymbolFilter("BCL2")
    expect_equal(value(fl), "BCL2")
    fl <- SymbolFilter(c(4, 5))
    expect_equal(value(fl), c("4", "5"))
    fl <- SymbolFilter(3)
    expect_equal(value(fl), "3")
    expect_error(SymbolFilter(NA))
    ## condition.
    expect_equal(condition(fl), "==")
    fl <- SymbolFilter("a", condition = "!=")
    expect_equal(condition(fl), "!=")
    expect_error(SymbolFilter("a", condition = "<"))
    expect_error(SymbolFilter("a", condition = ""))
    expect_error(SymbolFilter("a", condition = c("==", ">")))
    expect_error(SymbolFilter("a", condition = NULL))
    expect_error(SymbolFilter("a", condition = NA))
    expect_error(SymbolFilter("a", condition = 4))
})

test_that("GeneStartFilter as representative for integer filters", {
    gsf <- GeneStartFilter(10000, condition = ">")
    expect_equal(condition(gsf), ">")
    expect_error(GeneStartFilter("3"))
    expect_error(GeneStartFilter("B"))
    expect_error(GeneStartFilter(NA))
    expect_error(GeneStartFilter(NULL))
    expect_error(GeneStartFilter())
    ## Condition
    expect_error(GeneStartFilter(10000, condition = "startsWith"))
    expect_error(GeneStartFilter(10000, condition = "endsWith"))
    expect_error(GeneStartFilter(10000, condition = c("==", "<")))
})

test_that("GRangesFilter", {
    grf <- GRangesFilter(as("chr10:87869000-87876000", "GRanges"))
    expect_equal(condition(grf), "overlapping")
    expect_error(GRangesFilter(value = 3))
    expect_error(GRangesFilter(as("chr10:87869000-87876000", "GRanges"),
                               condition = "=="))
    grf <- GRangesFilter(as("chr10:87869000-87876000", "GRanges"),
                         condition = "within", feature = "tx")
    expect_equal(condition(grf), "within")
    expect_equal(grf@feature, "tx")
})
