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
    value(fl) <- "BCL2L11"
    expect_equal(value(fl), "BCL2L11")
    fl <- SymbolFilter(c(4, 5))
    expect_equal(value(fl), c("4", "5"))
    value(fl) <- 3
    expect_equal(value(fl), "3")
    expect_error(value(fl) <- NA)
    ## condition.
    expect_equal(condition(fl), "==")
    condition(fl) <- "!="
    expect_equal(condition(fl), "!=")
    expect_error(condition(fl) <- "<")
    expect_error(condition(fl) <- "")
    expect_error(condition(fl) <- c("==", ">"))
    expect_error(condition(fl) <- NULL)
    expect_error(condition(fl) <- NA)
    expect_error(condition(fl) <- 4)
})
