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
